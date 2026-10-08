/*
 * hornet -- C implementation of the Hornet allele-specific mapping tools.
 *
 * Produces exactly the same output as the Python scripts in ../mapping
 * (find_intersecting_snps.py, filter_remapped_reads.py, rmdup.py,
 * rmdup_pe.py): the same BAM records in the same order, the same fastq
 * content and the same run statistics.  rmdup/rmdup_pe reproduce CPython's
 * Mersenne Twister, so --seed N gives the same choice of reads as the Python
 * scripts run with random.seed(N).
 *
 * Usage:
 *   hornet find_intersecting_snps [-p] [-P] [-C chrom] [-t threads] in.bam snp_dir
 *   hornet filter_remapped_reads  [-p] [-t threads] to.remap.bam remap.bam keep.bam
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdarg.h>
#include <time.h>
#include <getopt.h>
#include <dirent.h>
#include <errno.h>
#include <fnmatch.h>
#include <pthread.h>
#include <unistd.h>
#include <fcntl.h>
#include <zlib.h>

#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/thread_pool.h>
#include <htslib/khash.h>
#include <htslib/vcf.h>

/* Maximum number of allelic combinations (including the original) per read
 * or pair: --max-seqs, 32 in --legacy mode as in the original code. */
static int max_seqs = 64;
#define MAX_SEQS_PER_READ max_seqs

/* Reads dropped first (default mode): unmapped, secondary, QC-fail,
 * supplementary; for pairs also reads whose mate is unmapped. */
#define FLAG_DROP_SINGLE (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FSUPPLEMENTARY)
#define FLAG_DROP_PAIRED (FLAG_DROP_SINGLE | BAM_FMUNMAP)

/* --legacy: reproduce the original Hornet output exactly, including its bugs
 * (phased-mode statistics, qualities of reverse-strand reads in the remap
 * fastqs, single-end filter_remapped_reads skipping reads). */
static int legacy = 0;

/* ------------------------------------------------------------------------ */
/* Utilities                                                                */
/* ------------------------------------------------------------------------ */

static void die(const char *fmt, ...)
{
    va_list ap;
    fflush(stdout);
    fprintf(stderr, "hornet: ");
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    exit(1);
}

static void *xmalloc(size_t n)
{
    void *p = malloc(n ? n : 1);
    if (!p) die("out of memory");
    return p;
}

static void *xrealloc(void *p, size_t n)
{
    p = realloc(p, n ? n : 1);
    if (!p) die("out of memory");
    return p;
}

static char *xstrdup(const char *s)
{
    char *p = strdup(s);
    if (!p) die("out of memory");
    return p;
}

#define GROW(ptr, n, m) do {                                  \
        if ((n) >= (m)) {                                      \
            (m) = (m) ? (m) * 2 : 16;                          \
            while ((m) <= (n)) (m) *= 2;                       \
            (ptr) = xrealloc((ptr), sizeof(*(ptr)) * (m));     \
        }                                                      \
    } while (0)

/* time.strftime("%b %d ") + time.strftime("%I:%M:%S") */
static const char *now_str(void)
{
    static char buf[64];
    time_t t = time(NULL);
    struct tm tm;
    localtime_r(&t, &tm);
    strftime(buf, sizeof buf, "%b %d %I:%M:%S", &tm);
    return buf;
}

/* "%.2f" % ((num / den)*100), "nan" on division by zero */
static void fmt_pct(char *buf, double num, double den)
{
    if (den == 0) strcpy(buf, "nan");
    else sprintf(buf, "%.2f", (num / den) * 100);
}

/* Python's int() for an ASCII token: [+-]digits with single '_' separators */
static int parse_pyint(const char *s, int64_t *out)
{
    int neg = 0;
    uint64_t v = 0;
    if (*s == '+' || *s == '-') neg = (*s++ == '-');
    if (*s < '0' || *s > '9') return -1;
    while (*s) {
        if (*s == '_') {
            s++;
            if (*s < '0' || *s > '9') return -1;
            continue;
        }
        if (*s < '0' || *s > '9') return -1;
        if (v > (uint64_t)INT64_MAX / 10) return -1;
        v = v * 10 + (uint64_t)(*s - '0');
        if (v > (uint64_t)INT64_MAX) return -1;
        s++;
    }
    *out = neg ? -(int64_t)v : (int64_t)v;
    return 0;
}

/* Whitespace as understood by Python's str.split() for ASCII text */
static inline int py_isspace(unsigned char c)
{
    return c == ' ' || (c >= '\t' && c <= '\r') || (c >= 0x1c && c <= 0x1f);
}

/* ------------------------------------------------------------------------ */
/* Line reader for (multi-member) gzip files with universal newlines        */
/* ------------------------------------------------------------------------ */

typedef struct {
    gzFile fp;
    unsigned char *buf;
    int len, pos, eof;
    char *line;
    size_t lsize;
} linereader_t;

static int lr_open(linereader_t *lr, const char *fname)
{
    memset(lr, 0, sizeof *lr);
    lr->fp = gzopen(fname, "rb");
    if (!lr->fp) return -1;
    gzbuffer(lr->fp, 1 << 20);
    lr->buf = xmalloc(1 << 20);
    lr->lsize = 256;
    lr->line = xmalloc(lr->lsize);
    return 0;
}

static inline int lr_getc(linereader_t *lr)
{
    if (lr->pos >= lr->len) {
        if (lr->eof) return -1;
        lr->len = gzread(lr->fp, lr->buf, 1 << 20);
        lr->pos = 0;
        if (lr->len <= 0) {
            int err;
            const char *msg = gzerror(lr->fp, &err);
            if (lr->len < 0 || (err != Z_OK && err != Z_STREAM_END))
                die("error reading gzip file: %s", msg);
            lr->eof = 1;
            lr->len = 0;
            return -1;
        }
    }
    return lr->buf[lr->pos++];
}

/* Returns line length (without terminator) or -1 at EOF. '\n', '\r\n' and
 * '\r' all end a line, like Python's text mode. */
static int lr_getline(linereader_t *lr)
{
    size_t n = 0;
    int c;
    for (;;) {
        c = lr_getc(lr);
        if (c < 0) {
            if (n == 0) return -1;
            break;
        }
        if (c == '\n') break;
        if (c == '\r') {
            int d = lr_getc(lr);
            if (d >= 0 && d != '\n') lr->pos--;
            break;
        }
        if (n + 2 > lr->lsize) {
            lr->lsize *= 2;
            lr->line = xrealloc(lr->line, lr->lsize);
        }
        lr->line[n++] = (char)c;
    }
    lr->line[n] = 0;
    return (int)n;
}

static void lr_close(linereader_t *lr)
{
    gzclose(lr->fp);
    free(lr->buf);
    free(lr->line);
}

/* Split in place on Python whitespace; returns number of tokens */
static int split_ws(char *s, char **tok, int maxtok)
{
    int n = 0;
    for (;;) {
        while (*s && py_isspace((unsigned char)*s)) s++;
        if (!*s) break;
        if (n < maxtok) tok[n] = s;
        n++;
        while (*s && !py_isspace((unsigned char)*s)) s++;
        if (!*s) break;
        *s++ = 0;
    }
    return n;
}

/* ------------------------------------------------------------------------ */
/* SNP table                                                                */
/* ------------------------------------------------------------------------ */

KHASH_MAP_INIT_STR(s2i, int)

typedef struct {
    int64_t *pos;      /* sorted, unique 0-based positions */
    uint32_t *code;    /* allele code for each position */
    int64_t n;
} snplist_t;

typedef struct {
    char **name;
    snplist_t *snps;
    int n, m;
    khash_t(s2i) *idx;
    /* interned allele strings (ref + alt) */
    char **alleles;
    int *alen;
    uint8_t *is_indel;
    int na, ma;
    khash_t(s2i) *aidx;
} snptab_t;

/* Interning table used while parsing one file */
typedef struct {
    char **str;
    int n, m;
    khash_t(s2i) *h;
} strtab_t;

static int strtab_get(strtab_t *t, const char *s)
{
    int ret;
    khint_t k = kh_get(s2i, t->h, s);
    if (k != kh_end(t->h)) return kh_val(t->h, k);
    GROW(t->str, t->n, t->m);
    t->str[t->n] = xstrdup(s);
    k = kh_put(s2i, t->h, t->str[t->n], &ret);
    kh_val(t->h, k) = t->n;
    return t->n++;
}

typedef struct {
    int64_t pos;
    uint32_t code;   /* local allele code */
    uint32_t chrom;  /* local chrom code (consolidated file only) */
} snprec_t;

typedef struct {
    const char *fname;
    int has_chrom;
    const char *chrom_only;
    snprec_t *rec;
    int64_t n, m;
    strtab_t alleles, chroms;
} loadjob_t;

static void *load_snp_file(void *arg)
{
    loadjob_t *j = arg;
    linereader_t lr;
    char *tok[5];
    char allele[1024];
    int ntok = j->has_chrom ? 4 : 3;
    int64_t lineno = 0;
    j->alleles.h = kh_init(s2i);
    j->chroms.h = kh_init(s2i);
    if (lr_open(&lr, j->fname) < 0) die("cannot open %s: %s", j->fname, strerror(errno));
    while (lr_getline(&lr) >= 0) {
        int64_t pos;
        int nt;
        lineno++;
        nt = split_ws(lr.line, tok, 5);
        if (nt != ntok)
            die("%s line %lld: expected %d columns, found %d", j->fname,
                (long long)lineno, ntok, nt);
        char **t = tok + (j->has_chrom ? 1 : 0);
        if (j->has_chrom && j->chrom_only && strcmp(tok[0], j->chrom_only) != 0)
            continue;
        if (parse_pyint(t[0], &pos) < 0)
            die("%s line %lld: invalid position '%s'", j->fname, (long long)lineno, t[0]);
        size_t lr_ = strlen(t[1]), la = strlen(t[2]);
        if (lr_ + la + 1 > sizeof allele)
            die("%s line %lld: allele too long", j->fname, (long long)lineno);
        memcpy(allele, t[1], lr_);
        memcpy(allele + lr_, t[2], la + 1);
        GROW(j->rec, j->n, j->m);
        j->rec[j->n].pos = pos - 1;
        j->rec[j->n].code = (uint32_t)strtab_get(&j->alleles, allele);
        j->rec[j->n].chrom = j->has_chrom ? (uint32_t)strtab_get(&j->chroms, tok[0]) : 0;
        j->n++;
    }
    lr_close(&lr);
    return NULL;
}

static int snptab_allele(snptab_t *T, const char *s)
{
    int ret;
    khint_t k = kh_get(s2i, T->aidx, s);
    if (k != kh_end(T->aidx)) return kh_val(T->aidx, k);
    GROW(T->alleles, T->na, T->ma);
    T->alen = xrealloc(T->alen, sizeof(int) * T->ma);
    T->is_indel = xrealloc(T->is_indel, T->ma);
    T->alleles[T->na] = xstrdup(s);
    T->alen[T->na] = (int)strlen(s);
    /* ('-' in alleles) or (max(len(i) for i in alleles) > 1): the second
     * test iterates over characters and is never true. */
    if (legacy)
        T->is_indel[T->na] = strchr(s, '-') != NULL;
    else  /* like WASP: anything but single-base ref and alt */
        T->is_indel[T->na] = strchr(s, '-') != NULL || T->alen[T->na] != 2;
    k = kh_put(s2i, T->aidx, T->alleles[T->na], &ret);
    kh_val(T->aidx, k) = T->na;
    return T->na++;
}

static int snptab_chrom(snptab_t *T, const char *name, int create)
{
    int ret;
    khint_t k = kh_get(s2i, T->idx, name);
    if (k != kh_end(T->idx)) return kh_val(T->idx, k);
    if (!create) return -1;
    GROW(T->name, T->n, T->m);
    T->snps = xrealloc(T->snps, sizeof(snplist_t) * T->m);
    T->name[T->n] = xstrdup(name);
    memset(&T->snps[T->n], 0, sizeof(snplist_t));
    k = kh_put(s2i, T->idx, T->name[T->n], &ret);
    kh_val(T->idx, k) = T->n;
    return T->n++;
}

typedef struct { int64_t pos; int64_t order; uint32_t code; } sortrec_t;

static int cmp_sortrec(const void *a, const void *b)
{
    const sortrec_t *x = a, *y = b;
    if (x->pos != y->pos) return x->pos < y->pos ? -1 : 1;
    return x->order < y->order ? -1 : (x->order > y->order);
}

/* Sort by position, keeping the last occurrence of duplicated positions
 * (dict assignment semantics). */
static void finalize_chrom(snplist_t *L, sortrec_t *r, int64_t n)
{
    int64_t i, k = 0;
    int sorted = 1;
    for (i = 1; i < n; i++)
        if (r[i].pos < r[i-1].pos) { sorted = 0; break; }
    if (!sorted) qsort(r, n, sizeof *r, cmp_sortrec);
    L->pos = xmalloc(sizeof(int64_t) * n);
    L->code = xmalloc(sizeof(uint32_t) * n);
    for (i = 0; i < n; i++) {
        if (i + 1 < n && r[i+1].pos == r[i].pos) continue;
        L->pos[k] = r[i].pos;
        L->code[k] = r[i].code;
        k++;
    }
    L->n = k;
}

static int has_suffix(const char *s, const char *suf)
{
    size_t a = strlen(s), b = strlen(suf);
    return a >= b && strcmp(s + a - b, suf) == 0;
}

static char *path_join(const char *dir, const char *name)
{
    size_t a = strlen(dir);
    char *p = xmalloc(a + strlen(name) + 2);
    if (a == 0) strcpy(p, name);
    else if (dir[a-1] == '/') sprintf(p, "%s%s", dir, name);
    else sprintf(p, "%s/%s", dir, name);
    return p;
}

static void add_job_records(snptab_t *T, loadjob_t **jobs, int njobs, int chrom)
{
    /* All jobs contribute to the same chromosome (per-chrom files). */
    int64_t n = 0, i, k = 0;
    int j;
    for (j = 0; j < njobs; j++) n += jobs[j]->n;
    if (n == 0) return;
    sortrec_t *r = xmalloc(sizeof(sortrec_t) * n);
    for (j = 0; j < njobs; j++) {
        loadjob_t *J = jobs[j];
        int *map = xmalloc(sizeof(int) * (J->alleles.n + 1));
        for (i = 0; i < J->alleles.n; i++) map[i] = snptab_allele(T, J->alleles.str[i]);
        for (i = 0; i < J->n; i++, k++) {
            r[k].pos = J->rec[i].pos;
            r[k].order = k;
            r[k].code = (uint32_t)map[J->rec[i].code];
        }
        free(map);
    }
    if (chrom < 0) die("internal error");
    finalize_chrom(&T->snps[chrom], r, n);
    free(r);
}

static void free_job(loadjob_t *J)
{
    int i;
    for (i = 0; i < J->alleles.n; i++) free(J->alleles.str[i]);
    for (i = 0; i < J->chroms.n; i++) free(J->chroms.str[i]);
    free(J->alleles.str);
    free(J->chroms.str);
    kh_destroy(s2i, J->alleles.h);
    kh_destroy(s2i, J->chroms.h);
    free(J->rec);
}

static snptab_t *get_snps(const char *snpdir, const char *chrom_only, int nthreads)
{
    snptab_t *T = calloc(1, sizeof *T);
    T->idx = kh_init(s2i);
    T->aidx = kh_init(s2i);
    char *all = path_join(snpdir, "all.txt.gz");
    if (access(all, F_OK) == 0) {
        printf("%s .... Loading snps from consolidated file.\n", now_str());
        fflush(stdout);
        loadjob_t J;
        memset(&J, 0, sizeof J);
        J.fname = all;
        J.has_chrom = 1;
        J.chrom_only = chrom_only;
        load_snp_file(&J);
        /* group records by chromosome, preserving file order */
        int nc = J.chroms.n, c;
        int64_t i, *cnt = calloc(nc + 1, sizeof(int64_t));
        for (i = 0; i < J.n; i++) cnt[J.rec[i].chrom]++;
        sortrec_t **bucket = xmalloc(sizeof(sortrec_t *) * (nc + 1));
        int64_t *fill = calloc(nc + 1, sizeof(int64_t));
        for (c = 0; c < nc; c++) bucket[c] = xmalloc(sizeof(sortrec_t) * (cnt[c] + 1));
        int *map = xmalloc(sizeof(int) * (J.alleles.n + 1));
        for (i = 0; i < J.alleles.n; i++) map[i] = snptab_allele(T, J.alleles.str[i]);
        for (i = 0; i < J.n; i++) {
            c = J.rec[i].chrom;
            sortrec_t *r = &bucket[c][fill[c]];
            r->pos = J.rec[i].pos;
            r->order = fill[c]++;
            r->code = (uint32_t)map[J.rec[i].code];
        }
        for (c = 0; c < nc; c++) {
            int ci = snptab_chrom(T, J.chroms.str[c], 1);
            finalize_chrom(&T->snps[ci], bucket[c], cnt[c]);
            free(bucket[c]);
        }
        free(bucket); free(fill); free(cnt); free(map);
        free_job(&J);
        free(all);
        return T;
    }
    free(all);

    /* glob(path.join(snpdir, '*.txt.gz')) in directory order */
    DIR *d = opendir(snpdir);
    loadjob_t *jobs = NULL;
    char **jchrom = NULL;
    int nj = 0, mj = 0, i;
    if (d) {
        struct dirent *e;
        while ((e = readdir(d)) != NULL) {
            if (e->d_name[0] == '.') continue;
            if (!has_suffix(e->d_name, ".txt.gz")) continue;
            /* path.basename(fname).split('.snps.txt.gz')[0] */
            char *chrom = xstrdup(e->d_name);
            char *p = strstr(chrom, ".snps.txt.gz");
            if (p) *p = 0;
            if (chrom_only && strcmp(chrom, chrom_only) != 0) { free(chrom); continue; }
            if (nj >= mj) {
                mj = mj ? mj * 2 : 32;
                jobs = xrealloc(jobs, sizeof(loadjob_t) * mj);
                jchrom = xrealloc(jchrom, sizeof(char *) * mj);
            }
            memset(&jobs[nj], 0, sizeof(loadjob_t));
            jobs[nj].fname = path_join(snpdir, e->d_name);
            jchrom[nj] = chrom;
            nj++;
        }
        closedir(d);
    }
    for (i = 0; i < nj; i++)
        printf("%s .... Loading snps from %s\n", now_str(), jobs[i].fname);
    fflush(stdout);

    /* Parse files in parallel */
    int nt = nthreads < 1 ? 1 : nthreads;
    if (nt > nj) nt = nj;
    if (nt <= 1) {
        for (i = 0; i < nj; i++) load_snp_file(&jobs[i]);
    } else {
        pthread_t *th = xmalloc(sizeof(pthread_t) * nj);
        int next = 0, running = 0, done = 0;
        int *started = calloc(nj, sizeof(int));
        /* simple batches of nt threads */
        while (done < nj) {
            int batch_start = next;
            running = 0;
            while (next < nj && running < nt) {
                if (pthread_create(&th[next], NULL, load_snp_file, &jobs[next]))
                    die("pthread_create failed");
                started[next] = 1;
                next++;
                running++;
            }
            for (i = batch_start; i < next; i++) { pthread_join(th[i], NULL); done++; }
        }
        free(th);
        free(started);
    }

    /* Merge, in glob order, files mapping to the same chromosome */
    int *used = calloc(nj + 1, sizeof(int));
    loadjob_t **grp = xmalloc(sizeof(loadjob_t *) * (nj + 1));
    for (i = 0; i < nj; i++) {
        int k, ng = 0;
        if (used[i]) continue;
        for (k = i; k < nj; k++)
            if (!used[k] && strcmp(jchrom[k], jchrom[i]) == 0) { used[k] = 1; grp[ng++] = &jobs[k]; }
        int64_t tot = 0;
        for (k = 0; k < ng; k++) tot += grp[k]->n;
        if (tot > 0) {
            int ci = snptab_chrom(T, jchrom[i], 1);
            add_job_records(T, grp, ng, ci);
        }
    }
    for (i = 0; i < nj; i++) { free_job(&jobs[i]); free((char *)jobs[i].fname); free(jchrom[i]); }
    free(jobs); free(jchrom); free(used); free(grp);
    if (T->n == 0) {
        printf("ERROR: Could not find appropriate SNP files\n");
        fflush(stdout);
        exit(1);
    }
    return T;
}


/* Heterozygous sites of one sample from a VCF/BCF. SNVs store the sample's
 * two alleles; every reference base covered by a heterozygous indel or other
 * non-SNV allele is an indel site ("-"), which wins over SNVs at the same
 * base. Homozygous and missing genotypes are ignored. */
typedef struct { sortrec_t *r; int64_t n, m; } recvec_t;

static int is_acgt(const char *a)
{
    return a[0] && !a[1] && (a[0] == 'A' || a[0] == 'C' || a[0] == 'G' || a[0] == 'T');
}

static snptab_t *get_snps_vcf(const char *vcf, const char *sample, const char *chrom_only)
{
    htsFile *fp = bcf_open(vcf, "r");
    if (!fp) die("cannot open %s", vcf);
    bcf_hdr_t *h = bcf_hdr_read(fp);
    if (!h) die("cannot read header of %s", vcf);
    int ns = bcf_hdr_nsamples(h);
    if (!sample) {
        if (ns != 1) die("%s has %d samples; choose one with --sample", vcf, ns);
        sample = h->samples[0];
    } else if (bcf_hdr_id2int(h, BCF_DT_SAMPLE, sample) < 0) {
        die("sample %s not found in %s", sample, vcf);
    }
    char *smp = xstrdup(sample);
    if (bcf_hdr_set_samples(h, smp, 0) != 0) die("cannot select sample %s", smp);
    printf("%s .... Loading heterozygous sites of sample %s from %s\n", now_str(), smp, vcf);
    fflush(stdout);

    snptab_t *T = calloc(1, sizeof *T);
    T->idx = kh_init(s2i);
    T->aidx = kh_init(s2i);
    int indel_code = snptab_allele(T, "-");
    recvec_t *vec = NULL;
    int nvec = 0;
    int32_t *gt = NULL;
    int ngt_arr = 0;
    int64_t order = 0;
    const int64_t INDEL_ORDER = (int64_t)1 << 40;
    bcf1_t *rec = bcf_init();
    int ret;
    while ((ret = bcf_read(fp, h, rec)) == 0) {
        const char *chrom = bcf_seqname(h, rec);
        if (chrom_only && strcmp(chrom, chrom_only) != 0) continue;
        int ngt = bcf_get_genotypes(h, rec, &gt, &ngt_arr);
        if (ngt != 2) continue;
        if (bcf_gt_is_missing(gt[0]) || gt[1] == bcf_int32_vector_end || bcf_gt_is_missing(gt[1]))
            continue;
        int a = bcf_gt_allele(gt[0]), b = bcf_gt_allele(gt[1]);
        if (a == b) continue;
        if (a > b) { int t = a; a = b; b = t; }
        bcf_unpack(rec, BCF_UN_STR);
        if (b >= rec->n_allele) continue;
        int ci = snptab_chrom(T, chrom, 1);
        if (ci >= nvec) {
            vec = xrealloc(vec, sizeof(recvec_t) * (ci + 1));
            memset(vec + nvec, 0, sizeof(recvec_t) * (ci + 1 - nvec));
            nvec = ci + 1;
        }
        recvec_t *V = &vec[ci];
        const char *A = rec->d.allele[a], *B = rec->d.allele[b];
        if (is_acgt(A) && is_acgt(B)) {
            char al[3] = { A[0], B[0], 0 };
            GROW(V->r, V->n, V->m);
            V->r[V->n].pos = rec->pos;
            V->r[V->n].order = order++;
            V->r[V->n].code = (uint32_t)snptab_allele(T, al);
            V->n++;
        } else {
            int64_t p, end = rec->pos + (rec->rlen > 0 ? rec->rlen : 1);
            for (p = rec->pos; p < end; p++) {
                GROW(V->r, V->n, V->m);
                V->r[V->n].pos = p;
                V->r[V->n].order = INDEL_ORDER + order++;
                V->r[V->n].code = (uint32_t)indel_code;
                V->n++;
            }
        }
    }
    if (ret < -1) die("error reading %s", vcf);
    int c;
    for (c = 0; c < nvec; c++) {
        if (vec[c].n) {
            qsort(vec[c].r, vec[c].n, sizeof(sortrec_t), cmp_sortrec);
            finalize_chrom(&T->snps[c], vec[c].r, vec[c].n);
        }
        free(vec[c].r);
    }
    free(vec);
    free(gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(h);
    bcf_close(fp);
    if (T->n == 0) {
        printf("ERROR: no heterozygous sites found for sample %s\n", smp);
        fflush(stdout);
        exit(1);
    }
    free(smp);
    return T;
}

/* ------------------------------------------------------------------------ */
/* find_intersecting_snps                                                   */
/* ------------------------------------------------------------------------ */

typedef struct {
    int64_t total, skipped, no_snps, has_snps, total_snps, ref_match, no_match,
            toss_indel, toss_manysnps, toss_anomalous_phase, not_proper_pair,
            unmapped, non_refalt_base, flagged;
} counts_t;

typedef struct {
    int64_t aidx;     /* index among the aligned (M/=/X) bases */
    int64_t qpos;     /* position in the read */
    int64_t rpos;     /* reference position */
    uint32_t code;
} site_t;

typedef struct { site_t *s; int n, m; int64_t n_aligned; } sites_t;

static inline int64_t lower_bound(const int64_t *a, int64_t n, int64_t lo, int64_t x)
{
    int64_t hi = n;
    while (lo < hi) {
        int64_t mid = lo + ((hi - lo) >> 1);
        if (a[mid] < x) lo = mid + 1; else hi = mid;
    }
    return lo;
}

/* Equivalent to scanning read.get_aligned_pairs(matches_only=True) (pysam
 * semantics: M/=/X are aligned, I/S/P advance the read, D/N the reference). */
static void find_sites(const bam1_t *b, const snplist_t *L, sites_t *S)
{
    S->n = 0;
    S->n_aligned = 0;
    int64_t start = b->core.pos;
    if (start < 0 || b->core.n_cigar == 0 || L == NULL || L->n == 0) return;
    const int64_t *P = L->pos;
    int64_t n = L->n;
    int64_t i = lower_bound(P, n, 0, start);
    if (i >= n) return;
    const uint32_t *cig = bam_get_cigar(b);
    int64_t pos = start, qpos = 0, aidx = 0;
    uint32_t k;
    for (k = 0; k < b->core.n_cigar; k++) {
        int op = bam_cigar_op(cig[k]);
        int64_t l = bam_cigar_oplen(cig[k]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            int64_t e = pos + l;
            if (i < n && P[i] < pos) i = lower_bound(P, n, i, pos);
            while (i < n && P[i] < e) {
                int64_t d = P[i] - pos;
                GROW(S->s, S->n, S->m);
                S->s[S->n].aidx = aidx + d;
                S->s[S->n].qpos = qpos + d;
                S->s[S->n].rpos = P[i];
                S->s[S->n].code = L->code[i];
                S->n++;
                i++;
            }
            pos = e;
            qpos += l;
            aidx += l;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CPAD) {
            qpos += l;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            pos += l;
        }
    }
    S->n_aligned = aidx;
}

/* A growable set of equal-length sequences */
typedef struct { char *buf; int len; int n; size_t cap; } seqset_t;  /* cap in bytes */

static inline char *ss_get(seqset_t *S, int i) { return S->buf + (size_t)i * (S->len + 1); }

static void ss_reset(seqset_t *S, int len)
{
    S->len = len;
    S->n = 0;
}

static char *ss_push(seqset_t *S)
{
    size_t need = (size_t)(S->n + 1) * (S->len + 1);
    if (need > S->cap) {
        S->cap = S->cap ? S->cap * 2 : 1024;
        while (S->cap < need) S->cap *= 2;
        S->buf = xrealloc(S->buf, S->cap);
    }
    return ss_get(S, S->n++);
}

static int ss_push_copy(seqset_t *S, int from)
{
    char *dst = ss_push(S);
    memcpy(dst, ss_get(S, from), S->len + 1);
    return S->n - 1;
}

/* Fill S with the read's query_sequence */
static void ss_init_read(seqset_t *S, const bam1_t *b)
{
    int L = b->core.l_qseq, i;
    ss_reset(S, L);
    char *s = ss_push(S);
    const uint8_t *q = bam_get_seq(b);
    for (i = 0; i < L; i++) s[i] = seq_nt16_str[bam_seqi(q, i)];
    s[L] = 0;
}

static inline char seq_at(seqset_t *S, int which, int64_t pos, const bam1_t *b)
{
    if (b->core.l_qseq == 0)
        die("read %s has no sequence (TypeError in the Python version)", bam_get_qname(b));
    if (pos < 0 || pos >= S->len)
        die("read %s: string index out of range (IndexError in the Python version)",
            bam_get_qname(b));
    return ss_get(S, which)[pos];
}

/* Single-end allelic versions (default mode): drop reads overlapping an
 * indel, then reads with more than max_seqs combinations. Unphased: every
 * combination of the alleles; phased: every SNP switched to the other
 * allele. Returns the number of sequences (0: drop, 1: keep). */
static int single_read_seqs_checked(const bam1_t *b, sites_t *ST, const snptab_t *T,
                                    counts_t *D, int phased, seqset_t *S)
{
    int k, a, j;
    int64_t prod = 1;
    for (k = 0; k < ST->n; k++)
        if (T->is_indel[ST->s[k].code]) { D->toss_indel++; return 0; }
    for (k = 0; k < ST->n; k++) {
        prod *= T->alen[ST->s[k].code];
        if (prod > max_seqs) { D->toss_manysnps++; return 0; }
    }
    ss_init_read(S, b);
    if (phased) {
        int n1 = ss_push_copy(S, 0);
        for (k = 0; k < ST->n; k++) {
            site_t *s = &ST->s[k];
            const char *al = T->alleles[s->code];
            int alen = T->alen[s->code];
            D->total_snps++;
            const char *f = memchr(al, seq_at(S, 0, s->qpos, b), alen);
            if (!f) { D->no_match++; D->non_refalt_base++; return 0; }
            D->ref_match++;
            ss_get(S, n1)[s->qpos] = al[1 - (f - al)];
        }
        D->has_snps++;
        return 2;
    }
    for (k = 0; k < ST->n; k++) {
        site_t *s = &ST->s[k];
        const char *al = T->alleles[s->code];
        int alen = T->alen[s->code];
        char rb = seq_at(S, 0, s->qpos, b);
        D->total_snps++;
        if (memchr(al, rb, alen)) {
            D->ref_match++;
            for (a = 0; a < alen; a++) {
                if (al[a] == rb) continue;
                int cur = S->n;
                for (j = 0; j < cur; j++) {
                    int t = ss_push_copy(S, j);
                    ss_get(S, t)[s->qpos] = al[a];
                }
            }
        } else {
            D->no_match++;
        }
    }
    if (S->n == 1) D->no_snps++; else D->has_snps++;
    return S->n;
}

/* get_read_seqs(): returns number of sequences (0: drop, 1: keep) */
static int single_read_seqs(const bam1_t *b, sites_t *ST, const snptab_t *T,
                            counts_t *D, seqset_t *S)
{
    int k;
    int64_t prev = -2;
    ss_init_read(S, b);
    for (k = 0; k < ST->n; k++) {
        site_t *s = &ST->s[k];
        if (S->n > MAX_SEQS_PER_READ) {
            /* The original tests the count at every aligned base, after the
             * indel test. */
            if (s->aidx != prev + 1) { D->toss_manysnps++; return 0; }
            if (T->is_indel[s->code]) { D->toss_indel++; return 0; }
            D->toss_manysnps++;
            return 0;
        }
        if (T->is_indel[s->code]) { D->toss_indel++; return 0; }
        prev = s->aidx;
        D->total_snps++;
        const char *al = T->alleles[s->code];
        int alen = T->alen[s->code], a, j;
        char rb = seq_at(S, 0, s->qpos, b);
        if (memchr(al, rb, alen)) {
            D->ref_match++;
            for (a = 0; a < alen; a++) {
                if (al[a] == rb) continue;
                int cur = S->n;
                for (j = 0; j < cur; j++) {
                    int t = ss_push_copy(S, j);
                    ss_get(S, t)[s->qpos] = al[a];
                }
            }
        } else {
            D->no_match++;
        }
    }
    if (S->n > MAX_SEQS_PER_READ && prev + 1 < ST->n_aligned) {
        D->toss_manysnps++;
        return 0;
    }
    if (S->n == 1) D->no_snps++; else D->has_snps++;
    return S->n;
}

typedef struct { int64_t rpos; uint32_t code; int64_t p1, p2; } dsnp_t;

static int cmp_dsnp(const void *a, const void *b)
{
    const dsnp_t *x = a, *y = b;
    return x->rpos < y->rpos ? -1 : (x->rpos > y->rpos);
}

/* get_dual_read_seqs(): returns number of sequences per mate (0: drop). */
static int dual_read_seqs(const bam1_t *b1, const bam1_t *b2, sites_t *S1,
                          sites_t *S2, const snptab_t *T, counts_t *D,
                          int phased, seqset_t *Q1, seqset_t *Q2,
                          dsnp_t **dbuf, int *dm)
{
    int k, i, n = 0;
    for (k = 0; k < S1->n; k++)
        if (T->is_indel[S1->s[k].code]) { D->toss_indel++; return 0; }
    for (k = 0; k < S2->n; k++)
        if (T->is_indel[S2->s[k].code]) { D->toss_indel++; return 0; }
    /* snps dict in insertion order */
    if (*dm < S1->n + S2->n) {
        *dm = S1->n + S2->n + 16;
        *dbuf = xrealloc(*dbuf, sizeof(dsnp_t) * *dm);
    }
    dsnp_t *E = *dbuf;
    for (k = 0; k < S1->n; k++) {
        E[n].rpos = S1->s[k].rpos; E[n].code = S1->s[k].code;
        E[n].p1 = S1->s[k].qpos; E[n].p2 = -1; n++;
    }
    for (k = 0; k < S2->n; k++) {
        for (i = 0; i < n; i++) if (E[i].rpos == S2->s[k].rpos) break;
        if (i == n) { E[n].rpos = S2->s[k].rpos; E[n].p1 = -1; n++; }
        E[i].code = S2->s[k].code;
        E[i].p2 = S2->s[k].qpos;
    }
    {
        int64_t prod = 1;
        for (i = 0; i < n; i++) {
            prod *= T->alen[E[i].code];
            if (prod > MAX_SEQS_PER_READ) { D->toss_manysnps++; return 0; }
        }
    }
    ss_init_read(Q1, b1);
    ss_init_read(Q2, b2);
    if (n == 0) { D->no_snps++; return 1; }
    if (T->alen[E[0].code] > 2) { D->toss_manysnps++; return 0; }

    if (phased) {
        dsnp_t *F = xmalloc(sizeof(dsnp_t) * n);
        memcpy(F, E, sizeof(dsnp_t) * n);
        qsort(F, n, sizeof(dsnp_t), cmp_dsnp);
        int n1 = ss_push_copy(Q1, 0), n2 = ss_push_copy(Q2, 0);
        for (i = 0; i < n; i++) {
            const char *al = T->alleles[F[i].code];
            int alen = T->alen[F[i].code];
            if (!legacy) {
                /* The original never counted SNPs in phased mode. */
                D->total_snps++;
                if ((F[i].p1 >= 0 && !memchr(al, seq_at(Q1, 0, F[i].p1, b1), alen)) ||
                    (F[i].p2 >= 0 && !memchr(al, seq_at(Q2, 0, F[i].p2, b2), alen))) {
                    D->no_match++;
                    D->non_refalt_base++;
                    free(F);
                    return 0;
                }
                D->ref_match++;
            }
            if (F[i].p1 >= 0) {
                char c = seq_at(Q1, 0, F[i].p1, b1);
                const char *f = memchr(al, c, alen);
                if (!f) { free(F); return 0; }   /* 'non_refalt_base' (not counted) */
                int idx = 1 - (int)(f - al);
                if (idx < 0) idx += alen;
                ss_get(Q1, n1)[F[i].p1] = al[idx];
            }
            if (F[i].p2 >= 0) {
                char c = seq_at(Q2, 0, F[i].p2, b2);
                const char *f = memchr(al, c, alen);
                if (!f) { free(F); return 0; }
                int idx = 1 - (int)(f - al);
                if (idx < 0) idx += alen;
                ss_get(Q2, n2)[F[i].p2] = al[idx];
            }
            if (F[i].p1 >= 0 && F[i].p2 >= 0 &&
                seq_at(Q1, 0, F[i].p1, b1) != seq_at(Q2, 0, F[i].p2, b2)) {
                D->toss_anomalous_phase++;
                free(F);
                return 0;
            }
        }
        free(F);
        D->has_snps++;
        return 2;
    }

    /* Unphased.  'cur' mirrors the Python variables seq1/seq2, which are
     * rebound by the inner zip() loops to the last pre-existing sequences. */
    int cur = 0;
    for (i = 0; i < n; i++) {
        const char *al = T->alleles[E[i].code];
        int alen = T->alen[E[i].code], a, j;
        int64_t p1 = E[i].p1, p2 = E[i].p2;
        int L = Q1->n;
        D->total_snps++;
        if (p1 < 0) {
            if (memchr(al, seq_at(Q2, cur, p2, b2), alen)) D->ref_match++; else D->no_match++;
            for (a = 0; a < alen; a++) {
                if (al[a] == seq_at(Q2, cur, p2, b2)) continue;
                for (j = 0; j < L; j++) {
                    ss_push_copy(Q1, j);
                    int t = ss_push_copy(Q2, j);
                    ss_get(Q2, t)[p2] = al[a];
                }
                cur = L - 1;
            }
        } else if (p2 < 0) {
            if (memchr(al, seq_at(Q1, cur, p1, b1), alen)) D->ref_match++; else D->no_match++;
            for (a = 0; a < alen; a++) {
                if (al[a] == seq_at(Q1, cur, p1, b1)) continue;
                for (j = 0; j < L; j++) {
                    int t = ss_push_copy(Q1, j);
                    ss_get(Q1, t)[p1] = al[a];
                    ss_push_copy(Q2, j);
                }
                cur = L - 1;
            }
        } else {
            if (seq_at(Q1, cur, p1, b1) != seq_at(Q2, cur, p2, b2)) {
                D->toss_anomalous_phase++;
                return 0;
            }
            if (memchr(al, seq_at(Q2, cur, p2, b2), alen)) D->ref_match++; else D->no_match++;
            for (a = 0; a < alen; a++) {
                if (al[a] == seq_at(Q2, cur, p2, b2)) continue;
                for (j = 0; j < L; j++) {
                    int t = ss_push_copy(Q1, j);
                    ss_get(Q1, t)[p1] = al[a];
                    t = ss_push_copy(Q2, j);
                    ss_get(Q2, t)[p2] = al[a];
                }
                cur = L - 1;
            }
        }
    }
    if (Q1->n == 1) D->no_snps++; else D->has_snps++;
    return Q1->n;
}

/* Buffered text output to a (BGZF-compressed, gzip-compatible) file */
typedef struct { BGZF *fp; char *buf; size_t n, m; } fqout_t;

static void fq_open(fqout_t *f, const char *fname, hts_tpool *pool)
{
    f->fp = bgzf_open(fname, "w6");
    if (!f->fp) die("cannot write %s", fname);
    if (pool) bgzf_thread_pool(f->fp, pool, 0);
    f->m = 1 << 20;
    f->n = 0;
    f->buf = xmalloc(f->m);
}

static void fq_flush(fqout_t *f)
{
    if (f->n && bgzf_write(f->fp, f->buf, f->n) < 0) die("write error");
    f->n = 0;
}

static inline void fq_put(fqout_t *f, const char *s, size_t l)
{
    if (f->n + l > f->m) {
        fq_flush(f);
        if (l > f->m) { f->m = l * 2; f->buf = xrealloc(f->buf, f->m); }
    }
    memcpy(f->buf + f->n, s, l);
    f->n += l;
}

static void fq_close(fqout_t *f)
{
    fq_flush(f);
    if (bgzf_close(f->fp) < 0) die("error closing fastq");
    free(f->buf);
}

static unsigned char RC[256];

static void init_rc(void)
{
    int i;
    for (i = 0; i < 256; i++) RC[i] = (unsigned char)i;
    RC['A'] = 'T'; RC['T'] = 'A'; RC['C'] = 'G'; RC['G'] = 'C';
}

/* read.qual: qualities + 33, or "None" when absent */
static void qual_string(const bam1_t *b, char **buf, size_t *m, size_t *len)
{
    int L = b->core.l_qseq, i;
    const uint8_t *q = bam_get_qual(b);
    if (*m < (size_t)L + 8) { *m = L + 8; *buf = xrealloc(*buf, *m); }
    if (L == 0 || q[0] == 0xff) { strcpy(*buf, "None"); *len = 4; return; }
    if ((b->core.flag & BAM_FREVERSE) && !legacy) {
        /* the sequence is reverse complemented, so reverse the qualities */
        for (i = 0; i < L; i++) (*buf)[i] = (char)(q[L - 1 - i] + 33);
    } else {
        for (i = 0; i < L; i++) (*buf)[i] = (char)(q[i] + 33);
    }
    *len = L;
}

static void write_fq_record(fqout_t *f, const char *loc, size_t ll, const char *seq,
                            int L, int rev, const char *qual, size_t ql, char **tmp, size_t *tm)
{
    fq_put(f, "@", 1);
    fq_put(f, loc, ll);
    fq_put(f, "\n", 1);
    if (rev) {
        int i;
        if (*tm < (size_t)L + 1) { *tm = L + 1; *tmp = xrealloc(*tmp, *tm); }
        for (i = 0; i < L; i++) (*tmp)[i] = (char)RC[(unsigned char)seq[L - 1 - i]];
        fq_put(f, *tmp, L);
    } else {
        fq_put(f, seq, L);
    }
    fq_put(f, "\n+", 2);
    fq_put(f, loc, ll);
    fq_put(f, "\n", 1);
    fq_put(f, qual, ql);
    fq_put(f, "\n", 1);
}

KHASH_MAP_INIT_STR(rd, bam1_t *)

typedef struct {
    samFile *keep, *remap, *dropped;
    sam_hdr_t *hdr;
    fqout_t fq[2];
    int nfq;
    char *qbuf[2], *tmp;
    size_t qm[2], tm;
    int64_t remap_num;
} outputs_t;

static void bam_out(samFile *fp, sam_hdr_t *h, const bam1_t *b)
{
    if (sam_write1(fp, h, b) < 0) die("error writing BAM");
}

/* write_read_seqs() */
static void write_seqs(outputs_t *O, const bam1_t **reads, seqset_t **seqs, int nreads, int nseq)
{
    int r, i;
    if (nseq == 0) {
        for (r = 0; r < nreads; r++) bam_out(O->dropped, O->hdr, reads[r]);
        return;
    }
    if (nseq == 1) {
        for (r = 0; r < nreads; r++) bam_out(O->keep, O->hdr, reads[r]);
        return;
    }
    for (r = 0; r < nreads; r++) bam_out(O->remap, O->hdr, reads[r]);
    int64_t left = reads[0]->core.pos, right = reads[0]->core.pos;
    for (r = 1; r < nreads; r++) {
        if (reads[r]->core.pos < left) left = reads[r]->core.pos;
        if (reads[r]->core.pos > right) right = reads[r]->core.pos;
    }
    const bam1_t *last = reads[nreads - 1];
    const char *rname = last->core.tid >= 0 ? sam_hdr_tid2name(O->hdr, last->core.tid) : "None";
    char loc[4096];
    int ll = snprintf(loc, sizeof loc, "%lld:%s:%lld:%lld:%d", (long long)O->remap_num,
                      rname, (long long)left, (long long)right, nseq - 1);
    if (ll >= (int)sizeof loc) die("reference name too long");
    size_t ql[2];
    for (r = 0; r < nreads; r++) qual_string(reads[r], &O->qbuf[r], &O->qm[r], &ql[r]);
    for (i = 1; i < nseq; i++)
        for (r = 0; r < nreads; r++)
            write_fq_record(&O->fq[r], loc, ll, ss_get(seqs[r], i), seqs[r]->len,
                            (reads[r]->core.flag & BAM_FREVERSE) != 0, O->qbuf[r], ql[r],
                            &O->tmp, &O->tm);
    O->remap_num++;
}

static samFile *open_bam_out(const char *fname, sam_hdr_t *h, hts_tpool *pool)
{
    samFile *fp = sam_open(fname, "wb");
    if (!fp) die("cannot write %s", fname);
    if (pool) {
        htsThreadPool tp = { pool, 0 };
        hts_set_opt(fp, HTS_OPT_THREAD_POOL, &tp);
    }
    if (sam_hdr_write(fp, h) < 0) die("cannot write header to %s", fname);
    return fp;
}

static char *strcat3(const char *a, const char *b)
{
    char *p = xmalloc(strlen(a) + strlen(b) + 2);
    sprintf(p, "%s.%s", a, b);
    return p;
}

static void print_find_stats(counts_t *R, int64_t n_unpaired, int is_paired, const char *keep_only,
                             int phased)
{
    char p[64];
    int64_t total_pairs;
    printf("\n%s .... Finished!\n\nRUN STATISTICS:\n", now_str());
    if (is_paired) {
        total_pairs = R->total / 2;
        printf("  Total input reads: %lld pairs.\n", (long long)total_pairs);
        fmt_pct(p, n_unpaired, total_pairs);
        printf("  Unpaired reads: %lld (%s%%)\n", (long long)n_unpaired, p);
    } else {
        total_pairs = R->total;
        printf("  Total input reads: %lld\n", (long long)total_pairs);
    }
    if (keep_only) printf("  Reads on skipped chromosomes %lld\n", (long long)R->skipped);
    if (!legacy) {
        fmt_pct(p, R->flagged, total_pairs);
        printf("  Reads dropped [unmapped/secondary/supplementary/QC-fail%s]: %lld (%s%%)\n",
               is_paired ? "/mate unmapped" : "", (long long)R->flagged, p);
    }
    fmt_pct(p, R->no_snps, total_pairs);
    printf("  Reads with no SNPs: %lld (%s%%)\n", (long long)R->no_snps, p);
    fmt_pct(p, R->has_snps, total_pairs);
    printf("  Reads overlapping SNPs: %lld (%s%%)\n", (long long)R->has_snps, p);
    printf("  Total SNPs covered: %lld\n", (long long)R->total_snps);
    fmt_pct(p, R->ref_match, R->total_snps);
    printf("\tReference SNP matches: %lld (%s%%)\n", (long long)R->ref_match, p);
    fmt_pct(p, R->no_match, R->total_snps);
    printf("\tNon-reference SNP matches: %lld (%s%%)\n", (long long)R->no_match, p);
    fmt_pct(p, R->toss_indel, total_pairs);
    printf("  Reads dropped [INDEL]: %lld (%s%%)\n", (long long)R->toss_indel, p);
    fmt_pct(p, R->toss_manysnps, total_pairs);
    printf("  Reads dropped [too many SNPs]: %lld (%s%%)\n", (long long)R->toss_manysnps, p);
    if (is_paired) {
        fmt_pct(p, R->toss_anomalous_phase, total_pairs);
        printf("  Reads dropped [anomalous pair]: %lld (%s%%)\n",
               (long long)R->toss_anomalous_phase, p);
    }
    if (phased && !legacy) {
        fmt_pct(p, R->non_refalt_base, total_pairs);
        printf("  Reads dropped [non ref/alt base]: %lld (%s%%)\n",
               (long long)R->non_refalt_base, p);
    }
    if (is_paired) {
        fmt_pct(p, R->not_proper_pair, total_pairs);
        printf("  Reads dropped [not proper pair]: %lld (%s%%)\n",
               (long long)R->not_proper_pair, p);
    }
    printf("\n");
    fflush(stdout);
}

/* simple pool of bam1_t records */
typedef struct { bam1_t **b; int n, m; } bpool_t;
static bam1_t *bp_get(bpool_t *p) { return p->n ? p->b[--p->n] : bam_init1(); }
static void bp_put(bpool_t *p, bam1_t *b) { GROW(p->b, p->n, p->m); p->b[p->n++] = b; }

static void usage_find(void)
{
    fprintf(stderr,
"usage: hornet find_intersecting_snps [options] infile snp_dir\n"
"       hornet find_intersecting_snps [options] --vcf FILE [--sample NAME] infile\n"
"  -p, --paired_end        reads are paired-end (default is single)\n"
"  -C, --limit-to-chrom    limit loading of SNPs to the specified chromosome\n"
"  -s, --sorted            input is coordinate sorted (accepted, unused)\n"
"  -P, --phased-snp-data   SNPs are phased (e.g. F1 hybrid): one swapped read per read\n"
"  -t, --threads           threads for BAM/fastq (de)compression and SNP loading (default 4)\n"
"      --progressbar       accepted for compatibility (ignored)\n"
"      --vcf FILE          use the heterozygous sites of a sample in a VCF/BCF\n"
"      --sample NAME       sample in --vcf (default: the only sample)\n"
"      --max-seqs N        drop reads (pairs) with more than N allelic combinations\n"
"                          (default 64, i.e. at most 6 SNPs)\n"
"      --legacy            reproduce the original output exactly, bugs included\n"
"                          (no read flag filters, --max-seqs 32)\n");
    exit(2);
}

static int cmd_find(int argc, char **argv)
{
    static struct option lo[] = {
        {"paired_end", no_argument, 0, 'p'},
        {"limit-to-chrom", required_argument, 0, 'C'},
        {"sorted", no_argument, 0, 's'},
        {"phased-snp-data", no_argument, 0, 'P'},
        {"progressbar", no_argument, 0, 1},
        {"threads", required_argument, 0, 't'},
        {"processes", required_argument, 0, 'j'},
        {"legacy", no_argument, 0, 2},
        {"vcf", required_argument, 0, 3},
        {"sample", required_argument, 0, 4},
        {"max-seqs", required_argument, 0, 5},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}};
    int c, is_paired = 0, phased = 0, threads = 4, max_seqs_set = 0;
    const char *keep_only = NULL, *vcf = NULL, *sample = NULL;
    while ((c = getopt_long(argc, argv, "pC:sPt:j:h", lo, NULL)) >= 0) {
        switch (c) {
        case 'p': is_paired = 1; break;
        case 'C': keep_only = optarg; break;
        case 's': break;
        case 'P': phased = 1; break;
        case 't': threads = atoi(optarg); break;
        case 'j': if (atoi(optarg) > threads) threads = atoi(optarg); break;
        case 1: break;
        case 2: legacy = 1; break;
        case 3: vcf = optarg; break;
        case 4: sample = optarg; break;
        case 5:
            max_seqs = atoi(optarg);
            max_seqs_set = 1;
            if (max_seqs < 1) die("--max-seqs must be >= 1");
            break;
        default: usage_find();
        }
    }
    if (argc - optind != (vcf ? 1 : 2)) usage_find();
    const char *infile = argv[optind], *snpdir = vcf ? NULL : argv[optind + 1];
    if (threads < 1) threads = 1;
    if (legacy) {
        if (vcf || max_seqs_set) die("--legacy cannot be combined with --vcf or --max-seqs");
        max_seqs = 32;
    }

    samFile *in = sam_open(infile, "r");
    if (!in) die("cannot open %s", infile);
    sam_hdr_t *hdr = sam_hdr_read(in);
    if (!hdr) die("cannot read header of %s", infile);

    snptab_t *T = vcf ? get_snps_vcf(vcf, sample, keep_only) : get_snps(snpdir, keep_only, threads);
    printf("%s .... Done with SNPs.\n", now_str());
    fflush(stdout);

    hts_tpool *pool = threads > 1 ? hts_tpool_init(threads) : NULL;
    if (pool) {
        htsThreadPool tp = { pool, 0 };
        hts_set_opt(in, HTS_OPT_THREAD_POOL, &tp);
    }

    /* fname.rsplit('.', 1)[0] */
    char *base = xstrdup(infile);
    char *dot = strrchr(base, '.');
    if (dot) *dot = 0;

    outputs_t O;
    memset(&O, 0, sizeof O);
    O.hdr = hdr;
    O.remap_num = 1;
    char *fn;
    O.keep = open_bam_out(fn = strcat3(base, "keep.bam"), hdr, pool); free(fn);
    O.remap = open_bam_out(fn = strcat3(base, "to.remap.bam"), hdr, pool); free(fn);
    O.dropped = open_bam_out(fn = strcat3(base, "dropped.bam"), hdr, pool); free(fn);
    if (is_paired) {
        fq_open(&O.fq[0], fn = strcat3(base, "remap.fq1.gz"), pool); free(fn);
        fq_open(&O.fq[1], fn = strcat3(base, "remap.fq2.gz"), pool); free(fn);
        O.nfq = 2;
    } else {
        fq_open(&O.fq[0], fn = strcat3(base, "remap.fq.gz"), pool); free(fn);
        O.nfq = 1;
    }

    /* SNP list for each tid of the BAM */
    int ntid = sam_hdr_nref(hdr), t;
    snplist_t **tabs = xmalloc(sizeof(snplist_t *) * (ntid + 1));
    for (t = 0; t < ntid; t++) {
        int ci = snptab_chrom(T, sam_hdr_tid2name(hdr, t), 0);
        tabs[t] = (ci >= 0 && T->snps[ci].n > 0) ? &T->snps[ci] : NULL;
    }
    const char *last_name = ntid > 0 ? sam_hdr_tid2name(hdr, ntid - 1) : NULL;

    counts_t D;
    memset(&D, 0, sizeof D);
    khash_t(rd) *unpaired[2] = { kh_init(rd), kh_init(rd) };
    bpool_t bp = { 0, 0, 0 };
    sites_t S1 = { 0 }, S2 = { 0 };
    seqset_t Q1 = { 0 }, Q2 = { 0 };
    dsnp_t *dbuf = NULL;
    int dm = 0, ret;
    bam1_t *b = bp_get(&bp);

    while ((ret = sam_read1(in, hdr, b)) >= 0) {
        int tid = b->core.tid;
        if (keep_only) {
            const char *cn;
            if (tid >= 0) cn = sam_hdr_tid2name(hdr, tid);
            else if (tid == -1 && last_name) cn = last_name;   /* chroms[-1] */
            else die("IndexError: tuple index out of range");
            if (strcmp(cn, keep_only) != 0) { D.skipped++; continue; }
        }
        D.total++;
        if (!is_paired) {
            if (legacy) {
                if (phased) die("NotImplementedError: Not yet implemented. It should work with paired end reads, though");
            } else if (b->core.flag & FLAG_DROP_SINGLE) {
                D.flagged++;
                continue;
            }
            find_sites(b, tid >= 0 ? tabs[tid] : NULL, &S1);
            if (S1.n == 0) {
                D.no_snps++;
                bam_out(O.keep, hdr, b);
                continue;
            }
            int ns = legacy ? single_read_seqs(b, &S1, T, &D, &Q1)
                            : single_read_seqs_checked(b, &S1, T, &D, phased, &Q1);
            const bam1_t *rr[1] = { b };
            seqset_t *qq[1] = { &Q1 };
            write_seqs(&O, rr, qq, 1, ns);
            continue;
        }
        uint16_t flag = b->core.flag;
        if (!legacy) {
            if (flag & FLAG_DROP_PAIRED) { D.flagged++; continue; }
            if (b->core.mtid != tid) { D.not_proper_pair++; continue; }
        }
        if (!(flag & BAM_FPROPER_PAIR)) { D.not_proper_pair++; continue; }
        int slot_self = (flag & BAM_FREAD2) ? 1 : 0;
        int slot_other = (flag & BAM_FREAD1) ? 1 : 0;
        const char *qname = bam_get_qname(b);
        khint_t k = kh_get(rd, unpaired[slot_other], qname);
        if (k == kh_end(unpaired[slot_other])) {
            khash_t(rd) *h = unpaired[slot_self];
            k = kh_get(rd, h, qname);
            if (k != kh_end(h)) {
                bp_put(&bp, kh_val(h, k));
            } else {
                k = kh_put(rd, h, xstrdup(qname), &ret);
            }
            kh_val(h, k) = b;
            b = bp_get(&bp);
            continue;
        }
        bam1_t *other = kh_val(unpaired[slot_other], k);
        free((char *)kh_key(unpaired[slot_other], k));
        kh_del(rd, unpaired[slot_other], k);
        bam1_t *both[2] = { NULL, NULL };
        both[slot_self] = b;
        both[slot_other] = other;
        if (!both[0] || !both[1])
            die("read %s is flagged as both or neither of read1/read2 "
                "(AttributeError in the Python version)", qname);
        bam1_t *r1 = both[0], *r2 = both[1];
        const bam1_t *rr[2] = { r1, r2 };
        if ((r1->core.flag & BAM_FUNMAP) || (r2->core.flag & BAM_FUNMAP)) {
            D.unmapped++;
            bam_out(O.dropped, hdr, r1);
            bam_out(O.dropped, hdr, r2);
        } else {
            int tid1 = r1->core.tid;
            snplist_t *L = tid1 >= 0 ? tabs[tid1] : NULL;
            find_sites(r1, L, &S1);
            find_sites(r2, L, &S2);
            if (S1.n == 0 && S2.n == 0) {
                D.no_snps++;
                bam_out(O.keep, hdr, r1);
                bam_out(O.keep, hdr, r2);
            } else {
                int ns = dual_read_seqs(r1, r2, &S1, &S2, T, &D, phased, &Q1, &Q2, &dbuf, &dm);
                seqset_t *qq[2] = { &Q1, &Q2 };
                write_seqs(&O, rr, qq, 2, ns);
            }
        }
        bp_put(&bp, other);
    }
    if (ret < -1) die("error reading %s", infile);

    int64_t n_unpaired = kh_size(unpaired[0]) + kh_size(unpaired[1]);
    if (sam_close(O.keep) < 0 || sam_close(O.remap) < 0 || sam_close(O.dropped) < 0)
        die("error closing output BAM");
    for (c = 0; c < O.nfq; c++) fq_close(&O.fq[c]);
    sam_close(in);
    if (pool) hts_tpool_destroy(pool);
    print_find_stats(&D, n_unpaired, is_paired, keep_only, phased);
    return 0;
}

/* ------------------------------------------------------------------------ */
/* filter_remapped_reads                                                    */
/* ------------------------------------------------------------------------ */

static int64_t field_int(const char *s, const char *qname)
{
    int64_t v;
    if (parse_pyint(s, &v) < 0)
        die("invalid read name %s (ValueError in the Python version)", qname);
    return v;
}

static int cmp_i64(const void *a, const void *b)
{
    int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
    return x < y ? -1 : (x > y);
}

static void usage_filter(void)
{
    fprintf(stderr,
"usage: hornet filter_remapped_reads [-p] [-t THREADS] [--legacy] to_remap_bam remap_bam keep_bam\n"
"  remap_bam can be in any order (--legacy: original algorithm, needs a name-sorted BAM)\n");
    exit(2);
}


/* WASP-like filter (default). For each remapped record of a read (pair):
 * no secondary/supplementary alignment (any discards the read), same
 * chromosome, start, strand and CIGAR (hard clips ignored) as the original,
 * for pairs a proper pair with the mate at the original mate position; the
 * number of such records must equal the number of versions written (times
 * two for pairs). The remapped BAM can be in any order. */
static uint64_t cigar_key(const bam1_t *b)
{
    uint64_t h = 1469598103934665603ULL;
    const uint32_t *c = bam_get_cigar(b);
    uint32_t k;
    for (k = 0; k < b->core.n_cigar; k++) {
        if (bam_cigar_op(c[k]) == BAM_CHARD_CLIP) continue;
        h = (h ^ c[k]) * 1099511628211ULL;
    }
    return h;
}

static int filter_strict(const char *to_remap_path, const char *remap_path,
                         const char *keep_path, int is_paired, hts_tpool *pool)
{
    htsThreadPool tp = { pool, 0 };
    int n = is_paired ? 2 : 1, ret;
    int64_t k = 0, nid, i;
    bam1_t *b = bam_init1();

    /* pass 1: original reads, n per remap ID */
    samFile *fp = sam_open(to_remap_path, "r");
    if (!fp) die("cannot open %s", to_remap_path);
    if (pool) hts_set_opt(fp, HTS_OPT_THREAD_POOL, &tp);
    sam_hdr_t *h1 = sam_hdr_read(fp);
    if (!h1) die("cannot read header of %s", to_remap_path);
    int32_t *otid = NULL;
    int64_t *opos = NULL;
    uint8_t *orev = NULL;
    uint64_t *ocig = NULL;
    int64_t m = 0;
    while ((ret = sam_read1(fp, h1, b)) >= 0) {
        int64_t j = k + n;   /* index (id * n + mate), ids start at 1 */
        if (j + n > m) {
            m = m ? m * 2 : 1 << 16;
            while (m < j + n) m *= 2;
            otid = xrealloc(otid, sizeof(int32_t) * (m / n + 1));
            opos = xrealloc(opos, sizeof(int64_t) * m);
            orev = xrealloc(orev, m);
            ocig = xrealloc(ocig, sizeof(uint64_t) * m);
        }
        if (k % n == 0) otid[k / n + 1] = b->core.tid;
        opos[j] = b->core.pos;
        orev[j] = (b->core.flag & BAM_FREVERSE) ? 1 : 0;
        ocig[j] = cigar_key(b);
        k++;
    }
    if (ret < -1) die("error reading %s", to_remap_path);
    if (k % n) die("%s has an odd number of reads", to_remap_path);
    nid = k / n;
    sam_close(fp);

    /* pass 2: remapped reads */
    samFile *rp = sam_open(remap_path, "r");
    if (!rp) die("cannot open %s", remap_path);
    if (pool) hts_set_opt(rp, HTS_OPT_THREAD_POOL, &tp);
    sam_hdr_t *h2 = sam_hdr_read(rp);
    if (!h2) die("cannot read header of %s", remap_path);
    int nref2 = sam_hdr_nref(h2), t;
    int *tidmap = xmalloc(sizeof(int) * (nref2 + 1));
    for (t = 0; t < nref2; t++) {
        tidmap[t] = sam_hdr_name2tid(h1, sam_hdr_tid2name(h2, t));
        if (tidmap[t] < 0) tidmap[t] = -2;
    }
    int32_t *expected = calloc(nid + 1, sizeof(int32_t));
    int32_t *correct = calloc(nid + 1, sizeof(int32_t));
    uint8_t *moved = calloc(nid + 1, 1), *multi = calloc(nid + 1, 1);
    char tok[64];
    while ((ret = sam_read1(rp, h2, b)) >= 0) {
        const char *q = bam_get_qname(b);
        const char *c1 = strchr(q, ':'), *cl = strrchr(q, ':');
        int64_t rid, num;
        size_t l;
        if (!c1) die("unexpected read name %s in %s (not from find_intersecting_snps?)", q, remap_path);
        l = c1 - q;
        if (l >= sizeof tok) die("unexpected read name %s", q);
        memcpy(tok, q, l); tok[l] = 0;
        if (parse_pyint(tok, &rid) < 0) die("unexpected read name %s in %s (not from find_intersecting_snps?)", q, remap_path);
        l = strcspn(cl + 1, "/");
        if (l >= sizeof tok) die("unexpected read name %s", q);
        memcpy(tok, cl + 1, l); tok[l] = 0;
        if (parse_pyint(tok, &num) < 0) die("unexpected read name %s in %s (not from find_intersecting_snps?)", q, remap_path);
        if (rid < 1 || rid > nid) die("read %s has no counterpart in %s", q, to_remap_path);
        expected[rid] = (int32_t)(num * n);
        uint16_t flag = b->core.flag;
        if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) { multi[rid] = 1; continue; }
        if (flag & BAM_FUNMAP) continue;
        int mate = (n == 2 && (flag & BAM_FREAD2)) ? 1 : 0;
        int64_t j = rid * n + mate;
        int tid = b->core.tid >= 0 ? tidmap[b->core.tid] : -3;
        int ok = tid == otid[rid] && b->core.pos == opos[j] &&
                 ((flag & BAM_FREVERSE) ? 1 : 0) == orev[j] && cigar_key(b) == ocig[j];
        if (ok && n == 2)
            ok = (flag & BAM_FPROPER_PAIR) && !(flag & BAM_FMUNMAP) && b->core.mtid >= 0 &&
                 tidmap[b->core.mtid] == tid && b->core.mpos == opos[rid * n + 1 - mate];
        if (ok) correct[rid]++; else moved[rid] = 1;
    }
    if (ret < -1) die("error reading %s", remap_path);
    sam_close(rp);

    /* pass 3: write kept reads */
    fp = sam_open(to_remap_path, "r");
    if (!fp) die("cannot open %s", to_remap_path);
    if (pool) hts_set_opt(fp, HTS_OPT_THREAD_POOL, &tp);
    sam_hdr_destroy(h1);
    h1 = sam_hdr_read(fp);
    samFile *keep = open_bam_out(keep_path, h1, pool);
    int64_t n_keep = 0, n_moved = 0, n_multi = 0, n_missing = 0;
    for (i = 0; sam_read1(fp, h1, b) >= 0; i++) {
        int64_t rid = i / n + 1;
        int last = (i + 1) % n == 0;
        if (multi[rid]) n_multi += last;
        else if (moved[rid]) n_moved += last;
        else if (expected[rid] == 0 || correct[rid] != expected[rid]) n_missing += last;
        else { bam_out(keep, h1, b); n_keep += last; }
    }
    if (sam_close(keep) < 0) die("error closing %s", keep_path);
    sam_close(fp);
    fprintf(stderr, "%s Finished\n\n", now_str());
    fprintf(stderr, "RUN STATISTICS:\n\tTotal remapped read (pair)s: %lld\n\t"
            "Read (pair)s remapped to the correct position: %lld %.2f%%\n"
            "\tDiscarded, remapped to another position, strand or CIGAR: %lld\n"
            "\tDiscarded, secondary/supplementary alignments: %lld\n"
            "\tDiscarded, unmapped or missing remapped reads: %lld\n",
            (long long)nid, (long long)n_keep, nid ? (double)n_keep / (double)nid * 100 : 0.0,
            (long long)n_moved, (long long)n_multi, (long long)n_missing);
    return 0;
}

static int cmd_filter(int argc, char **argv)
{
    static struct option lo[] = {
        {"threads", required_argument, 0, 't'}, {"legacy", no_argument, 0, 2},
        {"help", no_argument, 0, 'h'}, {0, 0, 0, 0}};
    int c, is_paired = 0, threads = 4;
    while ((c = getopt_long(argc, argv, "pt:h", lo, NULL)) >= 0) {
        switch (c) {
        case 'p': is_paired = 1; break;
        case 't': threads = atoi(optarg); break;
        case 2: legacy = 1; break;
        default: usage_filter();
        }
    }
    if (argc - optind != 3) usage_filter();
    if (threads < 1) threads = 1;
    fprintf(stderr, "%s Starting run\n", now_str());

    hts_tpool *pool = threads > 1 ? hts_tpool_init(threads) : NULL;
    if (!legacy) {
        filter_strict(argv[optind], argv[optind + 1], argv[optind + 2], is_paired, pool);
        if (pool) hts_tpool_destroy(pool);
        return 0;
    }
    htsThreadPool tp = { pool, 0 };
    samFile *to_remap = sam_open(argv[optind], "rb");
    if (!to_remap) die("cannot open %s", argv[optind]);
    samFile *remap = sam_open(argv[optind + 1], "rb");
    if (!remap) die("cannot open %s", argv[optind + 1]);
    if (pool) { hts_set_opt(to_remap, HTS_OPT_THREAD_POOL, &tp); hts_set_opt(remap, HTS_OPT_THREAD_POOL, &tp); }
    sam_hdr_t *h1 = sam_hdr_read(to_remap), *h2 = sam_hdr_read(remap);
    if (!h1 || !h2) die("cannot read BAM header");
    samFile *keep = open_bam_out(argv[optind + 2], h1, pool);

    int64_t *correct = NULL, *nums = NULL;
    int64_t nc = 0, mc = 0, nn = 0, mn = 0;
    int64_t prev_read_num = 0, counter = 1;
    bam1_t *b = bam_init1();
    char *qn = NULL;
    size_t qm = 0;
    int ret = sam_read1(remap, h2, b);
    if (ret < 0) die("StopIteration: %s has no reads", argv[optind + 1]);
    while (ret >= 0) {
        /* qname.strip().split(":") */
        const char *q = bam_get_qname(b);
        size_t ql = strlen(q);
        if (qm < ql + 1) { qm = ql + 1; qn = xrealloc(qn, qm); }
        memcpy(qn, q, ql + 1);
        char *f[64];
        int nf = 0;
        char *s = qn, *e;
        while (*s && py_isspace((unsigned char)*s)) s++;
        e = s + strlen(s);
        while (e > s && py_isspace((unsigned char)e[-1])) *--e = 0;
        f[nf++] = s;
        for (; *s; s++)
            if (*s == ':') { *s = 0; if (nf < 64) f[nf++] = s + 1; else die("too many fields in %s", q); }
        if (nf < 2) die("IndexError: read name %s", q);
        const char *chrm = f[1];
        int pf = ((b->core.flag & BAM_FREVERSE) && is_paired) ? 3 : 2;
        if (nf <= pf) die("IndexError: read name %s", q);
        int64_t pos = field_int(f[pf], q);
        int64_t read_num = field_int(f[0], q);
        int64_t num = field_int(f[nf - 1], q);

        if (read_num != prev_read_num) {
            if (read_num != counter) {
                int64_t skipped = read_num - counter, i;
                for (i = 0; i < skipped; i++) { GROW(nums, nn, mn); nums[nn++] = 1; }
                counter = read_num;
            }
            GROW(nums, nn, mn);
            nums[nn++] = is_paired ? num * 2 : num;
            prev_read_num = read_num;
            counter++;
        }
        if (b->core.tid != -1 && b->core.pos == pos &&
            strcmp(sam_hdr_tid2name(h2, b->core.tid), chrm) == 0) {
            const uint32_t *cig = bam_get_cigar(b);
            uint32_t k;
            int dels = 0;
            for (k = 0; k < b->core.n_cigar; k++) {
                int op = bam_cigar_op(cig[k]);
                if (op != BAM_CMATCH && op != BAM_CREF_SKIP && op != BAM_CSOFT_CLIP) dels++;
            }
            if (dels == 0) { GROW(correct, nc, mc); correct[nc++] = read_num; }
        }
        ret = sam_read1(remap, h2, b);
    }
    /* correct_maps.sort() -- the list is usually already sorted */
    {
        int64_t i;
        for (i = 1; i < nc; i++) if (correct[i] < correct[i-1]) break;
        if (i < nc) qsort(correct, nc, sizeof(int64_t), cmp_i64);
    }

    bam1_t *orig = bam_init1(), *second = bam_init1();
    if (sam_read1(to_remap, h1, orig) < 0) die("StopIteration: %s has no reads", argv[optind]);
    if (nn == 0) die("IndexError: list index out of range");
    int64_t orig_num = nums[0];
    int64_t line_num = 1, map_indx = 0, correct_n = 0, total_correct = 0;
    int end_of_file = 0;
    while (!end_of_file && map_indx < nc && line_num <= correct[nc - 1]) {
        if (line_num != correct[map_indx]) {
            if (correct_n == orig_num) {
                total_correct++;
                bam_out(keep, h1, orig);
                if (is_paired) {
                    if (sam_read1(to_remap, h1, orig) < 0)
                        die("ValueError: File ended unexpectedly (no pair found).");
                    bam_out(keep, h1, orig);
                }
            } else if (is_paired || legacy) {
                /* skip the mate of the discarded read (the original also did
                 * this for single-end data, dropping the next read) */
                if (sam_read1(to_remap, h1, second) < 0) { end_of_file = 1; break; }
            }
            if (sam_read1(to_remap, h1, orig) < 0) {
                end_of_file = 1;
            } else {
                if (line_num >= nn) die("IndexError: list index out of range");
                orig_num = nums[line_num];
            }
            line_num++;
            correct_n = 0;
        } else {
            correct_n++;
            map_indx++;
        }
    }
    if (correct_n == orig_num) {
        total_correct++;
        bam_out(keep, h1, orig);
        if (is_paired) {
            if (sam_read1(to_remap, h1, orig) < 0) {
                fprintf(stderr, "File ended unexpectedly (no pair found).");
                sam_close(keep);
                exit(0);
            }
            bam_out(keep, h1, orig);
        }
    }
    if (sam_close(keep) < 0) die("error closing %s", argv[optind + 2]);
    sam_close(to_remap);
    sam_close(remap);
    if (pool) hts_tpool_destroy(pool);
    fprintf(stderr, "%s Finished\n\n", now_str());
    fprintf(stderr, "RUN STATISTICS:\n\tTotal remapped read (pair)s: %lld\n\t"
            "Read (pair)s remapped to the correct position: %lld %.2f%%\n",
            (long long)line_num, (long long)total_correct,
            (double)total_correct / (double)line_num * 100);
    return 0;
}

/* ------------------------------------------------------------------------ */

static void usage(void)
{
    fprintf(stderr,
"hornet: unbiased allele-specific read mapping (C implementation)\n\n"
"usage: hornet <command> [options]\n\n"
"commands:\n"
"  find_intersecting_snps   find reads overlapping SNPs, write reads to remap\n"
"  filter_remapped_reads    keep reads that remapped exactly like the original\n\n"
"Remove duplicates with Picard MarkDuplicates --DUPLICATE_SCORING_STRATEGY RANDOM\n"
"(see README).\n");
    exit(2);
}

int main(int argc, char **argv)
{
    if (argc < 2) usage();
    const char *cmd = argv[1];
    init_rc();
    argv++;
    argc--;
    if (!strcmp(cmd, "find_intersecting_snps") || !strcmp(cmd, "find_intersecting_snps.py") ||
        !strcmp(cmd, "find"))
        return cmd_find(argc, argv);
    if (!strcmp(cmd, "filter_remapped_reads") || !strcmp(cmd, "filter_remapped_reads.py") ||
        !strcmp(cmd, "filter"))
        return cmd_filter(argc, argv);
    usage();
    return 2;
}

#!/bin/bash
# 1. With --legacy, the optimised Python scripts and the C binary must produce
#    exactly the same output as the original Hornet code (upstream e642ac0).
# 2. In the default mode (bugs fixed, read flag filters, --max-seqs, VCF
#    input, WASP-like filter) both must produce exactly the same output as
#    the simple reference implementations ref_find.py / ref_filter_strict.py.
#
#   tests/run_equivalence.sh [workdir]
#
# Needs: python3 with pysam (+ numpy, pandas optional), samtools, git.
# The C binary is tested if c/hornet has been built.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(dirname "$HERE")
WORK=${1:-$(mktemp -d)}
ORIG_COMMIT=e642ac0
mkdir -p "$WORK/orig"
for f in find_intersecting_snps.py filter_remapped_reads.py; do
    git -C "$REPO" show $ORIG_COMMIT:mapping/$f > "$WORK/orig/$f"
done
HORNET_C="$REPO/c/hornet"
[ -x "$HORNET_C" ] || HORNET_C=""
fail=0
pass=0

check() { # description, then a command that returns 0 on success
    local d=$1; shift
    if "$@"; then pass=$((pass+1)); echo "ok    $d"; else fail=$((fail+1)); echo "FAIL  $d"; fi
}
same_bam() { [ -s "$1" ] && [ -s "$2" ] && cmp -s <(samtools view --no-PG -h "$1") <(samtools view --no-PG -h "$2"); }
same_gz() { [ -s "$1" ] && [ -s "$2" ] && cmp -s <(zcat "$1") <(zcat "$2"); }
stats() { grep -v -E 'Loading|Done with|Finished|Starting' "$1"; }

# find_intersecting_snps: compare one implementation with the original
fis() { # data bam "flags" impl-name impl-cmd...
    local D=$1 B=$2 F=$3 name=$4; shift 4
    local b=${B%.*} W=$D/run_$(echo "$B $F $name" | tr -c 'a-zA-Z0-9\n' '_')
    rm -rf "$W"; mkdir -p "$W/orig" "$W/new"
    for x in orig new; do
        ln -s "$D/$B" "$W/$x/$B"
        [ -e "$D/$B.bai" ] && ln -s "$D/$B.bai" "$W/$x/$B.bai"
    done
    (cd "$W/orig" && python3 "$WORK/orig/find_intersecting_snps.py" $F $B "$D/snps" > stdout 2> stderr)
    (cd "$W/new" && "$@" $F $B "$D/snps" > stdout 2> stderr)
    local ok=0
    for f in keep.bam to.remap.bam dropped.bam; do same_bam "$W/orig/$b.$f" "$W/new/$b.$f" || ok=1; done
    for f in remap.fq.gz remap.fq1.gz remap.fq2.gz; do
        [ -e "$W/orig/$b.$f" ] || continue
        same_gz "$W/orig/$b.$f" "$W/new/$b.$f" || ok=1
    done
    stats "$W/orig/stdout" > "$W/o.txt"; stats "$W/new/stdout" > "$W/n.txt"
    if grep -q ZeroDivisionError "$W/orig/stderr"; then
        # The original crashes while printing the statistics in phased mode
        # (no SNPs are counted); compare what it managed to print.
        head -n "$(wc -l < "$W/o.txt")" "$W/n.txt" > "$W/n2.txt"; mv "$W/n2.txt" "$W/n.txt"
    fi
    cmp -s "$W/o.txt" "$W/n.txt" || ok=1
    # keep the remap outputs of the original for the filter test
    echo "$W/orig" > "$D/last_orig_$b$( [ -n "$F" ] && echo "$F" | tr -c 'a-zA-Z0-9\n' '_')"
    return $ok
}

echo "work dir: $WORK"
cd "$WORK"
python3 "$HERE/gen_test.py" data_a 1 20000 > /dev/null
python3 "$HERE/gen_test.py" data_b 7 20000 0 > /dev/null   # no read-name collisions
# consolidated SNP file with CRLF line endings
mkdir -p data_c/snps
for f in data_a/snps/*.snps.txt.gz; do
    c=$(basename "$f" .snps.txt.gz); zcat "$f" | awk -v c="$c" '{print c"\t"$0}'
done | sed 's/$/\r/' | gzip > data_c/snps/all.txt.gz
ln -sf "$WORK/data_a/pe.bam" data_c/pe.bam; ln -sf "$WORK/data_a/se.bam" data_c/se.bam

PYFIS="python3 $REPO/mapping/find_intersecting_snps.py"
impls=("python:$PYFIS --legacy")
[ -n "$HORNET_C" ] && impls+=("C:$HORNET_C find_intersecting_snps --legacy")

for impl in "${impls[@]}"; do
    name=${impl%%:*}; cmd=${impl#*:}
    for D in data_a data_b; do
        for B in pe.bam pe.sorted.bam; do
            check "$name find_intersecting_snps $D/$B -p"    fis "$WORK/$D" $B "-p" $name $cmd
            check "$name find_intersecting_snps $D/$B -p -P" fis "$WORK/$D" $B "-p -P" $name $cmd
        done
        for B in se.bam se.sorted.bam; do
            check "$name find_intersecting_snps $D/$B"       fis "$WORK/$D" $B "" $name $cmd
        done
    done
    check "$name find_intersecting_snps all.txt.gz -p"       fis "$WORK/data_c" pe.bam "-p" $name $cmd
    check "$name find_intersecting_snps all.txt.gz -p -C"    fis "$WORK/data_c" pe.bam "-p -C chr2" $name $cmd
    check "$name find_intersecting_snps -C chr1 (se)"        fis "$WORK/data_a" se.bam "-C chr1" $name $cmd
done
for B in pe.sorted.bam se.sorted.bam; do
    F="-p"; [ $B = se.sorted.bam ] && F=""
    check "python find_intersecting_snps -j 8 data_b/$B" fis "$WORK/data_b" $B "$F" pyj $PYFIS --legacy -j 8
done

# filter_remapped_reads
for mode in pe se; do
    D=$WORK/data_b; P=-p; [ $mode = se ] && P=""
    R=$(cat "$D/last_orig_$mode.sorted$( [ -n "$P" ] && echo "$P" | tr -c 'a-zA-Z0-9\n' '_')")
    if [ $mode = pe ]; then fq="$R/pe.sorted.remap.fq1.gz $R/pe.sorted.remap.fq2.gz"; else fq="$R/se.sorted.remap.fq.gz"; fi
    python3 "$HERE/sim_remap.py" "$D/$mode.sorted.bam" "$D/$mode.remap.bam" 3 $fq
    python3 "$WORK/orig/filter_remapped_reads.py" $P "$R/$mode.sorted.to.remap.bam" "$D/$mode.remap.bam" "$D/f.orig.bam" 2> "$D/f.orig.err"
    python3 "$REPO/mapping/filter_remapped_reads.py" --legacy $P "$R/$mode.sorted.to.remap.bam" "$D/$mode.remap.bam" "$D/f.py.bam" 2> "$D/f.py.err"
    check "python filter_remapped_reads $mode" eval 'same_bam "$D/f.orig.bam" "$D/f.py.bam" && cmp -s <(stats "$D/f.orig.err") <(stats "$D/f.py.err")'
    if [ -n "$HORNET_C" ]; then
        "$HORNET_C" filter_remapped_reads --legacy $P "$R/$mode.sorted.to.remap.bam" "$D/$mode.remap.bam" "$D/f.c.bam" 2> "$D/f.c.err"
        check "C filter_remapped_reads $mode" eval 'same_bam "$D/f.orig.bam" "$D/f.c.bam" && cmp -s <(stats "$D/f.orig.err") <(stats "$D/f.c.err")'
    fi
done

# ===========================================================================
# 2. Default mode vs the reference implementations
# ===========================================================================
LAST_W=
ref_fis() { # data bam "flags" text|vcf snp_dir|vcf_file [sample]
    local D=$1 B=$2 F=$3 kind=$4 src=$5 smp=${6:-}
    local b=${B%.*} x ok=0
    local W=$D/def_$(echo "$B $F $kind $(basename $src) $smp" | tr -c 'a-zA-Z0-9\n' '_')
    rm -rf "$W"; mkdir -p "$W/ref" "$W/py" "$W/c"
    for x in py c; do
        ln -s "$D/$B" "$W/$x/$B"
        [ -e "$D/$B.bai" ] && ln -s "$D/$B.bai" "$W/$x/$B.bai"
    done
    local RA HA
    if [ $kind = text ]; then RA="--snp-dir $src"; HA=""; else
        RA="--vcf $src"; HA="--vcf $src"
        [ -n "$smp" ] && { RA="$RA --sample $smp"; HA="$HA --sample $smp"; }
    fi
    python3 "$HERE/ref_find.py" $F $RA "$D/$B" "$W/ref/$b" || return 1
    if [ $kind = text ]; then
        (cd "$W/py" && $PYFIS $F $B $src > stdout 2> stderr)
        [ -n "$HORNET_C" ] && (cd "$W/c" && "$HORNET_C" find_intersecting_snps $F $B $src > stdout 2> stderr)
    else
        (cd "$W/py" && $PYFIS $F $HA $B > stdout 2> stderr)
        [ -n "$HORNET_C" ] && (cd "$W/c" && "$HORNET_C" find_intersecting_snps $F $HA $B > stdout 2> stderr)
    fi
    local impls=py; [ -n "$HORNET_C" ] && impls="py c"
    for x in $impls; do
        for f in keep to.remap dropped; do
            same_bam "$W/ref/$b.$f.sam" "$W/$x/$b.$f.bam" || { ok=1; echo "   $x $f.bam differs"; }
        done
        for f in remap.fq remap.fq1 remap.fq2; do
            [ -e "$W/ref/$b.$f" ] || continue
            cmp -s "$W/ref/$b.$f" <(zcat "$W/$x/$b.$f.gz") || { ok=1; echo "   $x $f differs"; }
        done
        stats "$W/$x/stdout" > "$W/$x.txt"
    done
    [ -n "$HORNET_C" ] && { cmp -s "$W/py.txt" "$W/c.txt" || { ok=1; echo "   statistics differ"; }; }
    LAST_W=$W
    return $ok
}

for D in data_a data_b; do
    S=$WORK/$D/snps
    check "default find $D/pe.bam -p"                ref_fis "$WORK/$D" pe.bam "-p" text $S
    check "default find $D/pe.bam -p -P"             ref_fis "$WORK/$D" pe.bam "-p -P" text $S
    check "default find $D/pe.sorted.bam -p"         ref_fis "$WORK/$D" pe.sorted.bam "-p" text $S
    check "default find $D/pe.bam -p --max-seqs 8"   ref_fis "$WORK/$D" pe.bam "-p --max-seqs 8" text $S
    check "default find $D/pe.bam -p -C chr3"        ref_fis "$WORK/$D" pe.bam "-p -C chr3" text $S
    check "default find $D/se.bam"                   ref_fis "$WORK/$D" se.bam "" text $S
    check "default find $D/se.bam -P"                ref_fis "$WORK/$D" se.bam "-P" text $S
    check "default find $D/se.bam --max-seqs 4"      ref_fis "$WORK/$D" se.bam "--max-seqs 4" text $S
done
check "default find all.txt.gz -p"                   ref_fis "$WORK/data_c" pe.bam "-p" text $WORK/data_c/snps

# VCF input: same as the equivalent SNP files
python3 "$HERE/make_vcf.py" data_a/snps data_a/pe.bam data_a/test.vcf.gz data_a/snps_from_vcf 5
for F in "-p" "-p -P" ""; do
    B=pe.bam; [ -z "$F" ] && B=se.bam
    check "default find --vcf $B $F" ref_fis "$WORK/data_a" $B "$F" vcf $WORK/data_a/test.vcf.gz sample1
    VW=$LAST_W
    check "default find snp files equivalent to --vcf, $B $F" ref_fis "$WORK/data_a" $B "$F" text $WORK/data_a/snps_from_vcf
    TW=$LAST_W
    check "  --vcf output = SNP-file output ($B $F)" eval '(b=${B%.*}; for f in keep to.remap dropped; do same_bam $VW/py/$b.$f.bam $TW/py/$b.$f.bam || exit 1; done; for f in $VW/py/*.fq*.gz; do same_gz $f $TW/py/$(basename $f) || exit 1; done)'
done

# Python per-chromosome parallel mode = serial
for B in pe.sorted.bam se.sorted.bam; do
    F="-p"; [ $B = se.sorted.bam ] && F=""
    check "default find -j 8 = serial, data_b/$B" eval '(
        W=$WORK/data_b/def_j_$B; rm -rf $W; mkdir -p $W; ln -s $WORK/data_b/$B $WORK/data_b/$B.bai $W/;
        (cd $W && $PYFIS $F -j 8 $B $WORK/data_b/snps > stdout 2> stderr);
        ref_fis "$WORK/data_b" $B "$F" text $WORK/data_b/snps > /dev/null; S=$LAST_W/py; b=${B%.*};
        for f in keep to.remap dropped; do same_bam $S/$b.$f.bam $W/$b.$f.bam || exit 1; done;
        for f in $S/*.fq*.gz; do same_gz $f $W/$(basename $f) || exit 1; done)'
done

# filter_remapped_reads (WASP-like) vs reference; remapped BAM in random order
for mode in pe se; do
    D=$WORK/data_b; P=-p; B=pe.bam; [ $mode = se ] && { P=""; B=se.bam; }
    ref_fis "$D" $B "$P" text $D/snps > /dev/null; R=$LAST_W/py; b=${B%.*}
    if [ $mode = pe ]; then fq="$R/$b.remap.fq1.gz $R/$b.remap.fq2.gz"; else fq="$R/$b.remap.fq.gz"; fi
    python3 "$HERE/sim_remap_strict.py" "$R/$b.to.remap.bam" "$D/strict.$mode.remap.bam" 11 $fq
    python3 "$HERE/ref_filter_strict.py" $P "$R/$b.to.remap.bam" "$D/strict.$mode.remap.bam" "$D/strict.$mode.ref.sam"
    python3 "$REPO/mapping/filter_remapped_reads.py" $P "$R/$b.to.remap.bam" "$D/strict.$mode.remap.bam" "$D/strict.$mode.py.bam" 2> "$D/strict.$mode.py.err"
    check "default filter_remapped_reads $mode (python) = reference" same_bam "$D/strict.$mode.ref.sam" "$D/strict.$mode.py.bam"
    if [ -n "$HORNET_C" ]; then
        "$HORNET_C" filter_remapped_reads -t 4 $P "$R/$b.to.remap.bam" "$D/strict.$mode.remap.bam" "$D/strict.$mode.c.bam" 2> "$D/strict.$mode.c.err"
        check "default filter_remapped_reads $mode (C) = reference" eval 'same_bam "$D/strict.$mode.ref.sam" "$D/strict.$mode.c.bam" && cmp -s <(stats "$D/strict.$mode.py.err") <(stats "$D/strict.$mode.c.err")'
    fi
    check "  ... and some reads were discarded for each reason ($mode)" eval '(
        e=$D/strict.$mode.py.err; for k in "another position" "secondary" "unmapped or missing"; do
        n=$(grep "$k" $e | grep -o "[0-9]*$"); [ "${n:-0}" -gt 0 ] || exit 1; done)'
done

echo
echo "$pass passed, $fail failed"
[ $fail = 0 ]

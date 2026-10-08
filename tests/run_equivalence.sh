#!/bin/bash
# 1. With --legacy, the optimised Python scripts and the C binary must produce
#    exactly the same output as the original Hornet code (upstream e642ac0).
# 2. By default (bugs fixed) Python and C must agree with each other, and the
#    differences from the original must be exactly the intended fixes.
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
for f in find_intersecting_snps.py filter_remapped_reads.py rmdup.py rmdup_pe.py; do
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
same_bam() { cmp -s <(samtools view --no-PG -h "$1") <(samtools view --no-PG -h "$2"); }
same_gz() { cmp -s <(zcat "$1") <(zcat "$2"); }
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

# ---------------------------------------------------------------------------
# Default (fixed) behaviour
# ---------------------------------------------------------------------------
fis_fixed() { # data bam "flags"
    local D=$1 B=$2 F=$3
    local b=${B%.*} W=$D/fixed_$(echo "$B $F" | tr -c 'a-zA-Z0-9\n' '_') x
    rm -rf "$W"; mkdir -p "$W/orig" "$W/py" "$W/c"
    for x in orig py c; do
        ln -s "$D/$B" "$W/$x/$B"
        [ -e "$D/$B.bai" ] && ln -s "$D/$B.bai" "$W/$x/$B.bai"
    done
    (cd "$W/orig" && python3 "$WORK/orig/find_intersecting_snps.py" $F $B "$D/snps" > stdout 2> stderr)
    (cd "$W/py" && $PYFIS $F $B "$D/snps" > stdout 2> stderr)
    local impl=py
    if [ -n "$HORNET_C" ]; then
        (cd "$W/c" && "$HORNET_C" find_intersecting_snps $F $B "$D/snps" > stdout 2> stderr)
        impl="py c"
    fi
    local ok=0 fq
    for x in $impl; do
        for f in keep.bam to.remap.bam dropped.bam; do same_bam "$W/orig/$b.$f" "$W/$x/$b.$f" || ok=1; done
        # remap fastqs: only the qualities of reverse-strand reads change
        if [ -e "$W/orig/$b.remap.fq1.gz" ]; then
            python3 "$HERE/check_fastq_quals.py" "$W/orig/$b.to.remap.bam" \
                "$W/orig/$b.remap.fq1.gz" "$W/$x/$b.remap.fq1.gz" \
                "$W/orig/$b.remap.fq2.gz" "$W/$x/$b.remap.fq2.gz" > /dev/null || ok=1
        else
            python3 "$HERE/check_fastq_quals.py" "$W/orig/$b.to.remap.bam" \
                "$W/orig/$b.remap.fq.gz" "$W/$x/$b.remap.fq.gz" > /dev/null || ok=1
        fi
        stats "$W/$x/stdout" > "$W/$x.txt"
    done
    stats "$W/orig/stdout" > "$W/orig.txt"
    if [ -n "$HORNET_C" ]; then
        for fq in "$W"/py/*.fq*.gz; do same_gz "$fq" "$W/c/$(basename "$fq")" || ok=1; done
        cmp -s "$W/py.txt" "$W/c.txt" || ok=1
    fi
    if echo "$F" | grep -q -- -P; then
        # phased statistics are now counted (the original divided by zero)
        grep -q "Total SNPs covered: 0$" "$W/py.txt" && ok=1
        grep -q "non ref/alt base" "$W/py.txt" || ok=1
    else
        cmp -s "$W/orig.txt" "$W/py.txt" || ok=1
    fi
    return $ok
}
for D in data_a data_b; do
    for B in pe.bam pe.sorted.bam; do
        check "fixed find_intersecting_snps $D/$B -p"    fis_fixed "$WORK/$D" $B "-p"
        check "fixed find_intersecting_snps $D/$B -p -P" fis_fixed "$WORK/$D" $B "-p -P"
    done
    check "fixed find_intersecting_snps $D/se.bam"       fis_fixed "$WORK/$D" se.bam ""
done
check "fixed find_intersecting_snps -j 8 (python) = serial" eval '
    W=$WORK/data_b/fixed_j; rm -rf $W; mkdir -p $W; ln -s $WORK/data_b/pe.sorted.bam* $W/;
    (cd $W && $PYFIS -p -j 8 pe.sorted.bam $WORK/data_b/snps > stdout 2>&1) &&
    S=$WORK/data_b/fixed_pe_sorted_bam__p/py &&
    same_bam $S/pe.sorted.keep.bam $W/pe.sorted.keep.bam &&
    same_bam $S/pe.sorted.to.remap.bam $W/pe.sorted.to.remap.bam &&
    same_gz $S/pe.sorted.remap.fq1.gz $W/pe.sorted.remap.fq1.gz &&
    same_gz $S/pe.sorted.remap.fq2.gz $W/pe.sorted.remap.fq2.gz'

# filter_remapped_reads: compare with a straightforward reference implementation
for mode in pe se; do
    D=$WORK/data_b; P=-p; [ $mode = se ] && P=""
    R=$(cat "$D/last_orig_$mode.sorted$( [ -n "$P" ] && echo "$P" | tr -c 'a-zA-Z0-9\n' '_')")
    python3 "$HERE/ref_filter.py" $P "$R/$mode.sorted.to.remap.bam" "$D/$mode.remap.bam" "$D/f.ref.sam"
    python3 "$REPO/mapping/filter_remapped_reads.py" $P "$R/$mode.sorted.to.remap.bam" "$D/$mode.remap.bam" "$D/f.pyfix.bam" 2> "$D/f.pyfix.err"
    check "fixed python filter_remapped_reads $mode = reference" same_bam "$D/f.ref.sam" "$D/f.pyfix.bam"
    if [ -n "$HORNET_C" ]; then
        "$HORNET_C" filter_remapped_reads $P "$R/$mode.sorted.to.remap.bam" "$D/$mode.remap.bam" "$D/f.cfix.bam" 2> "$D/f.cfix.err"
        check "fixed C filter_remapped_reads $mode = reference" eval 'same_bam "$D/f.ref.sam" "$D/f.cfix.bam" && cmp -s <(stats "$D/f.pyfix.err") <(stats "$D/f.cfix.err")'
    fi
    if [ $mode = pe ]; then
        python3 "$WORK/orig/filter_remapped_reads.py" -p "$R/$mode.sorted.to.remap.bam" "$D/$mode.remap.bam" "$D/f.orig.pe.bam" 2> /dev/null
        check "fixed filter_remapped_reads pe = original" same_bam "$D/f.orig.pe.bam" "$D/f.pyfix.bam"
    fi
done

# rmdup / rmdup_pe with the same random seed
python3 "$HERE/mkdup.py" data_b/pe.sorted.bam data_b/dup.pe.bam
python3 "$HERE/mkdup.py" data_b/se.sorted.bam data_b/dup.se.bam
for seed in 1 42 -3 123456789012; do
    for t in rmdup:se rmdup_pe:pe; do
        prog=${t%%:*}; m=${t#*:}
        python3 "$HERE/shim_run.py" $seed "$WORK/orig/$prog.py" data_b/dup.$m.bam data_b/$prog.orig.bam
        python3 "$REPO/mapping/$prog.py" --seed $seed data_b/dup.$m.bam data_b/$prog.py.bam
        check "python $prog --seed $seed" same_bam data_b/$prog.orig.bam data_b/$prog.py.bam
        if [ -n "$HORNET_C" ]; then
            "$HORNET_C" $prog --seed $seed data_b/dup.$m.bam data_b/$prog.c.bam
            check "C $prog --seed $seed" same_bam data_b/$prog.orig.bam data_b/$prog.c.bam
        fi
    done
done

echo
echo "$pass passed, $fail failed"
[ $fail = 0 ]

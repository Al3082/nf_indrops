/*
 * Pre-flight FASTQ sanity check.
 *
 * For each of R1, R2, R4:
 *   - total read count
 *   - integrity: total lines divisible by 4, seq and qual same length for every read
 *   - length distribution (first 100000 reads)
 *
 * Then barcode match rates:
 *   - CB1: R2 exact sequence vs cb_whitelist1
 *   - CB2: R4 first 8 bp vs cb_whitelist2
 *
 * WARNING lines are printed when integrity fails or match rate is below 5%.
 * Does not fail the pipeline — results are published for inspection.
 */
process check_fastq {
    tag "preflight ${label}"
    label "preflight"

    publishDir "${params.output_dir}/preflight", mode: 'copy'

    input:
        tuple val(label), path(r1), path(r2), path(r4)
        path cb_whitelist1
        path cb_whitelist2_sense
        path cb_whitelist2_rc

    output:
        path "${label}.preflight.txt"

    script:
    """
    check_file() {
        local tag=\$1 f=\$2 expected=\$3
        echo "[\$tag] \$(basename \$f)  (expected \$expected)"

        # Single pass: read count, line-count divisibility, seq/qual length parity
        zcat "\$f" | awk '
            NR%4==2 { seq_len = length(\$0); n++ }
            NR%4==0 { if (length(\$0) != seq_len) n_bad++ }
            END {
                printf "  Total reads : %d\\n", n
                if (NR%4 != 0)
                    printf "  CORRUPT: %d lines total (not divisible by 4 — truncated?)\\n", NR
                if (n_bad+0 > 0)
                    printf "  CORRUPT: %d reads have mismatched seq/qual lengths\\n", n_bad
                else if (NR%4 == 0)
                    printf "  Integrity   : OK\\n"
            }
        '

        # Length distribution (first 100000 reads)
        echo "  Length distribution (first 100000 reads):"
        zcat "\$f" | awk 'NR%4==2' | head -100000 \\
            | awk '{c[length(\$0)]++} END{for(l in c) print l, c[l]}' \\
            | sort -n | awk '{printf "    %d bp: %d reads\\n", \$1, \$2}'
        echo ""
    }

    {
    echo "=== Pre-flight check: ${label} ==="
    echo "Date: \$(date)"
    echo ""

    check_file R1 ${r1} "~95 bp, gene read"
    check_file R2 ${r2} "8 bp, CB1"
    check_file R4 ${r4} "14 bp, CB2+UMI"

    echo "--- CB1 match check (R2 exact sequence vs whitelist1) ---"
    zcat ${r2} | awk 'NR%4==2' | head -100000 > _r2_sample.txt
    n_r2=\$(wc -l < _r2_sample.txt)
    n_cb1=\$(grep -cFxf ${cb_whitelist1} _r2_sample.txt || true)
    pct_cb1=\$(awk -v nr="\$n_r2" -v nc="\$n_cb1" 'BEGIN{printf "%.1f", (nr>0) ? 100*nc/nr : 0}')
    echo "  Sampled : \$n_r2"
    echo "  Matches : \$n_cb1"
    echo "  Rate    : \${pct_cb1}%"
    awk -v p="\$pct_cb1" 'BEGIN{if (p+0 < 5) print "  WARNING: <5% match — R2 may not contain CB1 barcodes (wrong file or orientation?)"}'
    echo ""

    echo "--- CB1 match check (R2 exact sequence vs whitelist2) ---"
    zcat ${r2} | awk 'NR%4==2' | head -100000 > _r2_sample.txt
    n_r2=\$(wc -l < _r2_sample.txt)
    n_cb1=\$(grep -cFxf ${cb_whitelist2_sense} _r2_sample.txt || true) 
    pct_cb1=\$(awk -v nr="\$n_r2" -v nc="\$n_cb1" 'BEGIN{printf "%.1f", (nr>0) ? 100*nc/nr : 0}')
    echo "  Sampled : \$n_r2"
    echo "  Matches : \$n_cb1"
    echo "  Rate    : \${pct_cb1}%"
    awk -v p="\$pct_cb1" 'BEGIN{if (p+0 < 5) print "  WARNING: <5% match — R2 may not contain CB1 barcodes (wrong file or orientation?)"}'
    echo ""


    echo "--- CB2 match check (R4 first 8 bp vs whitelist2 — both orientations) ---"
    zcat ${r4} | awk 'NR%4==2{print substr(\$0,1,8)}' | head -100000 > _r4_cb2_sample.txt
    n_r4=\$(wc -l < _r4_cb2_sample.txt)

    n_cb2_sense=\$(grep -cFxf ${cb_whitelist2_sense} _r4_cb2_sample.txt || true)
    pct_cb2_sense=\$(awk -v nr="\$n_r4" -v nc="\$n_cb2_sense" 'BEGIN{printf "%.1f", (nr>0) ? 100*nc/nr : 0}')

    n_cb2_rc=\$(grep -cFxf ${cb_whitelist2_rc} _r4_cb2_sample.txt || true)
    pct_cb2_rc=\$(awk -v nr="\$n_r4" -v nc="\$n_cb2_rc" 'BEGIN{printf "%.1f", (nr>0) ? 100*nc/nr : 0}')

    echo "  Sampled          : \$n_r4"
    echo "  Sense matches    : \$n_cb2_sense  (\${pct_cb2_sense}%)"
    echo "  RC matches       : \$n_cb2_rc  (\${pct_cb2_rc}%)"
    awk -v s="\$pct_cb2_sense" -v r="\$pct_cb2_rc" 'BEGIN{if (s+0 >= r+0) print "  Best orientation : sense"; else print "  Best orientation : RC"}'
    awk -v s="\$pct_cb2_sense" -v r="\$pct_cb2_rc" 'BEGIN{if (s+0 < 5 && r+0 < 5) print "  WARNING: <5% match for both orientations — wrong file or incorrect whitelist?"}'

    } > ${label}.preflight.txt

    cat ${label}.preflight.txt >&2
    """
}

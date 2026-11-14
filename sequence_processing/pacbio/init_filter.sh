#!/bin/bash
# FASTQ Filter Script
# Filters sequences based on paired substring matching:
# - If sequence contains (substring_A AND substring_B) -> output as-is
# - If sequence contains (substring_C AND substring_D) -> output reverse complement
# Both go to a single output file
#
# Usage: ./filter_fastq_pairs.sh <input.fastq>
set -euo pipefail
INPUT_FASTQ="$1"
SUBSTRING_A="CGGGAAGGCAATGCAAG"
SUBSTRING_B="AATGGCGAATGGACGC"
SUBSTRING_C="CTTGCATTGCCTTCCCG"
SUBSTRING_D="GCGTCCATTCGCCATT"
OUTPUT_PREFIX="../pacbio-filtered"
OUTPUT="${OUTPUT_PREFIX}.fastq"
if [[ ! -f "$INPUT_FASTQ" ]]; then
    echo "Error: Input file '$INPUT_FASTQ' not found." >&2
    exit 1
fi
echo "Filtering: $INPUT_FASTQ"
echo "Forward: '$SUBSTRING_A' + '$SUBSTRING_B'"
echo "Reverse: '$SUBSTRING_C' + '$SUBSTRING_D' (will be reverse complemented)"
echo "Output: $OUTPUT"
# Use pigz if available for parallel decompression, otherwise zcat
if [[ "$INPUT_FASTQ" == *.gz ]]; then
    if command -v pigz &> /dev/null; then
        CAT_CMD="pigz -dc"
    else
        CAT_CMD="zcat"
    fi
else
    CAT_CMD="cat"
fi
$CAT_CMD "$INPUT_FASTQ" | awk -v subA="$SUBSTRING_A" -v subB="$SUBSTRING_B" \
                              -v subC="$SUBSTRING_C" -v subD="$SUBSTRING_D" \
                              -v outfile="$OUTPUT" '
BEGIN {
    total = 0; nAB = 0; nCD = 0; nBoth = 0
    # Build complement lookup table
    comp["A"] = "T"; comp["T"] = "A"
    comp["C"] = "G"; comp["G"] = "C"
    comp["a"] = "t"; comp["t"] = "a"
    comp["c"] = "g"; comp["g"] = "c"
    comp["N"] = "N"; comp["n"] = "n"
}

function revcomp(s,    i, c, result) {
    result = ""
    for (i = length(s); i >= 1; i--) {
        c = substr(s, i, 1)
        result = result (c in comp ? comp[c] : c)
    }
    return result
}

NR % 4 == 1 { header = $0 }
NR % 4 == 2 { seq = $0 }
NR % 4 == 3 { plus = $0 }
NR % 4 == 0 {
    total++
    qual = $0
    hasAB = (index(seq, subA) && index(seq, subB))
    hasCD = (index(seq, subC) && index(seq, subD))
    
    if (hasAB && hasCD) {
        nBoth++
    } else if (hasAB) {
        nAB++
        print header ORS seq ORS plus ORS qual > outfile
    } else if (hasCD) {
        nCD++
        # Reverse complement sequence and reverse quality string
        rc_seq = revcomp(seq)
        rev_qual = ""
        for (i = length(qual); i >= 1; i--) {
            rev_qual = rev_qual substr(qual, i, 1)
        }
        print header ORS rc_seq ORS plus ORS rev_qual > outfile
    }
}
END {
    printf "\nTotal: %d | Fwd: %d | Rev(rc): %d | Both: %d | Neither: %d\n", total, nAB, nCD, nBoth, total-nAB-nCD-nBoth
    if (nBoth > 0) printf "Note: %d reads matched both pairs (written to neither)\n", nBoth
}'

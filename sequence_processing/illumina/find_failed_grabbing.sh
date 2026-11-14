#!/bin/bash

# Check which jobs failed by looking for ERROR or OOM in err/out files
# Pass the original job array ID as argument: bash find_failed.sh 12345678

ORIGINAL_JOB_ID=$1
LOG_DIR="../logs/grab_barcodes"
FAILED_IDS=""

for i in $(seq 1 188); do
    errfile="${LOG_DIR}/grabbarcodes_${ORIGINAL_JOB_ID}_${i}.err"
    outfile="${LOG_DIR}/grabbarcodes_${ORIGINAL_JOB_ID}_${i}.out"

    failed=0
    # Also check err file for OOM signals
    if [ -f "$outfile" ]; then
        if grep -qi "ERROR|CANCELLED" "$errfile"; then
            failed=1
        fi
    fi

    # Also check err file for OOM signals
    if [ -f "$errfile" ]; then
        if grep -qi "oom\|out of memory\|exceeded memory\|Killed\|MemoryError\|CANCELLED" "$errfile"; then
            failed=1
        fi
    fi

    if [ $failed -eq 1 ]; then
        if [ -z "$FAILED_IDS" ]; then
            FAILED_IDS="$i"
        else
            FAILED_IDS="${FAILED_IDS},$i"
        fi
    fi
done

echo "Failed task IDs:"
echo "$FAILED_IDS"

# Write to file for use in resubmission
echo "$FAILED_IDS" > ../failed_task_ids.txt
echo ""
echo "Count: $(echo $FAILED_IDS | tr ',' '\n' | wc -l) failed jobs"
echo "Saved to failed_task_ids.txt"

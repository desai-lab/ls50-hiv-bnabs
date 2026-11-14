#!/usr/bin/env python3
"""
Copies all {experiment}_init_counts.csv files to an output directory as
{experiment}_init_counts_merged.csv, then streams through the Undetermined
files and appends matching rows (with the 'experiment' column stripped).

Designed for large files: never loads full files into memory.
"""

import csv
import io
import os
import shutil
import sys
from collections import defaultdict

# --- Configuration ---
INPUT_DIR = "."       # directory containing all the CSV files
OUTPUT_DIR = "merged"  # output directory for merged files

UNDETERMINED_FILES = [
    "Undetermined_S0_L003_init_counts.csv",
    "Undetermined_S0_L004_init_counts.csv",
]


def find_experiment_files(input_dir, undetermined_files):
    """Find all *_init_counts.csv files that are NOT Undetermined files."""
    undet_set = set(undetermined_files)
    experiment_files = []
    for fname in os.listdir(input_dir):
        if fname.endswith("_init_counts.csv") and fname not in undet_set:
            experiment_files.append(fname)
    return sorted(experiment_files)


def extract_experiment_name(filename):
    """Extract experiment name: 'LA_9_5_3_init_counts.csv' -> 'LA_9_5_3'."""
    suffix = "_init_counts.csv"
    if filename.endswith(suffix):
        return filename[: -len(suffix)]
    return None


def merged_filename(experiment_name):
    """Build merged output filename from experiment name."""
    return f"{experiment_name}_init_counts_merged.csv"


def file_ends_with_newline(filepath):
    """Check if a file ends with a newline character."""
    with open(filepath, "rb") as f:
        f.seek(0, 2)
        size = f.tell()
        if size == 0:
            return True
        f.seek(-1, 2)
        last_byte = f.read(1)
        return last_byte in (b"\n", b"\r")


def ensure_trailing_newline(filepath):
    """Append a newline to a file if it doesn't already end with one."""
    if not file_ends_with_newline(filepath):
        with open(filepath, "ab") as f:
            f.write(b"\n")
        return True
    return False


def stream_undetermined_file(undet_path, experiment_to_outpath, append_counts):
    """
    Stream through one Undetermined file, appending each row to the
    appropriate experiment output file.

    Uses a dict of open file handles (one per experiment encountered)
    to avoid repeatedly opening/closing files. The number of open handles
    equals the number of distinct experiments in the Undetermined file,
    which should be manageable.

    Args:
        undet_path: path to the Undetermined CSV
        experiment_to_outpath: dict mapping experiment name -> output file path
        append_counts: dict to accumulate per-experiment append counts
    """
    undet_name = os.path.basename(undet_path)
    open_handles = {}   # experiment_name -> (file_handle, csv_writer)
    unmatched = defaultdict(int)
    row_count = 0

    try:
        with open(undet_path, "r", newline="") as f:
            reader = csv.reader(f)
            header = next(reader)

            # Validate header
            expected = ["experiment", "umi1", "umi2", "bc", "counts"]
            if header != expected:
                print(f"WARNING: {undet_name} has unexpected header: {header}")
                print(f"  Expected: {expected}")

            exp_col = header.index("experiment")
            other_cols = [i for i in range(len(header)) if i != exp_col]

            for row in reader:
                if len(row) != len(header):
                    print(f"WARNING: Skipping malformed row in {undet_name}: {row}")
                    continue

                experiment = row[exp_col]
                row_count += 1

                if experiment not in experiment_to_outpath:
                    unmatched[experiment] += 1
                    continue

                # Open file handle on first encounter
                if experiment not in open_handles:
                    outpath = experiment_to_outpath[experiment]
                    fh = open(outpath, "a", newline="")
                    writer = csv.writer(fh)
                    open_handles[experiment] = (fh, writer)

                _, writer = open_handles[experiment]
                stripped_row = [row[i] for i in other_cols]
                writer.writerow(stripped_row)
                append_counts[experiment] += 1

    finally:
        # Close all open file handles
        for fh, _ in open_handles.values():
            fh.close()

    print(f"Streamed {row_count} rows from {undet_name}")

    if unmatched:
        total_unmatched = sum(unmatched.values())
        print(f"  WARNING: {total_unmatched} rows across {len(unmatched)} experiments "
              f"had no matching experiment file:")
        for exp in sorted(unmatched):
            print(f"    {exp} ({unmatched[exp]} rows)")

    return unmatched


def main():
    input_dir = INPUT_DIR
    output_dir = OUTPUT_DIR

    if len(sys.argv) >= 2:
        input_dir = sys.argv[1]
    if len(sys.argv) >= 3:
        output_dir = sys.argv[2]

    print(f"Input directory:  {input_dir}")
    print(f"Output directory: {output_dir}")

    # Find experiment files
    experiment_files = find_experiment_files(input_dir, UNDETERMINED_FILES)
    if not experiment_files:
        print("ERROR: No experiment files found!")
        sys.exit(1)
    print(f"Found {len(experiment_files)} experiment files.")

    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Copy experiment files to output dir with merged names,
    #         and ensure each ends with a newline.
    experiment_to_outpath = {}
    for fname in experiment_files:
        exp_name = extract_experiment_name(fname)
        src = os.path.join(input_dir, fname)
        dst = os.path.join(output_dir, merged_filename(exp_name))

        shutil.copy2(src, dst)
        if ensure_trailing_newline(dst):
            print(f"  NOTE: {fname} did not end with newline, added one.")

        experiment_to_outpath[exp_name] = dst

    print(f"Copied {len(experiment_files)} files to {output_dir}/")

    # Step 2: Stream through each Undetermined file and append rows.
    append_counts = defaultdict(int)
    all_unmatched = defaultdict(int)

    for undet_file in UNDETERMINED_FILES:
        undet_path = os.path.join(input_dir, undet_file)
        if not os.path.exists(undet_path):
            print(f"WARNING: {undet_file} not found, skipping.")
            continue
        unmatched = stream_undetermined_file(
            undet_path, experiment_to_outpath, append_counts
        )
        for exp, count in unmatched.items():
            all_unmatched[exp] += count

    # Summary
    total_appended = sum(append_counts.values())
    experiments_with_appends = sum(1 for c in append_counts.values() if c > 0)
    print(f"\nDone. Appended {total_appended} total rows across "
          f"{experiments_with_appends} experiment files.")

    if append_counts:
        print("\nPer-experiment append counts:")
        for exp in sorted(append_counts):
            print(f"  {merged_filename(exp)}: {append_counts[exp]} rows")

    if all_unmatched:
        total = sum(all_unmatched.values())
        print(f"\nWARNING: {total} rows across {len(all_unmatched)} experiments "
              f"in Undetermined files had no matching experiment file:")
        for exp in sorted(all_unmatched):
            print(f"  {exp} ({all_unmatched[exp]} rows)")


if __name__ == "__main__":
    main()


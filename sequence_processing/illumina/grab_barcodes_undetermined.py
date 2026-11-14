#!/usr/bin/env python3

import pandas as pd
import dnaio
import edlib
from rapidfuzz.distance import Levenshtein
import gzip
import subprocess
import shutil
import sys

# ---------- parameters ----------
MAX_IND_DIST = 1   # max edit distance for inline (r1/r2) indices
MAX_I57_DIST = 2   # max edit distance for i5/i7 indices

CONTAM_SEQ = 'AGGCCCG'
UPSTREAM = 'TGGCGTAATAGCGAAG'
DOWNSTREAM = 'CTGAATGGCGAATGGACGC'
UPSTREAM_L = len(UPSTREAM)


# ---------- cached fuzzy index matcher ----------

class IndexMatcher:
    """Caches fuzzy lookups so each distinct observed sequence is computed once."""

    def __init__(self, all_indices, max_dist):
        self.all_indices = list(all_indices)
        self.max_dist = max_dist
        self._cache = {}
        # Pre-seed exact matches
        for idx in self.all_indices:
            self._cache[idx] = idx

    def match(self, obs):
        """Return the unique matching index, or None if 0 or 2+ matches."""
        if obs in self._cache:
            return self._cache[obs]

        best = None
        for ref in self.all_indices:
            if Levenshtein.distance(obs, ref) <= self.max_dist:
                if best is not None:
                    # Ambiguous — two or more matches within max_dist
                    self._cache[obs] = None
                    return None
                best = ref

        self._cache[obs] = best
        return best


# ---------- sequence extraction ----------

def fuzzy_find(pattern, text, max_errors=1):
    result = edlib.align(
        pattern,
        text,
        mode="HW",       # infix/substring matching
        task="locations",
        k=max_errors,
    )
    return result["editDistance"], result["locations"]


def pull_seq_bw(seq, upstream, downstream, upstream_l):
    loc1 = seq.find(upstream)
    if loc1 == -1:
        ed, locs = fuzzy_find(upstream, seq, max_errors=1)
        if ed is None or ed > 1 or ed < 0:
            return None
        loc1 = locs[0][1] + 1
    else:
        loc1 += upstream_l

    loc2 = seq.find(downstream)
    if loc2 == -1:
        ed, locs = fuzzy_find(downstream, seq, max_errors=1)
        if ed is None or ed > 1 or ed < 0:
            return None
        loc2 = locs[0][0]

    return seq[loc1:loc2]


# ---------- I/O helpers ----------

def open_gzipped(path):
    """Open a gzipped file using pigz (multithreaded) if available, else gzip."""
    pigz = shutil.which('pigz')
    if pigz:
        proc = subprocess.Popen(
            [pigz, '-dc', path],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        return proc.stdout, proc
    else:
        return gzip.open(path, 'rb'), None


def extract_fastq_indices(fastq_gz_file, i7_matcher, i5_matcher):
    """
    Parse the undetermined FASTQ (gzipped) and extract matched i7/i5 per read.

    Returns dict: read_id -> (matched_i7, matched_i5)
    Only reads whose i7 AND i5 both unambiguously match are kept.
    We store the *matched* (canonical) indices so the main loop doesn't need to
    re-run fuzzy matching.
    """
    index_dict = {}

    stream, proc = open_gzipped(fastq_gz_file)
    try:
        with dnaio.open(stream, mode='r', fileformat='fastq') as reader:
            for record in reader:
                # dnaio exposes .name (up to first space) and .comment (rest of header)
                comment = record.comment or ''

                # Typical header: @READID 1:N:0:ACGT+TGCA
                # The indices field is the last colon-separated piece of the comment
                colon_pos = comment.rfind(':')
                if colon_pos == -1:
                    continue
                indices_field = comment[colon_pos + 1:]

                plus_pos = indices_field.find('+')
                if plus_pos == -1:
                    continue

                # Strip last 2 chars from each index (as in original script)
                raw_i7 = indices_field[:plus_pos - 2] if plus_pos >= 2 else indices_field[:plus_pos]
                raw_i5 = indices_field[plus_pos + 1:]
                raw_i5 = raw_i5[:-2] if len(raw_i5) >= 2 else raw_i5

                matched_i7 = i7_matcher.match(raw_i7)
                if matched_i7 is None:
                    continue
                matched_i5 = i5_matcher.match(raw_i5)
                if matched_i5 is None:
                    continue

                index_dict[record.name.split(' ')[0]] = (matched_i7, matched_i5)
    finally:
        if proc is not None:
            proc.stdout.close()
            proc.wait()

    return index_dict


# ---------- main ----------

def main():
    filename = sys.argv[1]
    filedir = sys.argv[2]
    fileloc = filedir + filename
    outdir = sys.argv[3]
    undetermined_r1 = sys.argv[4]
    plate = int(sys.argv[5])

    # --- metadata ---
    metadata = pd.read_csv('../seq_metadata.csv')

    i7_all = metadata['i7 seq rc'].unique().tolist()
    i5_all = metadata['i5 seq rc'].unique().tolist()

    # Build cached matchers
    i7_matcher = IndexMatcher(i7_all, MAX_I57_DIST)
    i5_matcher = IndexMatcher(i5_all, MAX_I57_DIST)

    # Filter to plate
    plate_meta = metadata[metadata['plate'] == plate]
    metadata_dict = plate_meta.set_index(['i7 seq rc', 'i5 seq rc']).T.to_dict()

    name = fileloc.split('/')[-1]

    # --- Step 1: extract i7/i5 from undetermined file (gzipped) ---
    print("Extracting indices from undetermined file...", flush=True)
    ind_dict = extract_fastq_indices(undetermined_r1, i7_matcher, i5_matcher)
    print(f"  Reads with valid i7/i5: {len(ind_dict):,}", flush=True)

    # --- Step 2: iterate main FASTQ (uncompressed) ---
    print("Processing main FASTQ...", flush=True)

    bc_counts = {}
    total_reads = 0
    rest_of_reads = 0
    reads_contam = 0

    with dnaio.open(fileloc, mode='r', fileformat='fastq') as reader:
        for record in reader:
            total_reads += 1

            read_id = record.name
            if read_id not in ind_dict:
                continue

            matched_i7, matched_i5 = ind_dict[read_id]

            # Plate lookup
            key_pair = (matched_i7, matched_i5)
            if key_pair not in metadata_dict:
                continue

            row = metadata_dict[key_pair]

            r1_ind = str(row['r1 index'])
            r2_ind = str(row['r2 index'])
            r1_ind_l = len(r1_ind)
            r2_ind_l = len(r2_ind)

            seq = record.sequence

            # Check inline indices
            r1_ind_obs = seq[8:r1_ind_l + 8]
            if r1_ind_obs != r1_ind and Levenshtein.distance(r1_ind_obs, r1_ind) > MAX_IND_DIST:
                continue

            r2_ind_obs = seq[-8 - r2_ind_l:-8]
            if r2_ind_obs != r2_ind and Levenshtein.distance(r2_ind_obs, r2_ind) > MAX_IND_DIST:
                continue

            # UMIs and barcode
            umi1 = seq[:8]
            umi2 = seq[-8:]
            bc = pull_seq_bw(seq, UPSTREAM, DOWNSTREAM, UPSTREAM_L)
            if bc is None:
                continue

            # Contamination check
            minibc = bc[:7]
            if Levenshtein.distance(minibc, CONTAM_SEQ) <= 1:
                reads_contam += 1
                continue

            new_key = (str(row['name']), umi1, umi2, bc)
            if new_key in bc_counts:
                bc_counts[new_key] += 1
            else:
                bc_counts[new_key] = 1

            rest_of_reads += 1

    # --- output ---
    print(f'Total Reads: {total_reads}')
    print(f'Passed Reads: {rest_of_reads}')
    print(f'Fraction Pass: {rest_of_reads / total_reads:.4f}' if total_reads > 0 else 'Fraction Pass: N/A')
    print(f'Contaminated Reads: {reads_contam}')

    outpath = f"{outdir}{name}_init_counts.csv"
    with open(outpath, "w") as f:
        f.write('experiment,umi1,umi2,bc,counts\n')
        for (exp, u1, u2, bc_val), count in bc_counts.items():
            f.write(f'{exp},{u1},{u2},{bc_val},{count}\n')

    print(f'Output written to {outpath}')


if __name__ == '__main__':
    main()


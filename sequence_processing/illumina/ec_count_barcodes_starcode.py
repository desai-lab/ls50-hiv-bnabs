"""
ec_count_barcodes_starcode.py
─────────────────────────────
Uses starcode (C, trie-based) for all Levenshtein-distance barcode operations,
and keeps lightweight Python Union-Find only for UMI deduplication (small groups).

Strategy overview:
─────────────────
Phase 1 — Barcode → super_bc assignment:
  a) Exact-match via dict lookup (handles the majority of 9.8M reads instantly)
  b) For remaining unmatched barcodes: write unique unmatched + all consensus BCs
     to a temp file, run starcode --connected-comp -d MAX_BC_DIST. Parse output
     to find which unmatched BCs cluster with a consensus BC → assign super_bc.
     Unmatched BCs that don't cluster with any consensus BC → "undetermined".

Phase 2 — UMI counting per super_bc:
  a) For the "undetermined" super_bc group: cluster raw barcodes with starcode
     -d MAX_BC_DIST, then within each BC cluster, count unique UMIs via Python
     Union-Find (groups are small).
  b) For all other super_bc groups: count unique UMIs via Python Union-Find.

Usage: python ec_count_barcodes_starcode.py <sample> <inputdir> <outputdir> <pacbiodir>

Requires: starcode binary available at $STARCODE or ~/software/starcode/starcode
"""

import pandas as pd
from rapidfuzz.distance import Levenshtein
import sys
import os
import subprocess
import tempfile
from collections import defaultdict

MAX_BC_DIST = 4
MAX_UMI_DIST = 2

# Locate starcode binary
STARCODE = os.environ.get('STARCODE',
           os.path.expanduser('~/software/starcode/starcode'))
THREADS = os.environ.get('SLURM_CPUS_PER_TASK', '8')


# ── Union-Find (only used for small per-group UMI dedup) ─────────────────────
class UnionFind:
    __slots__ = ('parent', 'rank')

    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x):
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return
        if self.rank[rx] < self.rank[ry]:
            rx, ry = ry, rx
        self.parent[ry] = rx
        if self.rank[rx] == self.rank[ry]:
            self.rank[rx] += 1

    def num_components(self):
        roots = set()
        for i in range(len(self.parent)):
            roots.add(self.find(i))
        return len(roots)


# ── Starcode helpers ─────────────────────────────────────────────────────────

def run_starcode(input_path, output_path, max_dist, threads=None):
    """Run starcode connected-components and return parsed clusters.

    Returns: list of (canonical_seq, count, [member_seqs])
    """
    if threads is None:
        threads = THREADS
    cmd = [
        STARCODE,
        '-d', str(max_dist),
        '-t', str(threads),
        '-i', input_path,
        '-o', output_path,
        '--connected-comp',
        '--print-clusters',
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"starcode stderr: {result.stderr}", file=sys.stderr)
        raise RuntimeError(f"starcode failed with exit code {result.returncode}")

    clusters = []
    with open(output_path, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            canonical = parts[0]
            count = int(parts[1])
            members = parts[2].split(',') if len(parts) > 2 else [canonical]
            clusters.append((canonical, count, members))
    return clusters


def starcode_cluster_to_seq_map(clusters):
    """Convert starcode cluster output to {sequence: cluster_id} mapping."""
    seq_to_cluster = {}
    for cid, (canonical, count, members) in enumerate(clusters):
        for seq in members:
            seq_to_cluster[seq] = cid
    return seq_to_cluster


# ── UMI dedup (unchanged from original — operates on small groups) ───────────

def ec_umis(u1_list, u2_list, max_dist):
    """Count unique UMIs using Union-Find."""
    n = len(u1_list)
    if n <= 1:
        return n

    uf = UnionFind(n)

    # Group by exact (umi1, umi2) pairs
    pair_groups = defaultdict(list)
    for i in range(n):
        pair_groups[(u1_list[i], u2_list[i])].append(i)

    unique_pairs = list(pair_groups.keys())

    # Merge identical pairs
    for indices in pair_groups.values():
        first = indices[0]
        for k in range(1, len(indices)):
            uf.union(first, indices[k])

    num_unique = len(unique_pairs)
    if num_unique <= 1:
        return 1

    # Pre-compute lengths for length-filter
    u1_lens = [len(p[0]) for p in unique_pairs]
    u2_lens = [len(p[1]) for p in unique_pairs]
    reps = [pair_groups[p][0] for p in unique_pairs]

    for i in range(num_unique):
        ri = reps[i]
        u1i, u2i = unique_pairs[i]
        l1i, l2i = u1_lens[i], u2_lens[i]
        for j in range(i + 1, num_unique):
            if uf.find(ri) == uf.find(reps[j]):
                continue
            if abs(l1i - u1_lens[j]) > max_dist:
                continue
            if abs(l2i - u2_lens[j]) > max_dist:
                continue
            u1j, u2j = unique_pairs[j]
            if Levenshtein.distance(u1i, u1j, score_cutoff=max_dist) <= max_dist:
                if Levenshtein.distance(u2i, u2j, score_cutoff=max_dist) <= max_dist:
                    uf.union(ri, reps[j])

    return uf.num_components()


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    fileloc0 = sys.argv[1]
    inputdir = sys.argv[2]
    outputdir = sys.argv[3]
    pacbiodir = sys.argv[4]

    # Verify starcode exists
    if not os.path.isfile(STARCODE):
        print(f"ERROR: starcode not found at {STARCODE}", file=sys.stderr)
        print("Set STARCODE env var or install to ~/software/starcode/starcode",
              file=sys.stderr)
        sys.exit(1)


    name0 = fileloc0.split('/')[-1].split('_init_')[0]
    name = '_'.join(name0.split('0_JAsce_')[1].split('_')[:-2])

    # Create a temp directory for intermediate files (cleaned up at end)
    tmpdir = tempfile.mkdtemp(prefix='ec_starcode_'+name+'_')
    print(f"Temp directory: {tmpdir}")

    try:
        _run(fileloc0, inputdir, outputdir, pacbiodir, tmpdir)
    finally:
        # Clean up all intermediate files
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)
        print("Cleaned up temp files.")


def _run(fileloc0, inputdir, outputdir, pacbiodir, tmpdir):
    # ── Load reference data ──────────────────────────────────────────────
    print("Loading reference data...")
    consensus_df = pd.read_csv(pacbiodir + 'consensus_seqs.csv')
    consensus_bcs = list(consensus_df['bc_cons'])
    consensus = consensus_df[['super_bc', 'bc_cons']].set_index('bc_cons').to_dict()['super_bc']
    consensus_set = set(consensus_bcs)

    all_barcodes = (pd.read_csv(pacbiodir + 'barcodes_pure.csv')[['super_bc', 'bc']]
                    .set_index('bc').to_dict()['super_bc'])

    name0 = fileloc0.split('/')[-1].split('_init_')[0]
    name = '_'.join(name0.split('0_JAsce_')[1].split('_')[:-2])

    print(f"Loading init_counts for {name}...")
    init_counts = pd.read_csv(inputdir + name + '_init_counts_merged.csv')
    n_reads = len(init_counts)
    print(f"  {n_reads:,} reads loaded")

    # ═════════════════════════════════════════════════════════════════════
    # PHASE 1: Assign super_bc to each read
    # ═════════════════════════════════════════════════════════════════════

    print("Phase 1: Barcode assignment...")

    # Step 1a: Exact matches (vectorized)
    bc_series = init_counts['bc']
    super_from_consensus = bc_series.map(consensus)
    super_from_all = bc_series.map(all_barcodes)
    superbarcodes = super_from_consensus.fillna(super_from_all)

    n_exact = superbarcodes.notna().sum()
    print(f"  Exact matches: {n_exact:,} / {n_reads:,} ({100*n_exact/n_reads:.1f}%)")

    # Step 1b: Fuzzy match remaining via starcode
    mask_unmatched = superbarcodes.isna()
    n_unmatched = mask_unmatched.sum()

    if n_unmatched > 0:
        unmatched_unique = bc_series[mask_unmatched].unique()
        # Filter out any that happen to match consensus or all_barcodes
        # (shouldn't happen, but be safe)
        unmatched_unique = [b for b in unmatched_unique
                           if b not in consensus and b not in all_barcodes]
        print(f"  Unique unmatched barcodes: {len(unmatched_unique):,}")

        # Write consensus BCs + unmatched BCs to a file for starcode.
        # We tag consensus BCs so we can identify them in the output.
        # Strategy: write each consensus BC once, then each unmatched BC once.
        # After starcode clusters them, any unmatched BC that clusters with
        # a consensus BC gets assigned that consensus BC's super_bc.
        sc_input = os.path.join(tmpdir, 'bc_fuzzy_input.txt')
        sc_output = os.path.join(tmpdir, 'bc_fuzzy_output.txt')

        print("  Writing starcode input for fuzzy matching...")
        with open(sc_input, 'w') as f:
            for bc in consensus_bcs:
                f.write(bc + '\n')
            for bc in unmatched_unique:
                f.write(bc + '\n')

        n_consensus = len(consensus_bcs)
        print(f"  Running starcode (d={MAX_BC_DIST}) on {n_consensus + len(unmatched_unique):,} sequences...")
        clusters = run_starcode(sc_input, sc_output, MAX_BC_DIST)

        # Build mapping: unmatched_bc -> super_bc
        ec_map = {}
        n_ec = 0
        n_na = 0
        for canonical, count, members in clusters:
            # Find which consensus BCs are in this cluster
            consensus_in_cluster = [m for m in members if m in consensus_set]
            non_consensus = [m for m in members if m not in consensus_set]

            if consensus_in_cluster and non_consensus:
                # Assign all non-consensus members to the super_bc of the
                # first consensus BC found (they're all in the same connected
                # component at distance ≤ MAX_BC_DIST from at least one
                # consensus BC).
                # If multiple consensus BCs are in the same cluster, use the
                # canonical (highest count) consensus BC's super_bc, which is
                # the most representative.
                # Pick the consensus BC with the highest implicit count
                # (starcode canonical is the highest-count member)
                assigned_sbc = None
                if canonical in consensus_set:
                    assigned_sbc = consensus[canonical]
                else:
                    # Pick the first consensus BC found
                    assigned_sbc = consensus[consensus_in_cluster[0]]

                for bc in non_consensus:
                    ec_map[bc] = assigned_sbc
                    n_ec += 1
            elif not consensus_in_cluster:
                # Cluster of only unmatched BCs — these are undetermined
                for bc in non_consensus:
                    ec_map[bc] = 'undetermined'
                    n_na += 1
            # else: cluster of only consensus BCs, no unmatched to assign

        # Also mark any unmatched BCs not found in starcode output as undetermined
        for bc in unmatched_unique:
            if bc not in ec_map:
                ec_map[bc] = 'undetermined'
                n_na += 1

        print(f"  Fuzzy matched: {n_ec:,}, undetermined: {n_na:,}")

        # Apply
        fuzzy_results = bc_series[mask_unmatched].map(ec_map)
        superbarcodes.loc[mask_unmatched] = fuzzy_results.values

    # Fill any remaining NaN as undetermined (shouldn't happen, but be safe)
    superbarcodes = superbarcodes.fillna('undetermined')

    # ═════════════════════════════════════════════════════════════════════
    # PHASE 2: UMI counting per super_bc
    # ═════════════════════════════════════════════════════════════════════

    print("Phase 2: UMI deduplication...")
    init_counts = init_counts.copy()
    init_counts['super_bc'] = superbarcodes.values

    grouped = init_counts.groupby('super_bc')
    n_groups = grouped.ngroups
    print(f"  {n_groups:,} super_bc groups")

    sbcs = []
    umi_counts = []

    for sbc, group in grouped:
        n_g = len(group)
        sbcs.append(sbc)

        if n_g == 1:
            umi_counts.append(1)
            continue

        if sbc == 'undetermined':
            # Cluster undetermined barcodes with starcode, then UMI dedup
            # within each barcode cluster
            bc_list = group['bc'].tolist()
            u1_list = group['umi1'].tolist()
            u2_list = group['umi2'].tolist()

            # Write undetermined BCs to file
            sc_undet_input = os.path.join(tmpdir, 'undet_bc_input.txt')
            sc_undet_output = os.path.join(tmpdir, 'undet_bc_output.txt')

            # Starcode needs one sequence per line; duplicates get counted.
            # We need to track which input row maps to which cluster.
            # Use --seq-id to get 1-based line numbers.
            with open(sc_undet_input, 'w') as f:
                for bc in bc_list:
                    f.write(bc + '\n')

            print(f"  Clustering {n_g:,} undetermined reads (d={MAX_BC_DIST})...")
            clusters = run_starcode(sc_undet_input, sc_undet_output, MAX_BC_DIST)

            # Build bc -> cluster_id mapping from starcode output
            bc_to_cluster = starcode_cluster_to_seq_map(clusters)

            # Assign cluster labels to each read
            cluster_labels = []
            next_fallback = len(clusters)
            for bc in bc_list:
                if bc in bc_to_cluster:
                    cluster_labels.append(bc_to_cluster[bc])
                else:
                    # Shouldn't happen, but assign unique cluster
                    cluster_labels.append(next_fallback)
                    next_fallback += 1

            # Group by cluster, then UMI dedup within each
            cluster_groups = defaultdict(list)
            for i, cid in enumerate(cluster_labels):
                cluster_groups[cid].append(i)

            all_counts = 0
            for indices in cluster_groups.values():
                u1 = [u1_list[i] for i in indices]
                u2 = [u2_list[i] for i in indices]
                all_counts += ec_umis(u1, u2, MAX_UMI_DIST)

            umi_counts.append(all_counts)
        else:
            c_i = ec_umis(group['umi1'].tolist(), group['umi2'].tolist(),
                         MAX_UMI_DIST)
            umi_counts.append(c_i)

    # ═════════════════════════════════════════════════════════════════════
    # OUTPUT
    # ═════════════════════════════════════════════════════════════════════

    outpath = outputdir + name + '_umi_counts2.csv'
    pd.DataFrame({'super_bc': sbcs, 'umi_counts': umi_counts}).to_csv(
        outpath, index=False)
    print(f"Done. Wrote {len(sbcs):,} rows to {outpath}")


if __name__ == '__main__':
    main()


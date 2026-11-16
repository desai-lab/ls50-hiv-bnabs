# Sequencing analysis code



- **`pacbio/`** — builds the *barcode → scFv genotype* lookup table from long-read
  PacBio HiFi amplicon sequencing of the pooled, downsampled library.
- **`illumina/`** — counts barcodes (UMI-deduplicated) in each FACS bin × antigen
  concentration × replicate sample from short-read Illumina sequencing, using the
  PacBio table to assign reads to variants.

Both pipelines are written as SLURM batch scripts (`.sbatch`) wrapping Python/shell
workers, and were run on the Harvard FASRC cluster (`desai` partition). Each
directory has a `pipeline_order.txt` giving the submission order.

---

## Requirements

- **Python** with `numpy`, `pandas`, `biopython`, `edlib`, `rapidfuzz`, `dnaio`,
  `networkx`, `numba`, `spoa` (conda envs `barcode_env` / `ec_env` on the cluster).
- **[starcode](https://github.com/gui11aume/starcode)** — expected at
  `$STARCODE` or `~/software/starcode/starcode`.
- **NGmerge** (read merging), **bam2fastq** from PacBio `pbtk` (BAM → FASTQ).
- `pigz` is used when available; otherwise `gzip`.

Paths to raw data, scratch, and output directories are hard-coded at the top of
each `.sbatch` file and in `make_count_tables.py`.

---

## `pacbio/` — barcode-to-genotype association

Run in the order listed in `pacbio/pipeline_order.txt`.

| Step | Script | What it does |
|---|---|---|
| 1 | `bamfastq_conv.sbatch` | `bam2fastq` on the HiFi BAM → `pacbio_reads.fastq`. |
| 2 | `init_filter.sbatch` → `init_filter.sh` | Orientation filter: an awk pass keeps reads containing both constant sequences of the forward pair, reverse-complements reads matching the reverse pair (quality string reversed with it), and discards reads matching both or neither. Output `../pacbio-filtered.fastq`. |
| 3 | `grab_barcodes.sbatch` → `grab_barcodes.py` | Extracts the DNA barcode, VH and VL regions plus their Phred strings from each read, using flanking constant sequences; exact match preferred, else Edlib infix alignment (`mode="HW"`) allowing ≤1 edit per flank. Reads missing any region are dropped. Output `../grab_barcodes.csv`. |
| 4 | `extract_barcodes.sbatch` → `extract_barcodes.py` | Restricts barcodes to 40–70 nt (`../grab_barcodes_filtered.csv`), writes the unique barcode list `../barcodes.txt` for starcode and a barcode-length histogram `../length_counts.csv`. |
| 5 | `run_starcode.sbatch` | Clusters unique barcodes with starcode in connected-component mode, max Levenshtein distance 6, `--print-clusters` → `../clustered_barcodes.txt`. |
| 6 | `merge_clusters.sbatch` → `merge_clusters.py` | Parses the starcode output into `super_bc` cluster IDs (`BC_1`, `BC_2`, …), joins them back onto every read (`../merged_clusters.csv`), and writes a randomly shuffled list of super-barcodes (`../shuffled_sbcs.csv`) so the array job below gets even chunks. |
| 7 | `submit_ec_pipeline.sh` | Launches the error-correction array job and schedules the merge job with `--dependency=afterok`. |
| 7a | `errorcorrection_array.sbatch` → `errorcorrection_parallel.py` | Performs error correction and computes consensus sequences for both barcodes and scFv sequences, split over 50 array tasks. Writes `chunks/chunk_<task>.csv` and `chunks/bcs_<task>.csv`. |
| 7b | `errorcorrection_merge.sbatch` → `errorcorrection_merge.py` | Concatenates the chunk files, reports how many super-barcodes passed, and writes the pipeline's final files: `../superbarcodes_all.csv`, `../superbarcodes_pure.csv` (purity = 1), `../consensus_seqs.csv` (consensus VH/VL/barcode/mini-barcode per super-barcode) and `../barcodes_pure.csv` (observed barcode → super-barcode map). |

---

## `illumina/` — Tite-Seq barcode counting

Run in the order listed in `illumina/pipeline_order.txt`. There are 188 demultiplexed
samples (antigen × replicate × concentration × FACS bin, plus expression-sort bins),
spread over two NovaSeq lanes.

| Step | Script | What it does |
|---|---|---|
| 1 | `run_ngmerge.sbatch` | NGmerge on each R1/R2 pair (array of 190; minimum overlap 80 nt, expected read length 150, gzipped output) into `merged_fastq/`. Sample names come from `sample_list.txt`, built by the first lines of `pipeline_order.txt`. |
| 2 | `gunzip_all.sbatch` | Decompresses the merged FASTQs in place. |
| 3 | `grab_barcodes.sbatch` → `grab_barcodes.py` | Per sample (array of 188, using `sample_list2.txt`, i.e. `sample_list.txt` minus the Undetermined files). Looks up the sample's expected inline indices in `../seq_metadata.csv` and requires both to match within Levenshtein distance 1; takes the first and last 8 nt of the merged read as the paired UMIs; extracts the barcode between its flanking constant sequences (exact, else Edlib ≤1 edit per flank); drops reads whose first 7 barcode bases are within 1 edit of the known contaminant 7-mer. Writes `(umi1, umi2, bc) → read count` as `{sample}_init_counts.csv` and prints pass/contaminant/wrong-index tallies. |
| 4 | `grab_barcodes_undetermined{1,2}.sbatch` → `grab_barcodes_undetermined.py` | Recovers reads that failed Illumina demultiplexing, one job per lane/plate. Parses i7/i5 from the original FASTQ headers (dropping the terminal 2 nt of each), matches them to the expected indices within 2 edits and only if the match is unique, maps the index pair to a sample on that plate via `../seq_metadata.csv`, then applies the same inline-index, UMI, barcode and contamination handling as step 3. Output carries an extra `experiment` column. Uses `dnaio` for parsing and `pigz` for decompression. |
| 5 | `merge_undetermined.sbatch` → `merge_undetermined.py` | Streams the two Undetermined tables row by row (never loading them fully) and appends each row to its sample's file, producing `{sample}_init_counts_merged.csv` in `merged_barcodes/`; reports per-sample append counts and any rows with no matching sample. |
| 6 | `ec_count_barcodes_starcode.sbatch` → `ec_count_barcodes_starcode.py` | The main counting step, one array task per sample. **Phase 1** assigns each read a `super_bc`: exact lookup against the PacBio consensus barcodes and then the observed barcodes; remaining unmatched barcodes are clustered together with the consensus barcodes by starcode (connected components, distance 4) and inherit the super-barcode of a consensus member of their component (canonical member preferred), or are labelled `undetermined`. **Phase 2** deduplicates the paired UMIs within each super-barcode group with a Union-Find over pairs linked when *both* UMIs are within 2 edits, so each connected component contributes one molecular count regardless of read multiplicity; the `undetermined` group is first sub-clustered by barcode (distance 4) before UMI dedup. Output `{sample}_umi_counts2.csv` (`super_bc`, `umi_counts`). |
| 7 | `make_count_tables.sbatch` → `make_count_tables.py` | Outer-joins the per-sample UMI counts on `super_bc` within each experiment group of `../seq_metadata.csv`, filling missing entries with 0, and writes one barcode × sample matrix per experiment: `{experiment}_counts.csv`. Falls back to a `umi_counts_subsample/` directory for samples that had to be rerun on subsampled input. |

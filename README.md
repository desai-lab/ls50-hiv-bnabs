# Combinatorial affinity landscapes in nine HIV-1 bnAbs

Code and processed data for "Epistasis generates antigen-specific tradeoffs and reshapes accessible maturation paths in HIV-1 bnAbs."

## Layout

```
data/                measurements, site metadata, and fitted landscapes (8.2 MB)
notebooks/           one notebook per figure or table
sequence_processing/ scripts to process PacBio and Illumina sequencing reads
figures/             generated on run
tables/              generated on run
```

&nbsp;

```
git clone <repo-url>
cd hiv-bnabs
pip install -r requirements.txt
jupyter lab notebooks/
```

Data

Files are flat and named according to the following antibody stems:

```
1NC9  CH235_12  N49P7  NIH4546  VRC01  VRC07  VRC13  VRC16  VRCPG20
```

The two antigens these antibodies are measured against are `CH505TF` and `SF162`, giving 54 files.

| File | Contents |
| :---- | :---- |
| `{stem}_kds.csv` | Tite-seq measurements, one row per genotype and replicate. Phenotype is `nlog10_Kd`; `site_0` through `site_9` hold G or M. |
| `{stem}_meta.csv` | The ten varied sites, with Kabat and IMGT numbering, region, and germline and mature residues. |
| `seq_metadata.csv` | Metadata for Illumina sequencing runs. Experiments LA, LB, LC, and LD correspond to SF162 (replicates 1 and 2) and CH505TF (replicates 1 and 2), respectively. |
| `{stem}_{antigen}_BEST_df_params.csv` | Coefficients from the LOGO-selected Tobit fit, with standard errors. |
| `{stem}_{antigen}_BEST_df_all.csv` | Predicted affinity for all 1024 genotypes, censored and latent. |

&nbsp;

Affinities are left-censored at the assay detection floor of 6.0 in units of \-log10 KD, with KD in molar. Larger values mean tighter binding. Genotypes pinned to the floor are measured lower bounds, not resolved affinities, and several notebooks jitter them into a band below the threshold so the pile-up stays readable.

## Notebooks

| Notebook | Produces |
| :---- | :---- |
| `[Fig 1_2] Summary Statistics` | Figures 1B, 2, S1 |
| `[Fig 3A] Heatmaps` | Figure 3A |
| `[Fig 3B] Variance Partitioning` | Figure 3B |
| `[Fig 4_5] 1NC9 + CH235.12 Focus Panels` | Figures 4A-E, 5A-E |
| `[Fig 6A_C] Antigen-Antigen Scatters` | Figures 6A, S20, S21, S22 |
| `[Fig 6B_D] Violins and Deltas` | Figure 6D and supplement violin variants |
| `[Fig 7] Simulations` | Figures 7A, 7B |
| `[Methods] Tobit Model Inference` | The fits in `data/`, reimplemented so the method reads in one file |
| `[Supplement] 1NC9 K19T` | Figure S8 |
| `[Supplement] Coefficients vs. Distance` | Figures S5A, S5B, Table S1 |
| `[Supplement] Exploratory Data Analysis` | Figure S2 and the genotype-count and censoring numbers |
| `[Supplement] Maturation Strip Plots` | Figures S6A, S6B |
| `[Supplement] NIH45-46 L45P` | Table S2 |
| `[Supplement] NIH45-46 Matched Pair Scatters` | Figures S17, S18, S19 |
| `[Supplement] Observed vs. Predicted` | Figures S3, S4 |
| `[Supplement] Site 54 Across Backgrounds` | Figure S12 |
| `[Supplement] Site Frequency By Affinity` | Figure S9 |
| `[Supplement] Site Statistics` | Effect sizes and significance behind Figures 1B and S1 |
| `[Supplement] VRC07 Violins` | Figures 6B, S13 |

&nbsp;

## Sequence processing

Scripts to process both PacBio reads (to associate scFv sequence to variant barcodes) as well as Illumina reads (to estimate frequencies of scFv variants across sorted bins in Tite-Seq experiments). See `sequence_processing/README.md` for further details.

## Model fitting

`[Methods] Tobit Model Inference` is a self-contained implementation of the censored regression described in the Methods: genotype encoding and left-censoring, the interaction-only polynomial basis, the Tobit log-likelihood with its analytic gradient and elastic-net penalty, penalized maximum likelihood, standard errors from the inverse observed information, prediction over all 1024 genotypes, leave-one-genotype-out cross-validation, and model selection across the order by penalty grid.
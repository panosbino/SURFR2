# SURFR2

Reference-free discovery of small-RNA sequences **specific to one condition**. SURFR2 counts k-mers per sample, normalises for sequencing depth, and compares any two conditions defined in a config file. It is built for everyday experiments, including very small designs:

| Design | What SURFR2 does |
|---|---|
| 1 vs 1 | **Descriptive mode.** Specificity filters only. Results are flagged as exploratory, because nothing can be tested without replicates. |
| ≥ 2 vs ≥ 2 | Specificity filters **plus** an edgeR quasi-likelihood test with BH FDR. |
| Paired samples or batches | Add a `block` column; it enters the test's design. |
| Several independent datasets | Add a `cohort` column; results can be required to replicate in each. |

## Quick start

```bash
cp config/config.example.yaml my.yaml          # edit: project, samplesheet, outdir, comparison, execution
cp config/samples.example.tsv my_samples.tsv   # one row per sample
bash run_SURFR2.sh -c my.yaml --dry-run        # validate; prints design, mode and thresholds per group
bash run_SURFR2.sh -c my.yaml                  # run (SLURM or local, set in the config)
```

## Inputs

**Samplesheet** (TSV, columns in any order):

| Column | Required | Meaning |
|---|---|---|
| `sample_id` | yes | unique; letters, digits, `.` `_` `-` |
| `condition` | yes | e.g. `treated` / `untreated`; conditions not named in the config are ignored |
| `file_type` | yes | `fastq` or `bam` |
| `path` | yes | relative paths resolve against the samplesheet's directory |
| `block` | no | pairing or batch (e.g. donor). Must be set for all or none of the compared samples. Ignored in descriptive mode. |
| `cohort` | no | independent datasets, compared separately. Default: all samples form one cohort. |

**Config**: see `config/config.example.yaml`, where every option is documented. The comparison is `comparison.case` vs `comparison.control`.

The validator checks the design before anything runs, and reports for each cohort:
- the group sizes
- whether it runs in *test* or *descriptive* mode, and why
- the resolved thresholds, e.g. "needs ≥ 3 case, ≤ 0 control detected"

It rejects designs that cannot work, such as a `block` that is confounded with condition when `test: edger` is forced.

## What a "case-specific" k-mer is

A k-mer passes when all of the following hold in each compared cohort:

1. **Detected in cases.** CPM ≥ `min_cpm_case` in at least `min_case_samples_detected` case samples (default: **all**).
2. **Abundant enough.** Mean case CPM ≥ `min_mean_cpm_case`.
3. **Absent from controls.** CPM ≥ `max_cpm_control` in at most `max_control_samples_detected` control samples (default: **none**).
4. **Large fold change.** log2((mean case + pc) / (mean control + pc)) ≥ `min_log2fc`.
5. **Significant** *(test mode only)*: edgeR FDR ≤ `statistics.fdr`.
6. **Absent from external controls**, if any are configured.

Sample-count thresholds accept `all`, an integer, or a fraction (e.g. `0.5`). Overlapping passing k-mers are then merged into sequences with dekupl-mergeTags.

**Specificity filters and the test answer different questions,** so SURFR2 requires both.
- A k-mer can be highly significant yet present in every control, e.g. 8× up.
- A k-mer can be perfectly specific yet present in too few samples for significance.

The regression tests include both cases.

## Statistics

- **Test.** edgeR quasi-likelihood F-test (`glmQLFit(robust = TRUE)` + `glmQLFTest`) of case vs control, with `~ block + condition` when blocks are given. Empirical-Bayes sharing of variability across k-mers makes 2–3 replicates per group workable.
- **Normalisation.** Tests use SURFR2's own effective library sizes (`norm.factors = 1`). TMM computed on the candidate set would be biased, because candidates are enriched for differences.
- **Condition-blind candidate selection.** The pre-filter counts samples regardless of condition. Selecting on case counts and then testing the same data would bias the p-values (Bourgon et al. 2010, *PNAS*).
- **One test per distinct count profile.** Overlapping k-mers of one molecule have identical counts; a 22-nt read yields six 17-mers. Tested separately, they act as pseudo-replicates: in testing they inflated edgeR's prior degrees of freedom about 6-fold, and they multiply the number of BH tests. SURFR2 tests each distinct profile once and maps the result back, which is lossless for p-values.
- **FDR** is controlled at the k-mer (profile) level, per cohort. Each merged sequence also gets an ACAT-combined p-value of its k-mers (Cauchy combination; Liu & Xie 2020, *JASA*), which remains valid under their strong correlation.

**Limits to keep in mind:**
- **n = 1 per group cannot be tested.** Descriptive mode says so in `run_summary.txt`, in the hover text and in the outputs.
- **Heterogeneous markers have low power.** For a k-mer present in only some case samples, the evidence can be weak even when it is perfectly specific. For example, 5 of 7 cases vs 0 of 6 controls gives p ≈ 0.01 by edgeR and by an exact test on detection alone. That k-mer will not pass at 5% FDR alongside hundreds of tests. If such markers matter, relax `min_case_samples_detected` and inspect the descriptive results.
- **Thresholds are defaults, not recommendations.** Choose CPM and fold-change thresholds for your library depth and question, and report them.

## Steps

| Step | Script | Unit |
|---|---|---|
| sample | `bin/surfr2_sample.sh` | one task per sample: BAM→FASTQ (drops secondary alignments, restores read orientation), miRTrace QC, KMC (`-ci1`, stranded) |
| matrix | `bin/surfr2_matrix.sh` | condition-blind candidate k-mers via KMC set operations; per-sample candidate counts; reference k-mers for `median_ratio` |
| filter | `bin/surfr2_filter.R` | normalisation, filters, edgeR, mergeTags, QC plots |
| report | `bin/surfr2_report.R` | scatterplots per cohort: PDF and interactive HTML |

Re-launching is safe. Each sample and the matrix step store a fingerprint of their inputs and parameters, and only changed work is redone. Use `--from filter` after changing filter or statistics settings, or `--from report` to redraw the plots.

## Outputs (`<outdir>/results/`)

| File | Content |
|---|---|
| `case_specific_sequences.tsv` | merged sequences: the main result. In test mode it adds `<cohort>_acat_pvalue` and `<cohort>_min_kmer_fdr`. |
| `case_specific_kmers.tsv` | per-cohort statistics of passing k-mers: detection counts, mean CPM, log2FC, p-value, FDR |
| `case_specific_kmers_{counts,cpm}.tsv` | k-mer × sample matrices |
| `candidate_kmer_stats.tsv.gz` | the same statistics plus a pass flag for every filter, for **all** candidates, to audit why a k-mer failed |
| `normalisation.tsv` | library sizes, size factors, effective library sizes |
| `run_summary.txt` | mode per cohort, test used, counts, warnings (**read this first**) |
| `plots/`, `interactive/` | QC plots, `scatter_<cohort>.pdf`, interactive `scatter_<cohort>.html` |

## Interactive results

**HTML report** (`results/interactive/scatter_<cohort>.html`). Open the file in a browser; no server is needed and it works offline.
- **Axes:** mean CPM in control (x) vs case (y), log10 with the configured pseudocount.
- **Background:** all candidate k-mers as an exact binned density (not hoverable).
- **Gold points:** case-specific k-mers. Hovering shows the k-mer, its merged sequence, and raw counts, mean CPM and detection rate for every cohort/condition. With ≤12 samples in total, it also lists each sample's count.
- **Dashed lines:** the `min_mean_cpm_case` and `min_log2fc` cut-offs. The diagonal is the exact decision boundary, because log2FC uses the same pseudocount as the axes.

**Shiny explorer** (`app/`) adds two things the HTML cannot do. You can click a k-mer to see per-sample counts and CPM, as a plot and a downloadable CSV. You can also search a known sequence (DNA or RNA, e.g. a mature miRNA) to highlight its k-mers. It needs only the `results/` folder, so the simplest route is to copy that folder to your laptop:

```bash
# laptop R needs: install.packages(c("shiny", "plotly", "dplyr", "readr", "ggplot2", "jsonlite"))
SURFR2_RESULTS=/path/to/results Rscript -e 'shiny::runApp("SURFR2/app", launch.browser = TRUE)'
```

To run it on Dardel instead, use an SSH tunnel:

```bash
# Dardel (login node is fine; the app is light)
SURFR2_RESULTS=<outdir>/results singularity exec -B /cfs/klemming <sif> \
  Rscript -e 'shiny::runApp("SURFR2/app", port = 8787, host = "127.0.0.1")'
# laptop
ssh -N -L 8787:127.0.0.1:8787 <user>@<the login node you started it on>   # then open http://localhost:8787
```

## Container change

Three additions to the Dockerfile:
1. In step 1, add `python3-yaml \` and `pandoc \` to the apt list. PyYAML is needed by the config validator, and pandoc makes the HTML reports single self-contained files.
2. Add `"plotly", "htmlwidgets", "shiny"` to the CRAN package vector in `install_r_packages.R`, and install edgeR from Bioconductor: `BiocManager::install("edgeR")`.
3. Extend the smoke tests:

```
python3 -c "import yaml" && \
pandoc --version | head -1 && \
Rscript --vanilla -e "for(p in c('edgeR','plotly','htmlwidgets','shiny','jsonlite','readr','dplyr')){library(p,character.only=TRUE,lib.loc=c('/opt/R/library',.libPaths()));cat(p,'OK\n')}" && \
```

Without pandoc, the report still runs, but each HTML file needs its `*_files/` folder next to it. Without plotly, the interactive plots are skipped with a warning and the PDFs are still produced.

## Validation

`bash tests/run_test.sh` generates synthetic data with planted truth, runs the full pipeline, and checks the results. Counts include biological noise (log-normal, CV ≈ 30%) and Poisson sampling, so the test is exercised realistically. There are five scenarios:

| Scenario | Design | Checks |
|---|---|---|
| `small` | 3 vs 3, paired | recovers the true sequence; rejects one significant but not specific, one specific but too low, and one in only 2 of 3 cases |
| `single` | 1 vs 1 | descriptive mode; flagged exploratory; no p-values |
| `null` | 3 vs 3, no true differences | no FDR discoveries; p-values not inflated (5.7% below 0.05) |
| `multi` | 2 cohorts + external controls, test on | replication and external-control logic; the underpowered marker is correctly not called |
| `multi_descriptive` | same data, test off | the filter-only path recovers it |

Also tested:
- exact read counts from BAMs, including removal of secondary alignments and restoration of reverse-strand reads
- resume and invalidation after changes
- SLURM job chains for every `--from` entry point, against a mocked scheduler
- the tree reduction used for candidate selection, against a brute-force union for 1–9 inputs

Not yet tested on real data or a real cluster. For a first real run, include a known positive control (e.g. a tissue-specific miRNA) and check that it is recovered.

## Fixes relative to SURFR1

These were verified against the KMC 3.2.4 and dekupl-mergeTags sources.

SURFR2 also corrects several issues found in SURFR1's scripts.

1. **mergeTags `-n` means unstranded merging** (`case 'n': stranded = 0`). SURFR1 used it on stranded libraries. SURFR2 passes `-n` only when `kmer.canonical: true`.
2. **mergeTags reads DE-kupl's column layout.** Value column 1 is read as a p-value, and the lowest becomes the contig representative. Value column 4 is read as log2FC, and it is read unconditionally, so a table with fewer than 4 value columns is read out of bounds. SURFR1's first value column was the raw TCGA cancer count, so each contig was represented by its least abundant k-mer. SURFR2 writes `rank, mean case CPM, mean control CPM, log2FC`.
3. **Enrichment ignored group size and depth.** `sample_ratio` was computed but never used.
4. **Wrong flag comments.** KMC `-b` is "no canonical form", not RAM-only (the flag itself was right). `dump -s` means sorted. `-cs4294967296` exceeds the 32-bit counter maximum of 2³²−1.
5. **Module loading inside the container.** `ml PDC` ran inside `singularity exec`, where Lmod is normally unavailable. SURFR2 loads modules on the host before `singularity exec`.
6. **Silent miRTrace failures.** miRTrace 1.0.1 can exit 0 after aborting, e.g. when PHRED auto-detection fails. SURFR2 checks both the output and the log. Set `qc.phred_offset` if auto-detection fails.


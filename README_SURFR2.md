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
| `tool_versions.txt`, `sessionInfo.txt` | exact tools (path, version, checksum), modules and R packages this run used |
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

## Environments: container or modules

Every step runs in one of two environments, set by `execution.environment`. SURFR2's own code (`bin/`, `app/`) always runs from this repository, so **changing SURFR2 code never requires rebuilding anything**. Only adding or upgrading a *tool or R package* does, and only in the container environment.

| | `container` | `modules` |
|---|---|---|
| Tools from | the SURFR2 image | environment modules, `path_prepend`, `tools.*` paths |
| R packages from | the image | `r_libs` (a personal library) |
| Adding a tool or package | rebuild the image | install it, or load a module |
| Reproducibility | fixed by the image | recorded per run in `tool_versions.txt` |
| Use for | production runs, publication | development |

The tool versions are defined in two scripts, `container/install_tools.sh` and `container/install_r_packages.R`. The Dockerfile runs them, and so does a development setup, so both environments get the same binaries.

**How a step runs.** The launcher writes `<run dir>/env.sh` (module loads, PATH additions, R library) and one job script per step in `<run dir>/jobs/`. Each job script is a login shell (`#!/bin/bash -l`), which initialises the module system in batch jobs, then sources `env.sh` and runs the step, inside the image in the container environment. A login shell resets PATH from the system profile, so tools must be declared through modules, `path_prepend` or `tools.*`; nothing is inherited from the shell you launched from. The job scripts can be inspected and resubmitted by hand.

**Tool check.** Before anything is submitted, the launcher runs `jobs/check.sh` in the same environment. It fails, with instructions, if a required tool or R package is missing; edgeR is required only when a cohort will be tested. It writes `tool_versions.txt` with each tool's resolved path, version and SHA-256, the loaded modules, R version and package versions. The filter step copies it, with R's `sessionInfo()`, into `results/`. In the modules environment this record is the only reliable account of what ran, because module defaults change; the validator warns about modules given without a version.

### Building the container

From the repository root on a machine with Docker:

```bash
docker buildx build --platform linux/amd64 -f container/Dockerfile -t surfr2:<version> --load container/
docker save surfr2:<version> -o surfr2_<version>.tar
```

On the cluster: `singularity build surfr2_<version>.sif docker-archive://surfr2_<version>.tar`, then set `execution.container` and `execution.modules: [PDC, singularity]`.

The image is based on `rocker/r-ver:4.4.3`, which pins R and installs CRAN packages from a dated snapshot. It also contains samtools 1.23.1, KMC 3.2.4 (official static binaries), miRTrace 1.0.1 on OpenJDK 17, mergeTags at a pinned commit, pigz, pandoc and PyYAML. Every download is checksum-verified, and the build fails if a smoke test fails.

### Development environment on Dardel

1. **Tools without modules.** KMC, miRTrace and mergeTags are unlikely to exist as modules. Install them once into a prefix:
   ```bash
   bash container/install_tools.sh /cfs/klemming/projects/snic/<project>/programs/surfr2-tools --with-pigz
   ```
   This needs `curl`, `unzip`, `gcc`, `make` and zlib headers; if gcc is missing, load `PrgEnv-gnu`.
2. **Modules for the rest.** Find versions with `module spider samtools`, `module spider R` and `module spider java`, and pin them with explicit versions.
3. **R packages.** Install into a personal library, one per R minor version, because packages are built against it:
   ```bash
   module load PDC R/<version>
   Rscript container/install_r_packages.R /cfs/klemming/projects/snic/<project>/programs/R-lib-<R version>
   ```
4. **Config.**
   ```yaml
   execution:
     environment: modules
     modules: [PDC/<ver>, samtools/<ver>, R/<ver>, java/<ver>]
     path_prepend: [/cfs/klemming/projects/snic/<project>/programs/surfr2-tools/bin]
     r_libs: /cfs/klemming/projects/snic/<project>/programs/R-lib-<R version>
   ```
5. **Check.** Run `bash run_SURFR2.sh -c my.yaml --dry-run`. The tool check prints what resolved, and names anything missing with the fix.

The validator itself needs PyYAML on the login node: `pip install --user pyyaml`, or pass `--container`.

To add a tool during development, install it, add it to `env.sh`'s inputs (a module, `path_prepend` or `tools.*`), and add a check to `bin/surfr2_check_tools.sh`. Before a production run, add it to `install_tools.sh`, `install_r_packages.R` or the Dockerfile and rebuild the image.

### Memory and miRTrace

The sample step runs miRTrace (Java) and then KMC, so its job needs `max(qc.mirtrace_memory_gb, kmer.kmc_memory_gb) + 1G`. Set `execution.slurm.sample.mem` explicitly; the validator checks it, and warns when it is missing.

SURFR2 does not use miRTrace's own Python wrapper, for two reasons found on Dardel:
- **It sizes Java's heap from the whole node.** It sets the heap to half of the node's physical RAM, not the job's allocation, and miRTrace grows into all of it (it sizes its hash table from the start-up heap). On a shared node, a library with many distinct sequences then exceeds the job's memory and is OOM-killed.
- **It discards Java's exit status,** so a killed run looks successful.

`container/install_tools.sh` installs a small launcher instead. It runs `java -Xms<N>g -Xmx<N>g -jar mirtrace.jar` with N = `qc.mirtrace_memory_gb` and passes the exit status through. The tool check refuses any other `mirtrace`.

Each sample also verifies that miRTrace finished. The QC-passed FASTA must contain exactly the reads miRTrace's own statistics (written at the end of a run) report as passing. A truncated output therefore stops the sample instead of passing silently. If a SLURM job reports an `oom_kill` event, it is marked `OUT_OF_MEMORY` and the later steps are cancelled; raise the memory settings and relaunch.

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

- both environments: the modules environment with real module loading (simulated Lmod) and an unknown module; the container environment with every step routed through the image (simulated Singularity); a missing tool or image stops the run before any job is submitted
- generated job scripts run as SLURM array tasks, with index offsets
- `install_tools.sh`: downloads, checksums (a tampered checksum aborts), builds and smoke tests; the official KMC binaries pass the full regression suite
- the Dockerfile lints clean (hadolint). It has not been built here (no Docker in this sandbox): build it once and run `bash tests/run_test.sh` against the image (`SURFR2_TEST_ENV=container SURFR2_TEST_CONTAINER=<sif>`).

Not yet tested on real data or a real cluster. For a first real run, include a known positive control (e.g. a tissue-specific miRNA) and check that it is recovered.

## Fixes relative to SURFR1

These were verified against the KMC 3.2.4 and dekupl-mergeTags sources.

SURFR2 also corrects several issues found in SURFR1's scripts.

1. **mergeTags `-n` means unstranded merging** (`case 'n': stranded = 0`). SURFR1 used it on stranded libraries. SURFR2 passes `-n` only when `kmer.canonical: true`.
2. **mergeTags reads DE-kupl's column layout.** Value column 1 is read as a p-value, and the lowest becomes the contig representative. Value column 4 is read as log2FC, and it is read unconditionally, so a table with fewer than 4 value columns is read out of bounds. SURFR1's first value column was the raw TCGA cancer count, so each contig was represented by its least abundant k-mer. SURFR2 writes `rank, mean case CPM, mean control CPM, log2FC`.
3. **Enrichment ignored group size and depth.** `sample_ratio` was computed but never used.
4. **Wrong flag comments.** KMC `-b` is "no canonical form", not RAM-only (the flag itself was right). `dump -s` means sorted. `-cs4294967296` exceeds the 32-bit counter maximum of 2³²−1.
5. **Module loading inside the container.** `ml PDC` ran inside `singularity exec`, where Lmod is normally unavailable. SURFR2 loads modules on the host before `singularity exec`.
6. **Silent miRTrace failures.** miRTrace's wrapper discards Java's exit status, so an aborted or killed run exits 0, and it sizes the heap from the whole node (see "Memory and miRTrace"). SURFR2 replaces the wrapper and verifies each run's output against miRTrace's own statistics. Set `qc.phred_offset` if quality-encoding auto-detection fails.


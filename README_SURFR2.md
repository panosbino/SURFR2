# SURFR2

Reference-free discovery of case-specific small-RNA sequences from **per-sample** k-mer counts, with **read-depth normalisation** and a **config-driven** comparison. SURFR2 succeeds SURFR (Kalogeropoulos et al., 2025), which counted k-mers on BAMs pooled per condition.

## Quick start

```bash
cp config/config.example.yaml my_config.yaml        # edit: project, outdir, comparison, slurm
cp config/samples.example.tsv my_samples.tsv        # one row per sample
bash run_SURFR2.sh -c my_config.yaml --dry-run      # validate + show jobs, submit nothing
bash run_SURFR2.sh -c my_config.yaml                # submit (run on the login node)
```

If the login node's `python3` lacks PyYAML, pass `--container <sandbox>` after `ml PDC singularity`. The launcher then runs the validator inside the container.

## Inputs

**Samplesheet** (TSV, exact header):

| sample_id | cohort | condition | file_type | path |
|---|---|---|---|---|
| TCGA-05-4244-01A | TCGA | cancer | bam | /path/TCGA-05-4244-01A.bam |
| C3L-00001-N | CPTAC | adjacent_normal | fastq | /path/C3L-00001-N.fastq.gz |

- `file_type` is `bam` or `fastq`. Relative paths resolve against the samplesheet's directory.
- `sample_id` must be unique.
- A condition not named in the config is ignored, with a warning.

**Config** (`config/config.example.yaml`, every key documented there):

- **Comparison:** `comparison.case` and `comparison.control` name conditions from the samplesheet.
- **External controls:** `external_controls` lists conditions pooled across cohorts in which a k-mer must be (near-)absent, e.g. non-cancer SRA libraries. This replaces SURFR1's pooled SRA table with per-sample, normalised data.
- **Replication:** `replicate_in: all` requires a k-mer to pass in every cohort that has at least `min_samples_per_group` case and control samples. This generalises SURFR1's TCGA ∩ CPTAC intersection.

The validator rejects unknown keys, underpowered cohorts, and pre-filter settings stricter than the final filter, before anything is submitted.

## Steps

| Step | Script | Unit | What it does |
|---|---|---|---|
| sample | `bin/surfr2_sample.sh` | SLURM array task / sample | BAM→FASTQ (drops secondary/supplementary, restores read orientation) → miRTrace QC → KMC (`-ci1`, stranded) |
| matrix | `bin/surfr2_matrix.sh` | 1 job | Candidate k-mers (≥ `min_count` reads in ≥ `min_case_samples` case samples), computed with KMC set operations; per-sample candidate counts; reference k-mers for `median_ratio` |
| filter | `bin/surfr2_filter.R` | 1 job | Normalisation, per-cohort filters, replication, external controls, dekupl-mergeTags, plots |

**Resumable.** Every sample and the matrix step write a parameter fingerprint to `.done`. Re-launching reuses work only if the inputs and parameters are unchanged. Changing an input file, `k`, or the QC settings recomputes the affected samples. Changing a pre-filter parameter recomputes only the matrix.

## Method

**Normalisation.** Each sample gets an effective library size *L*, and CPM = count / *L* × 10⁶.

- **`cpm`**: *L* = total k-mers counted in that sample.
- **`median_ratio`**: DESeq2 median-of-ratios size factors on k-mers present in every sample, rescaled to *L* = *s* · geomean(total k-mers). CPM thresholds therefore mean the same under both methods.

Prefer `median_ratio` when a few very abundant miRNAs differ between conditions. Such differences deflate every other CPM in one group, which is a composition bias that plain depth normalisation cannot remove. `plots/normalisation.pdf` shows how far the size factors depart from depth.

**Filters**, applied per replicate cohort on CPM:

- **Case prevalence:** the fraction of case samples with CPM ≥ `min_cpm_case` is at least `min_prevalence_case`.
- **Case abundance:** mean case CPM ≥ `min_mean_cpm_case`.
- **Control prevalence:** the fraction of control samples with CPM ≥ `max_cpm_control` is at most `max_prevalence_control`. A small allowance is reasonable because adjacent-normal tissue can contain tumour cells.
- **Fold change:** log2((mean_case + pc) / (mean_control + pc)) ≥ `min_log2fc`.
- **External controls:** prevalence ≤ `max_prevalence_external`.

**Thresholds are placeholders.** SURFR1's cut-offs (>200 counts, enrichment >40, <100 in adjacent) were set on raw pooled sums and do not transfer to per-sample CPM. They must be recalibrated, for example against a held-out cohort. The filters also provide no false-discovery control. See the limitations below.

## Outputs (`<outdir>/results/`)

| File | Content |
|---|---|
| `case_specific_sequences.tsv` | mergeTags contigs, the main result |
| `case_specific_kmers.tsv` | per-cohort statistics of the passing k-mers |
| `case_specific_kmers_{counts,cpm}.tsv` | k-mer × sample matrices (all samples) |
| `candidate_kmer_stats.tsv.gz` | statistics and per-filter pass flags for every candidate, for auditing why a k-mer failed |
| `normalisation.tsv` | library sizes, size factors, effective library sizes |
| `run_summary.txt` | counts and any warnings (always read this) |
| `plots/` | library sizes, normalisation, violin of passing k-mers, cohort Venn |

Per-sample miRTrace reports are in `<outdir>/samples/<id>/mirtrace/`. Each run's resolved config and logs are in `<outdir>/runs/<timestamp>/`.

## Container change

The config validator needs PyYAML. Add `python3-yaml \` to the apt list in step 1 of the Dockerfile, and add this line to the smoke tests:

```
python3 -c "import yaml" && \
```

## Fixes relative to SURFR1

These were verified against the KMC 3.2.4 and dekupl-mergeTags sources.

1. **mergeTags `-n` means unstranded merging** (`case 'n': stranded = 0`). SURFR1 used it on stranded libraries. SURFR2 passes `-n` only when `kmer.canonical: true`.
2. **mergeTags reads DE-kupl's column layout.** Value column 1 is read as a p-value, and the lowest becomes the contig representative. Value column 4 is read as log2FC, and it is read unconditionally, so a table with fewer than 4 value columns is read out of bounds. SURFR1's first value column was the raw TCGA cancer count, so each contig was represented by its least abundant k-mer. SURFR2 writes `rank, mean case CPM, mean control CPM, log2FC`.
3. **Enrichment ignored group size and depth.** `sample_ratio` was computed but never used.
4. **Wrong flag comments.** KMC `-b` is "no canonical form", not RAM-only (the flag itself was right). `dump -s` means sorted. `-cs4294967296` exceeds the 32-bit counter maximum of 2³²−1.
5. **Module loading inside the container.** `ml PDC` ran inside `singularity exec`, where Lmod is normally unavailable. SURFR2 loads modules on the host before `singularity exec`.
6. **Silent miRTrace failures.** miRTrace 1.0.1 can exit 0 after aborting, e.g. when PHRED auto-detection fails. SURFR2 checks both the output and the log. Set `qc.phred_offset` if auto-detection fails.

## Validation status and known limitations

**Tested** on a synthetic 33-sample dataset (2 cohorts plus external controls) with 5 planted sequences, of which exactly 2 should pass. Both were recovered, and the 3 decoys were rejected for the intended reasons (present in controls, present in one cohort only, present in external controls). Also tested:

- exact read counts, including exclusion of secondary alignments
- restoration of reverse-strand reads to forward orientation
- `cpm` and `median_ratio`, the empty-result path, and the pre-filter warning
- resume and invalidation
- SLURM submission logic against mocked `sbatch`/`scontrol`, including splitting arrays above MaxArraySize

**Not yet tested:** real GDC data and a real SLURM cluster.

Before production:

- **Run a small subset first.** Try e.g. 10 cases and 10 controls per cohort with `--from sample` and check the miRTrace reports.
- **Memory on Dardel.** `shared` gives memory per core, so `kmer.kmc_memory_gb` must fit in `sample.cpus` × memory-per-core.
- **Filter-step memory.** The R step holds roughly n_candidates × n_groups × 24 bytes. Check the candidate count in `matrix.log` before raising `filter` resources.
- **Library-size denominator.** Total k-mers means a read of length *L* contributes *L − k + 1* k-mers. If read-length distributions differ systematically between conditions (e.g. RNA degradation), `cpm` denominators shift. `median_ratio` is more robust to this.
- **No statistical test yet.** Filters give ranking and specificity, not error control. With per-sample counts, a per-k-mer test (e.g. Wilcoxon or a negative-binomial GLM on candidates, with BH correction) is now possible and is the natural next step.
- **Unpaired design.** Tumour/normal pairs from the same patient are treated as independent samples. A paired design would need a `subject_id` column.

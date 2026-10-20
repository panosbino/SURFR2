# SURFR2

**Reference-free comparison of small-RNA sequencing data between two conditions.**

SURFR2 counts every 17-nucleotide sequence (k-mer) in every sample, removes sequencing artifacts, normalises for library size, and reports the sequences that differ between your two conditions. It needs no genome and no miRNA annotation, so it also finds unannotated small RNAs, fragments of other RNAs, and sequences from other organisms (e.g. a pathogen).

Two analyses are available:

- **Differential expression (`dea`, default):** sequences significantly higher in *either* condition, tested with edgeR.
- **Condition-specific (`specific`):** sequences present in one condition and absent from the other.

This page takes you from your FASTQ files to results on **Dardel**. The full reference for every option is in [`README_SURFR2.md`](README_SURFR2.md).

---

## Is my data suitable?

| Requirement | Why |
|---|---|
| **Small-RNA sequencing** (miRNA-seq, EV small RNA, …) | SURFR2 is built for short inserts (~18–40 nt) |
| **Two conditions** | e.g. treated vs untreated, tumour vs normal |
| **≥ 2 samples per condition**, ideally ≥ 3 | Needed for the statistical test. With 1 vs 1 SURFR2 still runs, but only reports fold changes, and the results are marked exploratory. |
| **One species** | k-mers differ between species, so mixed species make every comparison noisy |
| **FASTQ (gzipped or not) or BAM**, single-end | Paired-end small-RNA data: use read 1 |

If your samples are **paired** (e.g. treated and untreated from the same donor or batch), SURFR2 can use that; see step 2.

---

## One-time setup

**1. Get the code** into your project space:

```bash
P=/cfs/klemming/projects/snic/<your project>/programs   # pick any folder you own
mkdir -p $P && cd $P
unzip SURFR2.zip                                         # creates $P/SURFR2
```

**2. Make sure Python can read the config** (on the login node):

```bash
python3 -c "import yaml" || pip install --user pyyaml
```

**3. Check that you can read the container image**, which holds all the tools (samtools, miRTrace, KMC, R, edgeR, …):

```bash
ls -l /cfs/klemming/projects/snic/naiss2024-6-235/programs/surfr2_0.1.sif
```

If you get "Permission denied", ask whoever sent you SURFR2 for access or for a copy of the image.

---

## Running SURFR2

### Step 1: Check how your library was made

Raw small-RNA reads usually still contain the sequencing adapter, and some kits add random bases or UMIs (molecular barcodes). SURFR2 must know which. Run this on **one** of your files:

```bash
bash $P/SURFR2/bin/surfr2_detect_adapter.sh /path/to/one_sample.fastq.gz
```

It reads the first 200,000 reads and ends with a line starting with `==>` that tells you what to put in the config:

| Verdict | What to do |
|---|---|
| "reads are already adapter-trimmed" | `adapter: null` |
| TruSeq / NEBNext adapter | `adapter: <the sequence it prints>` |
| QIAseq miRNA | `adapter: AACTGTAGGCACCATCAAT` and `umi_length: 12` |
| NEXTflex (random bases) | **Stop:** not supported yet. Please tell us. |
| no known adapter, long reads | Find the adapter in the kit's documentation or the paper's methods |

### Step 2: Make a samplesheet

A tab-separated file with one row per sample. Column order doesn't matter:

```
sample_id	condition	file_type	path
ctrl_1	untreated	fastq	/cfs/klemming/projects/snic/<project>/data/ctrl_1.fastq.gz
ctrl_2	untreated	fastq	/cfs/klemming/projects/snic/<project>/data/ctrl_2.fastq.gz
ctrl_3	untreated	fastq	/cfs/klemming/projects/snic/<project>/data/ctrl_3.fastq.gz
trt_1	treated	fastq	/cfs/klemming/projects/snic/<project>/data/trt_1.fastq.gz
trt_2	treated	fastq	/cfs/klemming/projects/snic/<project>/data/trt_2.fastq.gz
trt_3	treated	fastq	/cfs/klemming/projects/snic/<project>/data/trt_3.fastq.gz
```

- **`sample_id`:** unique; letters, digits, `.`, `_`, `-` only.
- **`condition`:** your group names. These names appear in all results.
- **`file_type`:** `fastq` or `bam`.
- **`path`:** absolute, or relative to the samplesheet.
- **Paired design (optional):** add a `block` column, e.g. the donor or batch (`donor1`, `donor1`, `donor2`, …). Each block should contain both conditions.

### Step 3: Write the config

Copy the Dardel template next to your samplesheet and edit the lines marked `<-- EDIT`:

```bash
cp $P/SURFR2/config/config.dardel.yaml my_config.yaml
```

The settings you must set:

```yaml
project: my_project                    # a short name
samplesheet: samples.tsv               # the file from step 2
outdir: /cfs/klemming/projects/snic/<project>/SURFR2_results/my_project
analysis: dea                          # or: specific

comparison:
  case: treated                        # exactly as in the 'condition' column
  control: untreated                   # fold changes are log2(case / control)

qc:
  adapter: null                        # from step 1
  umi_length: 0                        # 12 for QIAseq

execution:
  slurm:
    account: <your allocation>         # e.g. naiss2026-x-yyy
```

The template already contains the container settings for Dardel. If something is left unfilled (`<…>`), SURFR2 stops and tells you which line.

### Step 4: Check everything (dry run)

From the login node:

```bash
bash $P/SURFR2/run_SURFR2.sh -c my_config.yaml --dry-run
```

Nothing is submitted. The dry run:
- checks the config and samplesheet, and that every input file exists;
- prints your design, e.g. `treated=3  untreated=3  [test]`. `[test]` means the edgeR test will run; `[descriptive]` means it can't, usually because a group has a single sample;
- starts the container and checks every tool inside it;
- lists the jobs it would submit.

Fix anything it reports before continuing.

### Step 5: Run

```bash
bash $P/SURFR2/run_SURFR2.sh -c my_config.yaml
squeue -u $USER                        # watch the jobs
```

Run this on the login node, **not** with `sbatch`; the script submits the jobs itself. You get one array job processing every sample in parallel, followed by three short jobs. Six samples typically finish within an hour, depending on library size and queue.

Logs are in `<outdir>/runs/latest/logs/`, one per sample and one per step.

**If something fails:** fix the cause and run the same command again. Finished samples are reused, so only the failed work is redone.

### Step 6: Look at the results

Everything is in `<outdir>/results/`:

| File | What it is |
|---|---|
| **`run_summary.txt`** | **Read this first.** Design, settings, how many sequences were found, warnings |
| **`de_sequences.tsv`** | **The main result:** overlapping k-mers merged into sequences, with fold change, p-value, the condition they are higher in (`higher_in`), and the samples they are found in. With `analysis: specific`: `specific_sequences.tsv` |
| `de_kmers.tsv` | the individual significant k-mers, with all statistics |
| `de_kmers_counts.tsv`, `de_kmers_cpm.tsv` | their raw counts and normalised values (CPM) per sample |
| `kmer_stats.tsv.gz` | statistics for **every** tested k-mer, to look up why something was or wasn't called |
| `genome_mapping.tsv` | if genomes are configured: each sequence's best hit per genome (mismatches, hits, locus) |
| `normalisation.tsv` | per sample: reads used (after QC and artifact removal), artifact reads removed, size factor |
| `interactive/scatter_*.html`, `volcano_*.html` | interactive plots: open in any browser and hover over a point to see its counts in every sample |
| `plots/*.pdf` | the same plots for publication, plus QC plots |

**CPM** means reads containing the sequence per million reads. **`log2fc`** is log2(case / control): +1 means 2× higher in your case condition.

**Sanity check with a known sequence.** If you know a small RNA that should differ between your conditions (a published marker, a tissue-specific miRNA), check it explicitly:

```bash
SINGULARITY_TMPDIR=/tmp SURFR2_CONTROLS="miR-21=TAGCTTATCAGACTGATGTTGA" \
  singularity exec --no-mount bind-paths -B /cfs/klemming \
  /cfs/klemming/projects/snic/naiss2024-6-235/programs/surfr2_0.1.sif \
  Rscript $P/SURFR2/bin/surfr2_diagnose.R <outdir>/runs/latest
```

Use `name=SEQUENCE`, in DNA or RNA letters; several controls are comma-separated. Besides the control, it reports:
- whether your replicates agree with each other;
- whether samples group by condition;
- how many sequences are found in only one sample.

Plots go to `results/plots/diagnostics_main.pdf`.

**Explore interactively (optional).** Copy the `results/` folder to your laptop. With R and the packages `shiny`, `plotly`, `dplyr`, `readr`, `ggplot2` and `jsonlite` installed, run:

```bash
SURFR2_RESULTS=/path/to/results Rscript -e 'shiny::runApp("SURFR2/app", launch.browser = TRUE)'
```

Click any point to see its counts per sample, or search for a sequence.

---

## Optional: map the results to reference genomes

SURFR2 can map the sequences it finds to one or more genomes. The plots are then coloured by where each sequence maps: one colour per genome, one for "both genomes" (or "several" with three or more), and black for sequences that map to none. This shows at a glance, for example, which sequences come from the host and which from a pathogen.

Add to your config:

```yaml
genomes:
  - {name: human, bowtie_index: /path/to/bowtie1/hg38}      # an existing Bowtie 1 index, or
  - {name: toxoplasma, fasta: /path/to/toxoplasma.fasta}    # a FASTA: SURFR2 builds the index once
mapping:
  max_mismatches: 1          # 0 = exact matches only
```

- **Index:** the first run with a FASTA builds the index into `<outdir>/genomes/`. That takes minutes for small genomes and one to a few hours for human; later runs reuse it. Bowtie 2 indexes (`.bt2`) don't work, only Bowtie 1 (`.ebwt`).
- **Mapping:** sequences are matched end-to-end on both strands, with up to `max_mismatches` mismatches and no gaps. Only the best hits count: a sequence with an exact hit is reported as exact.
- **Results:** `de_sequences.tsv` gains, per genome, the mismatches of the best hit (0, 1, or empty if unmapped), the number of hits and the first locus, plus a `genome_category` column. `genome_mapping.tsv` lists the same in long format, and the run summary gives the counts.
- **Remapping only:** after adding a genome to an existing run, `run_SURFR2.sh -c my_config.yaml --from map` redoes just the mapping and the plots.

---

## Changing settings and re-running

**After changing anything, run the same command again.** SURFR2 records what each step was computed with. Samples and candidate k-mers are recomputed only if their inputs or settings changed; everything else is reused.

Two shortcuts skip the sample jobs entirely:

| You changed only… | Faster command |
|---|---|
| `dea` or `statistics` thresholds | `run_SURFR2.sh -c my_config.yaml --from filter` |
| nothing; you just want the plots redrawn | `run_SURFR2.sh -c my_config.yaml --from report` |

Don't use `--from filter` after changing `analysis` or the `specific` thresholds: those can change which k-mers are tested, and a full run picks that up.

**Useful options to know** (all explained in [`config/config.example.yaml`](config/config.example.yaml)):

- **`dea.min_abs_log2fc`** (default 1 = 2-fold): the minimum fold change. Use 2–3 for 1-vs-1 designs.
- **`statistics.fdr`** (default 0.05): the false discovery rate.
- **`normalization.method: median_ratio`:** use when one condition is dominated by a few extremely abundant small RNAs.
- **`artifacts.extra_fasta`:** extra sequences to remove, e.g. spike-ins. Adapters, Illumina primers and PhiX are removed automatically.

---

## Troubleshooting

| Message | Cause and fix |
|---|---|
| `CONFIG ERROR: …` | The config or samplesheet has a problem; the message names it. Unedited `<…>` placeholders are caught here too |
| `host python3 lacks PyYAML` | `pip install --user pyyaml` (setup step 2) |
| `tool check failed` | The container couldn't start or a tool is missing; see the log it names. If the message mentions mounting or extraction, make sure you're on Dardel and haven't changed `container_options` |
| Jobs were `CANCELLED`, a sample log says `oom_kill` | A sample ran out of memory: increase `mem` under `slurm.sample` and run again |
| `miRTrace reported an error` | Often a wrong adapter, or a quality-encoding problem: set `qc.phred_offset: 33` |
| Very few reads pass QC | Wrong or missing `qc.adapter`; repeat step 1 |
| `run_summary.txt` says EXPLORATORY | At least one condition has a single sample, so no statistical test was possible |
| Nothing significant | Run the sanity check in step 6: if replicates disagree strongly, the design (not SURFR2) is usually the cause |

For singularity commands you run yourself on Dardel, always add `SINGULARITY_TMPDIR=/tmp` and `--no-mount bind-paths`, as in the examples above.

---

## Feedback

SURFR2 is under active development, and your experience is the most useful input we can get. Please send:

1. Your `my_config.yaml` and samplesheet.
2. `results/run_summary.txt` and `results/tool_versions.txt`.
3. If something failed: the error message and the log file it points to.
4. Your impressions:
   - Where did you get stuck or confused?
   - Do the results make biological sense? Did known markers come out?
   - What's missing from the results that you'd need?
   - Was anything in the plots or tables unclear?

"This step was annoying" is useful feedback too.

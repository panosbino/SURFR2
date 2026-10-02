#!/usr/bin/env python3
"""Synthetic small-RNA datasets with planted ground truth for SURFR2 regression testing.

Usage: make_synthetic.py <scenario> ; writes data/, samples.tsv and truth.tsv into the
current directory. truth.tsv columns: name, sequence, expected call in the 'specific'
analysis with a test (pass|fail), the same in descriptive mode, expected call in the 'dea'
analysis (up = higher in the case condition, down, none), and why.

Counts are random: each sequence's abundance varies between samples (log-normal
biological noise, CV ~30%) and reads are Poisson-sampled. Without biological noise,
dispersion estimates collapse to ~0 and any test calls trivial differences significant.

Scenarios
  multi   two cohorts (8+6, 7+6) + 6 external controls; BAM and FASTQ inputs
  small   3 treated vs 3 untreated, paired by donor (block), single cohort
  single  1 treated vs 1 untreated (descriptive mode)
  null    3 vs 3, no true differences (calibration: expect ~no FDR discoveries)
  umi     the 'small' design as raw QIAseq reads: <insert><adapter><12-nt UMI>..., 101 nt,
          uneven PCR duplication per molecule, plus adapter dimers, reads without adapter
          and reads with a truncated UMI. data/<sample>.nmol = true number of molecules.
"""
import gzip
import os
import subprocess
import sys

import numpy as np

SCENARIO = sys.argv[1] if len(sys.argv) > 1 else "multi"
rng = np.random.default_rng(42)
os.makedirs("data", exist_ok=True)


def rnd(n):
    return "".join(rng.choice(list("ACGT"), n))


def rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


N_BASE = 300 if SCENARIO == "null" else 60
BASE = [(rnd(22), float(rng.lognormal(0, 1.5))) for _ in range(N_BASE)]  # shared "miRNAs"
BIO_SD = 0.3                                                             # log-scale noise


def sample_counts(depth, rel):
    """rel: {seq: expected reads per 1e5}; returns {seq: count} with biological noise + Poisson."""
    out = {}
    for s, r in rel.items():
        c = int(rng.poisson(depth * r / 1e5 * rng.lognormal(0, BIO_SD)))
        if c:
            out[s] = c
    return out


def rel_from(weights):
    tot = sum(weights.values())
    return {s: 1e5 * w / tot for s, w in weights.items()}


def write_sample(sid, counts, as_bam=False):
    if as_bam:  # aligned-style BAM: half the reads on '-' (SEQ reverse-complemented) + a secondary
        sam = f"data/{sid}.sam"
        with open(sam, "w") as f:
            f.write("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000000\n")
            j = 0
            for s, c in counts.items():
                for _ in range(c):
                    flag, seq = (16, rc(s)) if j % 2 else (0, s)
                    f.write(f"{sid}_{j}\t{flag}\tchr1\t100\t60\t{len(s)}M\t*\t0\t0\t{seq}\t{'I'*len(s)}\n")
                    j += 1
            s0 = next(iter(counts))
            f.write(f"{sid}_0\t256\tchr1\t500\t0\t{len(s0)}M\t*\t0\t0\t{s0}\t{'I'*len(s0)}\n")
        subprocess.run(["samtools", "view", "-b", "-o", f"data/{sid}.bam", sam], check=True)
        os.remove(sam)
        path, ftype = f"data/{sid}.bam", "bam"
    else:
        path, ftype = f"data/{sid}.fastq.gz", "fastq"
        with gzip.open(path, "wt", compresslevel=1) as f:
            j = 0
            for s, c in counts.items():
                rec = f"\n{s}\n+\n{'I'*len(s)}\n"
                for _ in range(c):
                    f.write(f"@{sid}_{j}{rec}")
                    j += 1
    with open(f"data/{sid}.nreads", "w") as f:
        f.write(str(sum(counts.values())))
    return ftype, path


QIA_ADAPTER = "AACTGTAGGCACCATCAAT"


def write_umi_sample(sid, counts):
    """QIAseq-style reads: each molecule gets a random UMI and an uneven number of PCR copies."""
    path = f"data/{sid}.fastq.gz"
    tail = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    n_mol = 0
    reads = []
    pairs = set()
    for s, c in counts.items():
        for _ in range(c):
            umi = rnd(12)
            pairs.add((s, umi))
            copies = max(1, int(rng.lognormal(1.2, 0.8)))      # uneven PCR duplication
            reads += [(s + QIA_ADAPTER + umi + tail + rnd(101))[:101]] * copies
            n_mol += 1
    n = len(reads)
    reads += [(QIA_ADAPTER + rnd(12) + tail + rnd(101))[:101] for _ in range(n // 30)]      # adapter dimers
    reads += [rnd(101) for _ in range(n // 50)]                                            # no adapter
    reads += [(rnd(75) + QIA_ADAPTER + rnd(12))[:101] for _ in range(n // 100)]            # UMI cut off
    order = rng.permutation(len(reads))
    with gzip.open(path, "wt", compresslevel=1) as f:
        for j, i in enumerate(order):
            f.write(f"@{sid}_{j}\n{reads[i]}\n+\n{'I' * 101}\n")
    with open(f"data/{sid}.nmol", "w") as f:       # molecules, and distinct (insert, UMI) pairs:
        f.write(f"{n_mol}\t{len(pairs)}")          # they differ by true UMI collisions
    with open(f"data/{sid}.nreads", "w") as f:
        f.write(str(len(reads)))
    return "fastq", path


rows, truth, planted = [], [], {}

if SCENARIO == "multi":
    header = "sample_id\tcohort\tcondition\tfile_type\tpath"
    ONCOMIR = rnd(22)   # strongly up in cases -> composition effect
    planted = {k: rnd(25) for k in ["S1_true", "S2_in_controls", "S3_cohortA_only",
                                    "S4_in_external", "S5_one_control"]}
    S = planted
    # S1 is truly specific but in only 5/7 cases of cohort B: exact test on detection gives
    # p ~ 0.016, edgeR p ~ 0.01 - too weak at 5% FDR. It passes the filters alone, and must
    # FAIL once the test is on. This documents the power limit for heterogeneous markers.
    truth = [("S1_true", "fail", "pass", "-", "specific, but 5/7 vs 0/6 in B is weak evidence (p~0.01)"),
             ("S2_in_controls", "fail", "fail", "-", "equally present in controls"),
             ("S3_cohortA_only", "fail", "fail", "-", "not replicated in cohort B"),
             ("S4_in_external", "fail", "fail", "-", "present in external controls"),
             ("S5_one_control", "pass", "pass", "-", "one control allowed by config")]
    for coh, cond, n in [("A", "cancer", 8), ("A", "adjacent_normal", 6), ("B", "cancer", 7),
                         ("B", "adjacent_normal", 6), ("SRA", "noncancer_sra", 6)]:
        for i in range(n):
            sid = f"{coh}_{cond}_{i}"
            case = cond == "cancer"
            rel = rel_from({**dict(BASE), ONCOMIR: 40.0 if case else 2.0})
            if case and i < int(0.75 * n):
                rel[S["S1_true"]] = 60
            if cond in ("cancer", "adjacent_normal"):
                rel[S["S2_in_controls"]] = 60
            if case and coh == "A":
                rel[S["S3_cohortA_only"]] = 60
            if case or cond == "noncancer_sra":
                rel[S["S4_in_external"]] = 60
            if case:
                rel[S["S5_one_control"]] = 60
            if cond == "adjacent_normal" and i == 0:
                rel[S["S5_one_control"]] = 20
            counts = sample_counts(int(rng.integers(20000, 100000)), rel)
            ftype, path = write_sample(sid, counts, as_bam=(coh == "B" and case))
            rows.append(f"{sid}\t{coh}\t{cond}\t{ftype}\t{path}")

elif SCENARIO in ("small", "null", "umi"):
    header = "sample_id\tcondition\tblock\tfile_type\tpath"
    planted = {k: rnd(25) for k in ["T1_true", "T2_signif_not_specific",
                                    "T3_specific_too_low", "T4_two_of_three", "T5_down"]}
    T = planted
    if SCENARIO in ("small", "umi"):
        truth = [("T1_true", "pass", "pass", "up", "all treated, no untreated"),
                 ("T2_signif_not_specific", "fail", "fail", "up", "8x up but present in all untreated"),
                 ("T3_specific_too_low", "fail", "fail", "none", "treated-only but too few reads to be a candidate"),
                 ("T4_two_of_three", "fail", "fail", "none", "in 2 of 3 treated: below the 3-sample candidate filter"),
                 ("T5_down", "fail", "fail", "down", "8x higher in untreated")]
    for d in range(3):
        donor = rng.lognormal(0, 0.2, N_BASE)        # donor effect shared by the pair
        for cond in ("treated", "untreated"):
            sid = f"{cond}_d{d}"
            rel = rel_from({s: w * donor[i] for i, (s, w) in enumerate(BASE)})
            if SCENARIO in ("small", "umi"):
                tr = cond == "treated"
                if tr:
                    rel[T["T1_true"]] = 60
                    rel[T["T3_specific_too_low"]] = 0.6
                    if d < 2:
                        rel[T["T4_two_of_three"]] = 60
                rel[T["T2_signif_not_specific"]] = 480 if tr else 60
                rel[T["T5_down"]] = 60 if tr else 480
            counts = sample_counts(int(rng.integers(40000, 100000)), rel)
            if SCENARIO == "umi":
                ftype, path = write_umi_sample(sid, counts)
            else:
                ftype, path = write_sample(sid, counts)
            rows.append(f"{sid}\t{cond}\tdonor{d}\t{ftype}\t{path}")

elif SCENARIO == "single":
    header = "sample_id\tcondition\tfile_type\tpath"
    planted = {"T1_true": rnd(25), "T2_in_both": rnd(25)}
    truth = [("T1_true", "pass", "pass", "up", "treated only"),
             ("T2_in_both", "fail", "fail", "none", "equally present in untreated")]
    for cond in ("treated", "untreated"):
        rel = rel_from(dict(BASE))
        if cond == "treated":
            rel[planted["T1_true"]] = 60
        rel[planted["T2_in_both"]] = 60
        ftype, path = write_sample(f"{cond}_1", sample_counts(60000, rel))
        rows.append(f"{cond}_1\t{cond}\t{ftype}\t{path}")
else:
    sys.exit(f"unknown scenario: {SCENARIO}")

with open("samples.tsv", "w") as f:
    f.write(header + "\n" + "\n".join(rows) + "\n")
with open("truth.tsv", "w") as f:
    for name, exp_test, exp_desc, exp_dea, why in truth:
        f.write(f"{name}\t{planted[name]}\t{exp_test}\t{exp_desc}\t{exp_dea}\t{why}\n")
print(f"{SCENARIO}: {len(rows)} samples, {len(truth)} planted sequences")

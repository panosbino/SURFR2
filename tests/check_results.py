#!/usr/bin/env python3
"""Compare a SURFR2 run on synthetic data with its planted truth. Exit 1 on any failure.

Usage: check_results.py <scenario> <outdir>
truth.tsv columns: name, sequence, specific (with test), specific (descriptive), dea, why
"""
import csv, gzip, os, re, sys

scenario, res = sys.argv[1], os.path.join(sys.argv[2], "results")
fails = []


def check(ok, msg):
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        fails.append(msg)


def tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


summary = open(os.path.join(res, "run_summary.txt")).read()
analysis = "dea" if "analysis: differential expression" in summary else "specific"
descriptive = "(test)" not in summary          # every compared cohort ran without a test
prefix = "de" if analysis == "dea" else "specific"
cmp_ = re.search(r"comparison: (\S+) vs (\S+)", summary)
CASE, CTRL = cmp_.group(1), cmp_.group(2)
print(f"     analysis: {analysis}; mode: {'descriptive' if descriptive else 'test'}; {CASE} vs {CTRL}")

truth = [l.rstrip("\n").split("\t") for l in open("truth.tsv") if l.strip()]
stats = tsv(os.path.join(res, "kmer_stats.tsv.gz"))
seq_file = os.path.join(res, f"{prefix}_sequences.tsv")
seqs = tsv(seq_file) if os.path.exists(seq_file) else []
called = {r["contig"]: r for r in seqs}

# ---- planted truth -----------------------------------------------------------
if analysis == "specific":
    for name, seq, exp_test, exp_desc, _, why in truth:
        expect = exp_desc if descriptive else exp_test
        found = seq in called
        check(found == (expect == "pass"), f"{name}: {'recovered' if found else 'rejected'} (expected {expect}: {why})")
    n_expect = sum(t[3 if descriptive else 2] == "pass" for t in truth)
    check(len(seqs) == n_expect, f"{len(seqs)} sequence(s) reported, expected {n_expect}")
else:
    for name, seq, _, _, exp_dea, why in truth:
        got = "none" if seq not in called else ("up" if called[seq]["higher_in"] == CASE else "down")
        check(got == exp_dea, f"{name}: called {got} (expected {exp_dea}: {why})")
    planted = {t[1] for t in truth}
    extra = [c for c in called if c not in planted]
    print(f"     {len(extra)} called sequence(s) outside the planted set")
    if scenario in ("small_dea", "umi", "artifacts"):
        check(len(extra) <= 1, "at most one chance call among background sequences at FDR 5%")
    if scenario == "null_dea":
        check(len(seqs) == 0 and not any(r["pass_all"] == "TRUE" for r in stats),
              "no differentially expressed k-mers without true differences")

# ---- results are labelled with condition names, never case/control -------------
cols = list(stats[0].keys()) + (list(seqs[0].keys()) if seqs else [])
bad = [c for c in cols if re.search(r"(^|_)(case|control)(_|$)", c)]
check(not bad, f"no 'case'/'control' in column names {bad if bad else ''}".rstrip())
check(all(f"found_in_{c}" in stats[0] for c in (CASE, CTRL)), f"found_in_{CASE} and found_in_{CTRL} columns present")

# ---- found_in lists are exact for a planted sequence ------------------------------
samples = list(csv.DictReader(open(os.path.join(sys.argv[2], "runs", "latest", "samples.tsv")), delimiter="\t"))
by_cond = {}
for s in samples:
    by_cond.setdefault(s["condition"], []).append(s["sample_id"])
t1 = next((t for t in truth if t[0] == "T1_true"), None)
if t1:
    row = next(r for r in stats if r["kmer"] in t1[1])
    check(row[f"found_in_{CASE}"] == ",".join(by_cond[CASE]) and row[f"found_in_{CTRL}"] == "",
          f"T1 found_in lists exact ({CASE}: all; {CTRL}: none)")

# ---- library size is QC-passed reads (cpm normalisation) -------------------------
if "normalisation: cpm" in summary:
    norm = tsv(os.path.join(res, "normalisation.tsv"))
    check(all(float(r["eff_libsize"]) == float(r["qc_reads"]) for r in norm), "library size = QC-passed reads")

# ---- plots -------------------------------------------------------------------------
html = os.listdir(os.path.join(res, "interactive")) if os.path.isdir(os.path.join(res, "interactive")) else []
check(any(h.startswith("scatter_") for h in html), "interactive scatter written")
if descriptive:
    check(not any(h.startswith("volcano_") for h in html), "no volcano plot without a test")
    check("EXPLORATORY" in summary, "run summary flags results as EXPLORATORY")
    check(not any(k.endswith("_fdr") for k in stats[0]), "no p-values reported in descriptive mode")
else:
    check(any(h.startswith("volcano_") for h in html), "interactive volcano written")
    check("EXPLORATORY" not in summary, "run is not flagged exploratory")

# ---- scenario-specific ---------------------------------------------------------------
if scenario == "artifacts":
    for smp in samples:
        sid = smp["sample_id"]
        st = dict(l.rstrip("\n").split("\t") for l in open(os.path.join(sys.argv[2], "samples", sid, "artifact_stats.tsv")))
        injected = int(open(f"data/{sid}.nart").read())
        check(int(st["artifact_reads"]) == injected,
              f"{sid}: {st['artifact_reads']} artifact reads removed of {st['qc_passed_reads']} (injected: {injected})")
    rc = lambda x: x[::-1].translate(str.maketrans("ACGT", "TGCA"))
    srcs = ["".join(l.strip() for l in open("contaminant.fa") if not l.startswith(">"))]
    if os.environ.get("SURFR2_PHIX"):
        op = gzip.open if os.environ["SURFR2_PHIX"].endswith(".gz") else open
        srcs.append("".join(l.strip() for l in op(os.environ["SURFR2_PHIX"], "rt") if not l.startswith(">")).upper())
    art = {s_[i:i + 17] for src in srcs for s_ in (src, rc(src)) for i in range(len(s_) - 16)}
    left = sum(r["kmer"] in art for r in stats)
    check(left == 0, f"no artifact k-mer among the {len(stats)} candidates ({left} found)")
    print(f"     artifact references checked: {len(srcs)} ({'with' if len(srcs) > 1 else 'without'} PhiX)")
if scenario == "umi":
    for smp in samples:
        sid = smp["sample_id"]
        st = dict(l.rstrip("\n").split("\t") for l in open(os.path.join(sys.argv[2], "samples", sid, "umi_stats.tsv")))
        n_mol, n_pairs = map(int, open(f"data/{sid}.nmol").read().split())
        check(int(st["unique_molecules"]) == n_pairs,
              f"{sid}: {st['unique_molecules']} unique molecules from {st['input_reads']} reads "
              f"(true distinct insert+UMI pairs: {n_pairs}; molecules {n_mol}, {n_mol - n_pairs} UMI collisions)")
if scenario == "small":
    t2 = next(t[1] for t in truth if t[0] == "T2_signif_not_specific")
    rows = [r for r in stats if r["kmer"] in t2]
    check(rows and all(float(r["main_fdr"]) <= 0.05 for r in rows), "T2 k-mers are statistically significant ...")
    check(rows and all(r["pass_all"] == "FALSE" for r in rows), "... yet excluded by the specificity filters")
if scenario in ("null", "null_dea"):
    tested = [r for r in stats if r.get("main_pvalue") not in (None, "", "NA")]
    n_sig = sum(float(r["main_fdr"]) <= 0.05 for r in tested)
    frac = sum(float(r["main_pvalue"]) < 0.05 for r in tested) / len(tested)
    print(f"     null: {len(tested)} tested k-mers, {n_sig} at FDR <= 0.05, {frac:.3f} with p < 0.05")
    check(n_sig == 0, "no FDR discoveries without true differences")
    check(frac <= 0.15, "fraction p < 0.05 is not inflated (<= 0.15; k-mers are correlated in groups)")
sys.exit(1 if fails else 0)

#!/usr/bin/env python3
"""Compare a SURFR2 run on synthetic data with its planted truth. Exit 1 on any failure."""
import csv, gzip, os, sys

scenario, res = sys.argv[1], os.path.join(sys.argv[2], "results")
fails = []
def check(ok, msg):
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok: fails.append(msg)
def tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

summary = open(os.path.join(res, "run_summary.txt")).read()
descriptive = "(test)" not in summary          # every compared cohort ran without a test
truth = [l.rstrip("\n").split("\t") for l in open("truth.tsv") if l.strip()]
print(f"     mode: {'descriptive' if descriptive else 'test'}")
seq_file = os.path.join(res, "case_specific_sequences.tsv")
contigs = {r["contig"] for r in tsv(seq_file)} if os.path.exists(seq_file) else set()
for name, seq, exp_test, exp_desc, why in truth:
    expect = exp_desc if descriptive else exp_test
    found = seq in contigs
    check(found == (expect == "pass"), f"{name}: {'recovered' if found else 'rejected'} (expected {expect}: {why})")
n_expect = sum(t[3 if descriptive else 2] == "pass" for t in truth)
check(len(contigs) == n_expect, f"{len(contigs)} sequence(s) reported, expected {n_expect}")

stats = tsv(os.path.join(res, "candidate_kmer_stats.tsv.gz"))
if scenario in ("multi", "multi_descriptive", "small"):
    for f in os.listdir(os.path.join(res, "interactive")):
        check(os.path.getsize(os.path.join(res, "interactive", f)) > 0, f"interactive/{f} written")
    if scenario != "multi_descriptive":
        check("EXPLORATORY" not in summary, "run is not flagged exploratory")
if scenario == "small":
    t2 = dict((t[0], t[1]) for t in truth)["T2_signif_not_specific"]
    rows = [r for r in stats if r["kmer"] in t2]
    check(rows and all(float(r["main_fdr"]) <= 0.05 for r in rows),
          "T2 k-mers are statistically significant (FDR <= 0.05) ...")
    check(rows and all(r["pass_all"] == "FALSE" for r in rows), "... yet excluded by the specificity filters")
if scenario in ("single", "multi_descriptive"):
    check("EXPLORATORY" in summary, "run summary flags results as EXPLORATORY")
    check(not any(k.endswith("_fdr") for k in stats[0]), "no p-values reported in descriptive mode")
if scenario == "null":
    tested = [r for r in stats if r["main_pvalue"] not in ("", "NA")]
    n_sig = sum(float(r["main_fdr"]) <= 0.05 for r in tested)
    frac = sum(float(r["main_pvalue"]) < 0.05 for r in tested) / len(tested)
    print(f"     null: {len(tested)} tested k-mers, {n_sig} at FDR <= 0.05, {frac:.3f} with p < 0.05")
    check(n_sig == 0, "no FDR discoveries without true differences")
    check(frac <= 0.15, "fraction p < 0.05 is not inflated (<= 0.15; k-mers are correlated in groups)")
sys.exit(1 if fails else 0)

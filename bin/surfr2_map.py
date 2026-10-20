#!/usr/bin/env python3
"""
SURFR2 - map the called (merged) sequences to one or more reference genomes.

Usage: surfr2_map.py <run_config_dir>

For every genome in the config (resolved_config.json -> genomes):
  * a Bowtie 1 index is used as given (bowtie_index), or built once from the FASTA into
    <outdir>/genomes/<name>/ and reused while its fingerprint matches;
  * the merged sequences (results/<de|specific>_sequences.tsv) are mapped end-to-end
    with at most mapping.max_mismatches mismatches (bowtie -v, no gaps), both strands,
    reporting only the best mismatch level (--best --strata), up to max_hits + 1 hits.

Adds to the sequences table, per genome <g>:
    <g>_mismatches   0 / 1 / ... for the best hit, empty if unmapped
    <g>_hits         number of best hits (">max_hits" when capped)
    <g>_locus        first best hit as chrom:position(1-based):strand
and overall:
    genome_category  <g> | both (two genomes configured, both hit) |
                     multiple (more than two configured, two or more hit) | unmapped
Writes results/genome_mapping.tsv (one row per sequence and genome) and appends a
mapping block to results/run_summary.txt.
"""

import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] [map] {msg}", file=sys.stderr, flush=True)


def die(msg):
    log(f"ERROR: {msg}")
    sys.exit(1)


def run(cmd, **kw):
    log("  $ " + " ".join(cmd))
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, **kw)
    if r.returncode != 0:
        die(f"{cmd[0]} failed (exit {r.returncode}):\n{r.stderr.strip()[-2000:]}")
    return r


def bowtie_version():
    out = run(["bowtie", "--version"]).stdout.splitlines()
    return out[0].split()[-1] if out else "unknown"


def ensure_index(genome, outdir, threads, bt_version):
    """Return the index prefix for a genome, building and caching it from FASTA if needed."""
    if "bowtie_index" in genome:
        return genome["bowtie_index"]
    fasta = genome["fasta"]
    st = os.stat(fasta)
    fingerprint = f"fasta={os.path.realpath(fasta)} size={st.st_size} mtime={int(st.st_mtime)} bowtie={bt_version}"
    gdir = os.path.join(outdir, "genomes", genome["name"])
    prefix = os.path.join(gdir, "index")
    done = os.path.join(gdir, ".done")
    if os.path.exists(done) and open(done).read().strip() == fingerprint:
        log(f"{genome['name']}: index up to date")
        return prefix
    log(f"{genome['name']}: building Bowtie index from {fasta} (can take hours for large genomes)")
    shutil.rmtree(gdir, ignore_errors=True)
    os.makedirs(gdir)
    run(["bowtie-build", "--threads", str(threads), "-q", fasta, prefix])
    with open(done + ".part", "w") as fh:
        fh.write(fingerprint + "\n")
    os.replace(done + ".part", done)
    return prefix


def map_sequences(prefix, fasta_in, max_mm, max_hits, threads):
    """Bowtie -v: end-to-end, no gaps, <= max_mm mismatches, best stratum only."""
    r = run(["bowtie", "-f", "-v", str(max_mm), "-k", str(max_hits + 1), "--best", "--strata",
             "-p", str(threads), "-x", prefix, fasta_in])          # Bowtie 1 prints aligned reads only
    hits = {}
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 4:
            continue
        name, strand, chrom, offset = f[0], f[1], f[2], int(f[3])
        mm = len(f[7].split(",")) if len(f) > 7 and f[7] else 0
        h = hits.setdefault(name, {"n": 0, "mm": mm, "locus": f"{chrom.split()[0]}:{offset + 1}:{strand}"})
        h["n"] += 1
        h["mm"] = min(h["mm"], mm)
    return hits


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: surfr2_map.py <run_config_dir>")
    run_dir = sys.argv[1]
    cfg = json.load(open(os.path.join(run_dir, "resolved_config.json")))
    genomes = cfg.get("genomes") or []
    if not genomes:
        log("no genomes configured - nothing to do")
        return
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK") or os.cpu_count() or 1)
    max_mm = int(cfg["mapping"]["max_mismatches"])
    max_hits = int(cfg["mapping"]["max_hits"])
    prefix = "specific" if cfg.get("analysis") == "specific" else "de"
    res_dir = os.path.join(cfg["outdir"], "results")
    seq_file = os.path.join(res_dir, f"{prefix}_sequences.tsv")
    if not os.path.exists(seq_file):
        log(f"{seq_file} not found (no sequences were called) - nothing to map")
        return
    with open(seq_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols, rows = list(reader.fieldnames), list(reader)
    if not rows:
        log("no sequences to map")
        return
    # drop columns of an earlier mapping run, so re-running is idempotent
    names = [g["name"] for g in genomes]
    stale = [c for c in cols if c == "genome_category" or any(c.startswith(f"{n}_") and
             c[len(n) + 1:] in ("mismatches", "hits", "locus") for n in names)]
    cols = [c for c in cols if c not in stale]

    bt_version = bowtie_version()
    log(f"Bowtie {bt_version}; {len(rows)} sequences; {len(genomes)} genome(s); <= {max_mm} mismatch(es)")
    tmp = tempfile.mkdtemp(prefix="surfr2_map_", dir=res_dir)
    try:
        fa = os.path.join(tmp, "sequences.fa")
        with open(fa, "w") as fh:
            for i, r in enumerate(rows):
                fh.write(f">s{i}\n{r['contig']}\n")
        long_rows = []
        for g in genomes:
            idx = ensure_index(g, cfg["outdir"], threads, bt_version)
            hits = map_sequences(idx, fa, max_mm, max_hits, threads)
            n_mapped = 0
            for i, r in enumerate(rows):
                h = hits.get(f"s{i}")
                if h:
                    n_mapped += 1
                    n = f">{max_hits}" if h["n"] > max_hits else str(h["n"])
                    r[f"{g['name']}_mismatches"], r[f"{g['name']}_hits"], r[f"{g['name']}_locus"] = str(h["mm"]), n, h["locus"]
                else:
                    r[f"{g['name']}_mismatches"] = r[f"{g['name']}_hits"] = r[f"{g['name']}_locus"] = ""
                long_rows.append([r["contig"], g["name"], r[f"{g['name']}_mismatches"],
                                  r[f"{g['name']}_hits"], r[f"{g['name']}_locus"]])
            log(f"{g['name']}: {n_mapped} of {len(rows)} sequences map with <= {max_mm} mismatch(es)")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    counts = {}
    for r in rows:
        hit = [n for n in names if r[f"{n}_mismatches"] != ""]
        if not hit:
            cat = "unmapped"
        elif len(hit) == 1:
            cat = hit[0]
        else:
            cat = "both" if len(names) == 2 else "multiple"
        r["genome_category"] = cat
        counts[cat] = counts.get(cat, 0) + 1

    out_cols = cols + [f"{n}_{x}" for n in names for x in ("mismatches", "hits", "locus")] + ["genome_category"]
    tmp_out = seq_file + ".part"
    with open(tmp_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp_out, seq_file)
    with open(os.path.join(res_dir, "genome_mapping.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["contig", "genome", "mismatches", "hits", "locus"])
        w.writerows(long_rows)

    order = names + ["both", "multiple", "unmapped"]
    summary = ["", f"genome mapping (Bowtie {bt_version}, <= {max_mm} mismatch(es), best hits only):"]
    summary += [f"  {c}: {counts[c]} sequence(s)" for c in order if c in counts]
    for g in genomes:
        exact = sum(r[f"{g['name']}_mismatches"] == "0" for r in rows)
        mm = sum(r[f"{g['name']}_mismatches"] not in ("", "0") for r in rows)
        summary.append(f"  {g['name']}: {exact} exact, {mm} with mismatches "
                       f"({'index ' + g['bowtie_index'] if 'bowtie_index' in g else 'from ' + g['fasta']})")
    # replace (not append) the mapping block, so '--from map' re-runs stay clean
    sfile = os.path.join(res_dir, "run_summary.txt")
    old = open(sfile).read().split("\n\ngenome mapping (")[0].rstrip("\n") if os.path.exists(sfile) else ""
    with open(sfile, "w") as fh:
        fh.write(old + "\n" + "\n".join(summary) + "\n")
    log("done: " + ", ".join(f"{c}={counts[c]}" for c in order if c in counts))


if __name__ == "__main__":
    main()

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
and, for a genome with a GTF (genomes: [{..., gtf: annotation.gtf[.gz]}]):
    <g>_gene_orientation  sense | antisense | intergenic (mapped, no gene overlapped);
                          empty if unmapped. Sense wins: antisense genes are reported
                          only when no best hit overlaps a gene on its own strand.
    <g>_gene_name         overlapped genes (gene_name, else gene_id), ';'-separated,
    <g>_gene_id           ordered by how many best hits overlap them, then by gene
    <g>_gene_biotype      length (smallest = most specific first, e.g. a miRNA before
                          its host gene); at most 5, then '+N more'. Biotype is
                          gene_biotype (Ensembl) or gene_type (GENCODE).
    <g>_annotated_hits    best hits overlapping a reported gene / best hits examined
  Every best hit is annotated (up to max_hits + 1), not only the reported locus, so
  multi-copy small RNAs (miRNA paralogues, tRNAs, snoRNAs) are annotated by all copies.
  Genes are spans from the GTF 'gene' lines, or from the extent of a gene_id's other
  features when the GTF has no 'gene' lines. Chromosome names are matched with and
  without a 'chr' prefix (and M <-> MT).
and overall:
    genome_category  <g> | both (two genomes configured, both hit) |
                     multiple (more than two configured, two or more hit) | unmapped
Writes results/genome_mapping.tsv (one row per sequence and genome) and appends a
mapping block to results/run_summary.txt.
"""

import bisect
import csv
import gzip
import json
import os
import re
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
    """Bowtie -v: end-to-end, no gaps, <= max_mm mismatches, best stratum only.
    Returns {name: {n, mm, locus, all: [(chrom, start0, strand), ...]}}."""
    r = run(["bowtie", "-f", "-v", str(max_mm), "-k", str(max_hits + 1), "--best", "--strata",
             "-p", str(threads), "-x", prefix, fasta_in])          # Bowtie 1 prints aligned reads only
    hits = {}
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 4:
            continue
        name, strand, chrom, offset = f[0], f[1], f[2], int(f[3])
        mm = len(f[7].split(",")) if len(f) > 7 and f[7] else 0
        chrom = chrom.split()[0]
        h = hits.setdefault(name, {"n": 0, "mm": mm, "locus": f"{chrom}:{offset + 1}:{strand}", "all": []})
        h["n"] += 1
        h["mm"] = min(h["mm"], mm)
        h["all"].append((chrom, offset, strand))
    return hits


# -----------------------------------------------------------------------------
# GTF annotation
# -----------------------------------------------------------------------------
GTF_ATTR = re.compile(r'(\S+) "([^"]*)"')
BIN = 100_000
ANNOT_COLS = ("gene_orientation", "gene_name", "gene_id", "gene_biotype", "annotated_hits")
MAX_GENES = 5


class GeneIndex:
    """Gene spans from a GTF, binned per chromosome for overlap queries."""

    def __init__(self, path):
        opener = gzip.open if path.endswith(".gz") else open
        from_gene_lines, from_features = {}, {}
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                attrs = dict(GTF_ATTR.findall(f[8]))
                gid = attrs.get("gene_id")
                if not gid:
                    continue
                chrom, start, end, strand = f[0], int(f[3]), int(f[4]), f[6]
                name = attrs.get("gene_name") or attrs.get("gene") or gid
                biotype = attrs.get("gene_biotype") or attrs.get("gene_type") or ""
                if f[2] == "gene":
                    from_gene_lines[(chrom, gid)] = [chrom, start, end, strand, gid, name, biotype]
                else:
                    g = from_features.get((chrom, gid))
                    if g is None:
                        from_features[(chrom, gid)] = [chrom, start, end, strand, gid, name, biotype]
                    else:
                        g[1], g[2] = min(g[1], start), max(g[2], end)
                        if not g[6]:
                            g[6] = biotype
        for key, g in from_features.items():          # genes without a 'gene' line
            from_gene_lines.setdefault(key, g)
        self.genes = list(from_gene_lines.values())   # [chrom, start1, end1, strand, id, name, biotype]
        self.bins = {}
        for i, (chrom, s, e, *_r) in enumerate(self.genes):
            b = self.bins.setdefault(chrom, {})
            for k in range((s - 1) // BIN, (e - 1) // BIN + 1):
                b.setdefault(k, []).append(i)
        self.chroms = set(self.bins)

    def resolve_chrom(self, chrom):
        """GTF chromosome name for a genome (FASTA) chromosome name, or None."""
        if chrom in self.chroms:
            return chrom
        alts = [chrom[3:]] if chrom.startswith("chr") else ["chr" + chrom]
        alts += {"M": ["MT", "chrMT"], "MT": ["M", "chrM"], "chrM": ["MT"], "chrMT": ["M"]}.get(chrom, [])
        return next((a for a in alts if a in self.chroms), None)

    def overlaps(self, chrom, start0, length):
        """Indices of genes overlapping the 0-based half-open interval [start0, start0 + length)."""
        c = self.resolve_chrom(chrom)
        if c is None:
            return []
        s1, e1 = start0 + 1, start0 + length                 # 1-based inclusive, like the GTF
        b = self.bins[c]
        found = set()
        for k in range((s1 - 1) // BIN, (e1 - 1) // BIN + 1):
            for i in b.get(k, ()):
                g = self.genes[i]
                if g[1] <= e1 and s1 <= g[2]:
                    found.add(i)
        return found


def annotate(index, seq, hit):
    """Annotation columns (dict keyed by ANNOT_COLS) for one mapped sequence."""
    sense, anti = {}, {}                     # gene index -> number of best hits overlapping it
    n_sense_hits = n_anti_hits = 0
    for chrom, start0, strand in hit["all"]:
        s_here = a_here = False
        for i in index.overlaps(chrom, start0, len(seq)):
            g_strand = index.genes[i][3]
            if g_strand == strand or g_strand not in "+-":     # unstranded genes count as sense
                sense[i] = sense.get(i, 0) + 1; s_here = True
            else:
                anti[i] = anti.get(i, 0) + 1; a_here = True
        n_sense_hits += s_here
        n_anti_hits += a_here
    if sense:
        chosen, orient, n_annot = sense, "sense", n_sense_hits
    elif anti:
        chosen, orient, n_annot = anti, "antisense", n_anti_hits
    else:
        return {"gene_orientation": "intergenic", "gene_name": "", "gene_id": "", "gene_biotype": "",
                "annotated_hits": f"0/{len(hit['all'])}"}
    order = sorted(chosen, key=lambda i: (-chosen[i], index.genes[i][2] - index.genes[i][1], index.genes[i][5]))
    def joined(col, uniq=False):
        vals = [index.genes[i][col] for i in order]
        if uniq:
            vals = list(dict.fromkeys(v for v in vals if v))
        more = len(vals) - MAX_GENES
        return ";".join(vals[:MAX_GENES]) + (f";+{more} more" if more > 0 else "")
    return {"gene_orientation": orient, "gene_name": joined(5), "gene_id": joined(4),
            "gene_biotype": joined(6, uniq=True), "annotated_hits": f"{n_annot}/{len(hit['all'])}"}


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
             c[len(n) + 1:] in ("mismatches", "hits", "locus") + ANNOT_COLS for n in names)]
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
            gn = g["name"]
            idx = ensure_index(g, cfg["outdir"], threads, bt_version)
            hits = map_sequences(idx, fa, max_mm, max_hits, threads)
            index = None
            if g.get("gtf"):
                log(f"{gn}: reading genes from {g['gtf']}")
                index = GeneIndex(g["gtf"])
                log(f"{gn}: {len(index.genes):,} genes on {len(index.chroms)} chromosome(s)")
                hit_chroms = {c for h in hits.values() for c, _s, _st in h["all"]}
                unmatched = sorted(c for c in hit_chroms if index.resolve_chrom(c) is None)
                if hit_chroms and len(unmatched) == len(hit_chroms):
                    log(f"WARNING: {gn}: none of the hit chromosomes ({', '.join(sorted(hit_chroms)[:5])} ...) "
                        f"occur in the GTF (e.g. {', '.join(sorted(index.chroms)[:5])}): wrong GTF for this genome?")
                elif unmatched:
                    log(f"{gn}: hit chromosomes not in the GTF (left unannotated): {', '.join(unmatched[:10])}"
                        + (" ..." if len(unmatched) > 10 else ""))
            n_mapped = 0
            for i, r in enumerate(rows):
                h = hits.get(f"s{i}")
                if h:
                    n_mapped += 1
                    n = f">{max_hits}" if h["n"] > max_hits else str(h["n"])
                    r[f"{gn}_mismatches"], r[f"{gn}_hits"], r[f"{gn}_locus"] = str(h["mm"]), n, h["locus"]
                else:
                    r[f"{gn}_mismatches"] = r[f"{gn}_hits"] = r[f"{gn}_locus"] = ""
                if index is not None:
                    ann = annotate(index, r["contig"], h) if h else dict.fromkeys(ANNOT_COLS, "")
                    for c in ANNOT_COLS:
                        r[f"{gn}_{c}"] = ann[c]
                long_rows.append([r["contig"], gn, r[f"{gn}_mismatches"], r[f"{gn}_hits"], r[f"{gn}_locus"]]
                                 + [r.get(f"{gn}_{c}", "") for c in ANNOT_COLS])
            log(f"{gn}: {n_mapped} of {len(rows)} sequences map with <= {max_mm} mismatch(es)")
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

    annotated = [g["name"] for g in genomes if g.get("gtf")]
    out_cols = cols + [f"{n}_{x}" for n in names for x in ("mismatches", "hits", "locus")
                       + (ANNOT_COLS if n in annotated else ())] + ["genome_category"]
    tmp_out = seq_file + ".part"
    with open(tmp_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp_out, seq_file)
    with open(os.path.join(res_dir, "genome_mapping.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["contig", "genome", "mismatches", "hits", "locus"] + list(ANNOT_COLS))
        w.writerows(long_rows)

    order = names + ["both", "multiple", "unmapped"]
    summary = ["", f"genome mapping (Bowtie {bt_version}, <= {max_mm} mismatch(es), best hits only):"]
    summary += [f"  {c}: {counts[c]} sequence(s)" for c in order if c in counts]
    for g in genomes:
        exact = sum(r[f"{g['name']}_mismatches"] == "0" for r in rows)
        mm = sum(r[f"{g['name']}_mismatches"] not in ("", "0") for r in rows)
        summary.append(f"  {g['name']}: {exact} exact, {mm} with mismatches "
                       f"({'index ' + g['bowtie_index'] if 'bowtie_index' in g else 'from ' + g['fasta']})")
    for gn in annotated:
        mapped = [r for r in rows if r[f"{gn}_mismatches"] != ""]
        orient = {o: sum(r[f"{gn}_gene_orientation"] == o for r in mapped) for o in ("sense", "antisense", "intergenic")}
        summary.append(f"  {gn} annotation ({next(g['gtf'] for g in genomes if g['name'] == gn)}): "
                       + ", ".join(f"{v} {k}" for k, v in orient.items()))
        bt = {}
        for r in mapped:
            if r[f"{gn}_gene_orientation"] == "sense":
                b = r[f"{gn}_gene_biotype"].split(";")[0] or "(no biotype)"
                bt[b] = bt.get(b, 0) + 1
        if bt:
            top = sorted(bt.items(), key=lambda kv: -kv[1])
            summary.append(f"    sense, by biotype of the first gene: "
                           + ", ".join(f"{k} {v}" for k, v in top[:8]) + (", ..." if len(top) > 8 else ""))
    # replace (not append) the mapping block, so '--from map' re-runs stay clean
    sfile = os.path.join(res_dir, "run_summary.txt")
    old = open(sfile).read().split("\n\ngenome mapping (")[0].rstrip("\n") if os.path.exists(sfile) else ""
    with open(sfile, "w") as fh:
        fh.write(old + "\n" + "\n".join(summary) + "\n")
    log("done: " + ", ".join(f"{c}={counts[c]}" for c in order if c in counts))


if __name__ == "__main__":
    main()

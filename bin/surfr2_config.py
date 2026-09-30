#!/usr/bin/env python3
"""
SURFR2 - configuration validator and resolver.

Reads the user's YAML config and TSV samplesheet, validates them, and writes a
resolved run configuration that every downstream step consumes:

  <out>/resolved_config.json  - full config with defaults applied, absolute paths
  <out>/samples.tsv           - validated samplesheet, 1-based 'index' + 'role'
  <out>/params.env            - shell-quoted KEY=VALUE subset for bash steps

Validation is deliberately strict: an error here costs seconds, the same error
discovered after a 500-task array job costs a day.

Usage:
  surfr2_config.py CONFIG.yaml --out DIR [--no-file-check] [--container PATH]
"""

import argparse
import copy
import csv
import json
import math
import os
import re
import shlex
import sys

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("ERROR: PyYAML is required (apt: python3-yaml, pip: pyyaml).")

SAMPLE_COLUMNS = ["sample_id", "cohort", "condition", "file_type", "path"]
FILE_TYPES = {"bam": (".bam",), "fastq": (".fastq", ".fq", ".fastq.gz", ".fq.gz")}
ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

DEFAULTS = {
    "project": None,
    "samplesheet": None,
    "outdir": None,
    "comparison": {"case": None, "control": None,
                   "external_controls": [], "replicate_in": "all"},
    "kmer": {"k": 17, "canonical": False, "kmc_memory_gb": 8},
    "qc": {"mirtrace_species": "hsa", "adapter": None, "phred_offset": None},
    "candidates": {"min_count": 5, "min_case_samples": 10},
    "normalization": {"method": "cpm", "max_reference_kmers": 100000},
    "filters": {
        "min_samples_per_group": 5,
        "min_cpm_case": 1.0,
        "min_prevalence_case": 0.2,
        "min_mean_cpm_case": 2.0,
        "max_cpm_control": 0.5,
        "max_prevalence_control": 0.05,
        "min_log2fc": 3.0,
        "pseudocount_cpm": 0.1,
        "max_prevalence_external": 0.01,
    },
    "mergetags": {"enabled": True, "min_overlap": 8},
    "tools": {
        "samtools": "samtools",
        "pigz": "pigz",
        "mirtrace": "/opt/mirtrace/mirtrace",
        "kmc": "/opt/kmc/bin/kmc",
        "kmc_tools": "/opt/kmc/bin/kmc_tools",
        "mergetags": "/opt/dekupl/bin/mergeTags",
        "rscript": "Rscript",
    },
    "execution": {
        "executor": "slurm",
        "container": None,
        "bind": [],
        "keep_intermediates": False,
        "slurm": {
            "account": None,
            "modules": [],
            "sample": {"partition": "shared", "cpus": 8, "time": "04:00:00",
                       "array_throttle": 50},
            "matrix": {"partition": "memory", "cpus": 16, "time": "12:00:00"},
            "filter": {"partition": "memory", "cpus": 4, "time": "06:00:00"},
        },
    },
}


class ConfigError(Exception):
    pass


def merge_defaults(defaults, user, path=""):
    """Recursively overlay user values on defaults; reject unknown keys (typo guard)."""
    if user is None:
        return copy.deepcopy(defaults)
    if not isinstance(user, dict):
        raise ConfigError(f"'{path or 'config'}' must be a mapping")
    out = copy.deepcopy(defaults)
    for key, val in user.items():
        full = f"{path}.{key}" if path else key
        if key not in defaults:
            raise ConfigError(f"Unknown config key '{full}' (typo?)")
        if isinstance(defaults[key], dict) and defaults[key] and key != "tools":
            out[key] = merge_defaults(defaults[key], val, full)
        elif key == "tools":
            if not isinstance(val, dict):
                raise ConfigError("'tools' must be a mapping")
            unknown = set(val) - set(defaults["tools"])
            if unknown:
                raise ConfigError(f"Unknown tools: {sorted(unknown)}")
            out[key].update(val)
        else:
            out[key] = val
    return out


def req(cond, msg):
    if not cond:
        raise ConfigError(msg)


def check_number(cfg, section, key, lo=None, hi=None, integer=False):
    v = cfg[section][key]
    name = f"{section}.{key}"
    req(isinstance(v, (int, float)) and not isinstance(v, bool), f"{name} must be a number")
    if integer:
        req(float(v).is_integer(), f"{name} must be an integer")
    if lo is not None:
        req(v >= lo, f"{name} must be >= {lo} (got {v})")
    if hi is not None:
        req(v <= hi, f"{name} must be <= {hi} (got {v})")


def resolve(path, base):
    path = os.path.expanduser(str(path))
    return os.path.normpath(path if os.path.isabs(path) else os.path.join(base, path))


def validate_config(cfg):
    for key in ("project", "samplesheet", "outdir"):
        req(cfg[key], f"'{key}' is required")
    req(ID_RE.match(str(cfg["project"])), "project must match [A-Za-z0-9._-]+")

    cmp_ = cfg["comparison"]
    req(cmp_["case"] and cmp_["control"], "comparison.case and comparison.control are required")
    req(cmp_["case"] != cmp_["control"], "comparison.case and comparison.control must differ")
    req(isinstance(cmp_["external_controls"], list), "comparison.external_controls must be a list")
    overlap = {cmp_["case"], cmp_["control"]} & set(cmp_["external_controls"])
    req(not overlap, f"external_controls overlaps case/control: {sorted(overlap)}")
    ri = cmp_["replicate_in"]
    req(ri == "all" or (isinstance(ri, list) and ri), "comparison.replicate_in must be 'all' or a non-empty list")

    check_number(cfg, "kmer", "k", 9, 31, integer=True)
    check_number(cfg, "kmer", "kmc_memory_gb", 1, integer=True)
    req(isinstance(cfg["kmer"]["canonical"], bool), "kmer.canonical must be true/false")
    if cfg["kmer"]["canonical"]:
        print("WARNING: kmer.canonical=true merges a k-mer with its reverse complement. "
              "Small-RNA libraries are stranded; this is almost certainly not what you want.",
              file=sys.stderr)

    po = cfg["qc"]["phred_offset"]
    req(po in (None, 33, 64), "qc.phred_offset must be null (auto-detect), 33 or 64")
    ad = cfg["qc"]["adapter"]
    req(ad is None or re.match(r"^[ACGTNacgtn]+$", str(ad)), "qc.adapter must be null or a DNA sequence")

    check_number(cfg, "candidates", "min_count", 1, integer=True)
    check_number(cfg, "candidates", "min_case_samples", 1, integer=True)

    req(cfg["normalization"]["method"] in ("cpm", "median_ratio"),
        "normalization.method must be 'cpm' or 'median_ratio'")
    check_number(cfg, "normalization", "max_reference_kmers", 1000, integer=True)

    f = "filters"
    check_number(cfg, f, "min_samples_per_group", 2, integer=True)
    for key in ("min_cpm_case", "min_mean_cpm_case", "max_cpm_control"):
        check_number(cfg, f, key, 0)
    check_number(cfg, f, "min_log2fc", 0)
    req(cfg[f]["min_log2fc"] > 0, "filters.min_log2fc must be > 0 (case-specific = up in case; "
                                  "mergeTags also treats log2FC <= 0 as a separate 'down' set)")
    check_number(cfg, f, "pseudocount_cpm", 1e-9)
    for key in ("min_prevalence_case", "max_prevalence_control", "max_prevalence_external"):
        check_number(cfg, f, key, 0, 1)
    req(cfg[f]["max_cpm_control"] <= cfg[f]["min_cpm_case"],
        "filters.max_cpm_control > filters.min_cpm_case makes 'absent in control' weaker "
        "than 'present in case' - almost certainly a mistake")

    check_number(cfg, "mergetags", "min_overlap", 1, cfg["kmer"]["k"] - 1, integer=True)

    ex = cfg["execution"]
    req(ex["executor"] in ("slurm", "local"), "execution.executor must be 'slurm' or 'local'")
    req(isinstance(ex["bind"], list), "execution.bind must be a list")
    if ex["executor"] == "slurm":
        req(ex["slurm"]["account"], "execution.slurm.account is required for the slurm executor")
        for step in ("sample", "matrix", "filter"):
            s = ex["slurm"][step]
            req(isinstance(s.get("cpus"), int) and s["cpus"] >= 1, f"slurm.{step}.cpus must be int >= 1")
            req(re.match(r"^(\d+-)?\d{1,2}:\d{2}:\d{2}$", str(s.get("time", ""))),
                f"slurm.{step}.time must look like HH:MM:SS or D-HH:MM:SS")


def read_samplesheet(path, check_files):
    req(os.path.isfile(path), f"Samplesheet not found: {path}")
    base = os.path.dirname(path)
    rows, seen = [], set()
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        req(header is not None, "Samplesheet is empty")
        header = [h.strip() for h in header]
        req(header == SAMPLE_COLUMNS,
            f"Samplesheet header must be exactly: {chr(9).join(SAMPLE_COLUMNS)} (tab-separated); got {header}")
        for lineno, fields in enumerate(reader, start=2):
            if not fields or all(not x.strip() for x in fields) or fields[0].startswith("#"):
                continue
            req(len(fields) == len(SAMPLE_COLUMNS),
                f"Samplesheet line {lineno}: expected {len(SAMPLE_COLUMNS)} tab-separated fields, got {len(fields)}")
            r = dict(zip(SAMPLE_COLUMNS, (x.strip() for x in fields)))
            for col in ("sample_id", "cohort", "condition"):
                req(ID_RE.match(r[col]), f"Samplesheet line {lineno}: invalid {col} '{r[col]}' "
                                         "(allowed: letters, digits, . _ -)")
            req(r["sample_id"] not in seen, f"Samplesheet line {lineno}: duplicate sample_id '{r['sample_id']}'")
            seen.add(r["sample_id"])
            ft = r["file_type"].lower()
            req(ft in FILE_TYPES, f"Samplesheet line {lineno}: file_type must be one of {sorted(FILE_TYPES)}")
            r["file_type"] = ft
            r["path"] = resolve(r["path"], base)
            req(not re.search(r"[\s,'\"]", r["path"]),
                f"Samplesheet line {lineno}: path must not contain whitespace, commas or quotes")
            req(r["path"].lower().endswith(FILE_TYPES[ft]),
                f"Samplesheet line {lineno}: extension of '{r['path']}' does not match file_type '{ft}'")
            if check_files:
                req(os.path.isfile(r["path"]), f"Samplesheet line {lineno}: file not found: {r['path']}")
            rows.append(r)
    req(rows, "Samplesheet contains no samples")
    return rows


def assign_roles(cfg, rows):
    cmp_ = cfg["comparison"]
    role_of = {cmp_["case"]: "case", cmp_["control"]: "control"}
    role_of.update({c: "external" for c in cmp_["external_controls"]})
    ignored = sorted({r["condition"] for r in rows} - set(role_of))
    if ignored:
        print(f"WARNING: samples with conditions {ignored} are not part of this comparison "
              "and will be ignored.", file=sys.stderr)
    kept = [dict(r, role=role_of[r["condition"]]) for r in rows if r["condition"] in role_of]
    present = {r["condition"] for r in kept}
    for cond in role_of:
        req(cond in present, f"Condition '{cond}' from the config has no samples in the samplesheet")
    for i, r in enumerate(kept, start=1):
        r["index"] = i
    return kept


def check_design(cfg, rows):
    """Cohort-level checks, and that the raw-count pre-filter is looser than the final filter."""
    f = cfg["filters"]
    n = {}
    for r in rows:
        if r["role"] in ("case", "control"):
            n.setdefault(r["cohort"], {"case": 0, "control": 0})[r["role"]] += 1

    usable = sorted(c for c, v in n.items() if v["case"] >= f["min_samples_per_group"]
                    and v["control"] >= f["min_samples_per_group"])
    ri = cfg["comparison"]["replicate_in"]
    if ri == "all":
        cohorts = usable
        req(cohorts, f"No cohort has >= {f['min_samples_per_group']} case AND control samples")
        dropped = sorted(set(n) - set(usable))
        if dropped:
            print(f"WARNING: cohorts {dropped} lack enough case or control samples and are "
                  "excluded from replication (their samples still count toward candidates).",
                  file=sys.stderr)
    else:
        cohorts = sorted(ri)
        for c in cohorts:
            req(c in n, f"replicate_in cohort '{c}' not found among case/control samples")
            req(c in usable, f"replicate_in cohort '{c}' has {n[c]} samples; need >= "
                             f"{f['min_samples_per_group']} of both case and control")
    if len(cohorts) == 1:
        print(f"WARNING: only one cohort ({cohorts[0]}) - no independent replication. "
              "SURFR1's cross-cohort intersection is not reproduced.", file=sys.stderr)

    # A k-mer passing the final filter in every replicate cohort is expressed in at least
    # sum_c ceil(p * n_case_c) case samples. The pooled pre-filter must not demand more.
    p = f["min_prevalence_case"]
    implied = sum(max(1, math.ceil(p * n[c]["case"] - 1e-9)) for c in cohorts)
    mcs = cfg["candidates"]["min_case_samples"]
    req(mcs <= implied,
        f"candidates.min_case_samples={mcs} is stricter than the final filter implies "
        f"(min_prevalence_case={p} over replicate cohorts => >= {implied} case samples). "
        f"Set candidates.min_case_samples <= {implied}, or k-mers that would pass could be lost.")
    return n, cohorts


def write_outputs(cfg, rows, out):
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "resolved_config.json"), "w") as fh:
        json.dump(cfg, fh, indent=2, sort_keys=True)

    cols = ["index", "sample_id", "cohort", "condition", "role", "file_type", "path"]
    with open(os.path.join(out, "samples.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([r[c] for c in cols])

    ex, t = cfg["execution"], cfg["tools"]
    env = {
        "PROJECT": cfg["project"],
        "OUTDIR": cfg["outdir"],
        "RUN_CONFIG_DIR": out,
        "N_SAMPLES": len(rows),
        "N_CASE": sum(r["role"] == "case" for r in rows),
        "K": cfg["kmer"]["k"],
        "CANONICAL": str(cfg["kmer"]["canonical"]).lower(),
        "KMC_MEMORY_GB": cfg["kmer"]["kmc_memory_gb"],
        "MIRTRACE_SPECIES": cfg["qc"]["mirtrace_species"],
        "MIRTRACE_ADAPTER": cfg["qc"]["adapter"] or "",
        "MIRTRACE_PHRED": cfg["qc"]["phred_offset"] or "",
        "CAND_MIN_COUNT": cfg["candidates"]["min_count"],
        "CAND_MIN_CASE_SAMPLES": cfg["candidates"]["min_case_samples"],
        "NORM_METHOD": cfg["normalization"]["method"],
        "KEEP_INTERMEDIATES": str(ex["keep_intermediates"]).lower(),
        "EXECUTOR": ex["executor"],
        "CONTAINER": ex["container"] or "",
        "BIND": ",".join(ex["bind"]),
        "SLURM_ACCOUNT": ex["slurm"]["account"] or "",
        "SLURM_MODULES": " ".join(ex["slurm"]["modules"]),
        "SAMTOOLS": t["samtools"], "PIGZ": t["pigz"], "MIRTRACE": t["mirtrace"],
        "KMC": t["kmc"], "KMC_TOOLS": t["kmc_tools"], "RSCRIPT": t["rscript"],
    }
    for step in ("sample", "matrix", "filter"):
        s = ex["slurm"][step]
        env[f"SLURM_{step.upper()}_PARTITION"] = s.get("partition", "")
        env[f"SLURM_{step.upper()}_CPUS"] = s.get("cpus", 1)
        env[f"SLURM_{step.upper()}_TIME"] = s.get("time", "")
    env["SLURM_ARRAY_THROTTLE"] = ex["slurm"]["sample"].get("array_throttle", 50)

    with open(os.path.join(out, "params.env"), "w") as fh:
        fh.write("# Generated by surfr2_config.py - do not edit; re-run the launcher instead.\n")
        for k, v in env.items():
            fh.write(f"{k}={shlex.quote(str(v))}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("--out", required=True, help="directory for resolved run configuration")
    ap.add_argument("--no-file-check", action="store_true", help="skip checking that input files exist")
    ap.add_argument("--container", help="override execution.container")
    a = ap.parse_args()

    try:
        cfg_path = os.path.abspath(a.config)
        req(os.path.isfile(cfg_path), f"Config not found: {cfg_path}")
        with open(cfg_path) as fh:
            user = yaml.safe_load(fh)
        cfg = merge_defaults(DEFAULTS, user)
        validate_config(cfg)
        base = os.path.dirname(cfg_path)
        cfg["samplesheet"] = resolve(cfg["samplesheet"], base)
        cfg["outdir"] = resolve(cfg["outdir"], base)
        if a.container is not None:
            cfg["execution"]["container"] = a.container or None
        if cfg["execution"]["container"]:
            cfg["execution"]["container"] = resolve(cfg["execution"]["container"], base)
        cfg["config_file"] = cfg_path
        unsafe = re.compile(r"[\s'\"\\$`]")
        for label, val in [("outdir", cfg["outdir"])] + [(f"tools.{k}", v) for k, v in cfg["tools"].items()]:
            req(not unsafe.search(str(val)), f"{label} must not contain whitespace, quotes, '$', '`' or '\\': {val}")

        rows = assign_roles(cfg, read_samplesheet(cfg["samplesheet"], not a.no_file_check))
        counts, cohorts = check_design(cfg, rows)
        cfg["replicate_cohorts"] = cohorts
        write_outputs(cfg, rows, os.path.abspath(a.out))
    except ConfigError as e:
        sys.exit(f"CONFIG ERROR: {e}")

    ext = {}
    for r in rows:
        if r["role"] == "external":
            ext[r["condition"]] = ext.get(r["condition"], 0) + 1
    print(f"Project {cfg['project']}: {len(rows)} samples", file=sys.stderr)
    for c in sorted(counts):
        tag = "replicate" if c in cohorts else "candidates only"
        print(f"  {c:<12} case={counts[c]['case']:<5} control={counts[c]['control']:<5} [{tag}]", file=sys.stderr)
    for cond, k in sorted(ext.items()):
        print(f"  external '{cond}': {k}", file=sys.stderr)


if __name__ == "__main__":
    main()

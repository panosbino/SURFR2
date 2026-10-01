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

REQUIRED_COLUMNS = ["sample_id", "condition", "file_type", "path"]
OPTIONAL_COLUMNS = ["cohort", "block"]
DEFAULT_COHORT = "main"
MOVED_KEYS = {"execution.slurm.modules": "execution.modules"}
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
    "candidates": {"min_count": 3, "min_samples": "auto"},
    "normalization": {"method": "cpm", "max_reference_kmers": 100000},
    "filters": {
        "min_samples_per_group": 1,
        # sample-count thresholds: "all", an integer, or a fraction in (0, 1)
        "min_case_samples_detected": "all",
        "max_control_samples_detected": 0,
        "max_external_samples_detected": 0,
        "min_cpm_case": 1.0,
        "max_cpm_control": 0.5,
        "min_mean_cpm_case": 2.0,
        "min_log2fc": 2.0,
        "pseudocount_cpm": 0.1,
    },
    "statistics": {"test": "auto", "fdr": 0.05},
    "mergetags": {"enabled": True, "min_overlap": 8},
    # Plain names are resolved through PATH: inside the container (PATH set by the image)
    # or after module loads / path_prepend. Give an absolute path to pin a specific binary.
    "tools": {
        "samtools": "samtools",
        "pigz": "pigz",
        "mirtrace": "mirtrace",
        "kmc": "kmc",
        "kmc_tools": "kmc_tools",
        "mergetags": "mergeTags",
        "rscript": "Rscript",
    },
    "execution": {
        "executor": "slurm",
        # container: tools from a Singularity image; modules: tools from environment
        # modules and/or fixed binary paths (development - no image rebuilds)
        "environment": "container",
        "container": None,
        "bind": [],
        "modules": [],          # loaded before every step, in both environments
        "path_prepend": [],     # modules environment: directories put first on PATH
        "r_libs": None,         # modules environment: personal R library
        "keep_intermediates": False,
        "slurm": {
            "account": None,
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
            if full in MOVED_KEYS:
                raise ConfigError(f"'{full}' has moved to '{MOVED_KEYS[full]}'")
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


def check_count_spec(v, name, allow_zero):
    """'all', an integer (>= 0 or >= 1), or a fraction strictly between 0 and 1."""
    if v == "all":
        return
    ok_int = isinstance(v, int) and not isinstance(v, bool) and v >= (0 if allow_zero else 1)
    ok_frac = isinstance(v, float) and 0 < v < 1
    req(ok_int or ok_frac, f"{name} must be 'all', an integer{' >= 0' if allow_zero else ' >= 1'}, "
                           f"or a fraction between 0 and 1 (got {v!r})")


def resolve_min(spec, n):
    """Minimum number of samples required (case side)."""
    if spec == "all":
        return n
    if isinstance(spec, float):
        return max(1, math.ceil(spec * n - 1e-9))
    return spec


def resolve_max(spec, n):
    """Maximum number of samples allowed (control/external side)."""
    if spec == "all":
        return n
    if isinstance(spec, float):
        return math.floor(spec * n + 1e-9)
    return min(spec, n)


def matrix_rank(rows, tol=1e-9):
    """Rank by Gaussian elimination (design matrices are tiny; avoids a numpy dependency)."""
    m = [list(map(float, r)) for r in rows]
    rank, ncol = 0, len(m[0]) if m else 0
    for c in range(ncol):
        piv = max(range(rank, len(m)), key=lambda i: abs(m[i][c]), default=None)
        if piv is None or abs(m[piv][c]) < tol:
            continue
        m[rank], m[piv] = m[piv], m[rank]
        for i in range(len(m)):
            if i != rank and abs(m[i][c]) > tol:
                f = m[i][c] / m[rank][c]
                m[i] = [x - f * y for x, y in zip(m[i], m[rank])]
        rank += 1
    return rank


def design_rank(conditions, blocks):
    """Rank and column count of the model matrix ~ [block +] condition."""
    levels = sorted(set(blocks))[1:] if blocks is not None else []
    rows = [[1.0] + [1.0 if b == lv else 0.0 for lv in levels] + [1.0 if c == "case" else 0.0]
            for c, b in zip(conditions, blocks if blocks is not None else [None] * len(conditions))]
    return matrix_rank(rows), 2 + len(levels)


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
    ms = cfg["candidates"]["min_samples"]
    req(ms == "auto" or (isinstance(ms, int) and not isinstance(ms, bool) and ms >= 1),
        "candidates.min_samples must be 'auto' or an integer >= 1")

    req(cfg["normalization"]["method"] in ("cpm", "median_ratio"),
        "normalization.method must be 'cpm' or 'median_ratio'")
    check_number(cfg, "normalization", "max_reference_kmers", 1000, integer=True)

    f = "filters"
    check_number(cfg, f, "min_samples_per_group", 1, integer=True)
    for key in ("min_cpm_case", "min_mean_cpm_case", "max_cpm_control"):
        check_number(cfg, f, key, 0)
    check_number(cfg, f, "pseudocount_cpm", 1e-9)
    check_number(cfg, f, "min_log2fc", 0)
    req(cfg[f]["min_log2fc"] > 0, "filters.min_log2fc must be > 0 (case-specific = up in case; "
                                  "mergeTags also treats log2FC <= 0 as a separate 'down' set)")
    req(cfg[f]["max_cpm_control"] <= cfg[f]["min_cpm_case"],
        "filters.max_cpm_control > filters.min_cpm_case makes 'absent in control' weaker "
        "than 'present in case' - almost certainly a mistake")
    check_count_spec(cfg[f]["min_case_samples_detected"], "filters.min_case_samples_detected", allow_zero=False)
    check_count_spec(cfg[f]["max_control_samples_detected"], "filters.max_control_samples_detected", allow_zero=True)
    check_count_spec(cfg[f]["max_external_samples_detected"], "filters.max_external_samples_detected", allow_zero=True)

    st = cfg["statistics"]
    req(st["test"] in ("auto", "edger", "none"), "statistics.test must be 'auto', 'edger' or 'none'")
    check_number(cfg, "statistics", "fdr", 1e-12, 1)

    check_number(cfg, "mergetags", "min_overlap", 1, cfg["kmer"]["k"] - 1, integer=True)

    ex = cfg["execution"]
    req(ex["executor"] in ("slurm", "local"), "execution.executor must be 'slurm' or 'local'")
    req(ex["environment"] in ("container", "modules"), "execution.environment must be 'container' or 'modules'")
    for key in ("bind", "modules", "path_prepend"):
        req(isinstance(ex[key], list) and all(isinstance(x, str) and x for x in ex[key]),
            f"execution.{key} must be a list of non-empty strings")
    for m in ex["modules"]:
        req(re.match(r"^[A-Za-z0-9._/+-]+$", m), f"execution.modules: invalid module name '{m}'")
    unpinned = [m for m in ex["modules"] if "/" not in m]
    if unpinned and ex["environment"] == "modules":
        print(f"WARNING: modules without a version {unpinned} load the site default, which can change "
              "between runs; pin them (e.g. R/4.4.1). tool_versions.txt records what was loaded.",
              file=sys.stderr)
    if ex["environment"] == "container":
        req(ex["container"], "execution.environment 'container' needs execution.container "
                             "(or use environment: modules)")
        req(not ex["path_prepend"] and not ex["r_libs"],
            "execution.path_prepend and execution.r_libs apply to the 'modules' environment only; "
            "inside a container they would override the image's tools")
    elif ex["container"]:
        print("WARNING: execution.container is ignored in the 'modules' environment.", file=sys.stderr)
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
        missing = [c for c in REQUIRED_COLUMNS if c not in header]
        unknown = [c for c in header if c not in REQUIRED_COLUMNS + OPTIONAL_COLUMNS]
        req(not missing, f"Samplesheet is missing required column(s) {missing}; required: "
                         f"{REQUIRED_COLUMNS}, optional: {OPTIONAL_COLUMNS} (tab-separated)")
        req(not unknown, f"Samplesheet has unknown column(s) {unknown}; allowed: "
                         f"{REQUIRED_COLUMNS + OPTIONAL_COLUMNS}")
        req(len(set(header)) == len(header), "Samplesheet header has duplicate column names")
        for lineno, fields in enumerate(reader, start=2):
            if not fields or all(not x.strip() for x in fields) or fields[0].startswith("#"):
                continue
            req(len(fields) == len(header),
                f"Samplesheet line {lineno}: expected {len(header)} tab-separated fields, got {len(fields)}")
            r = dict(zip(header, (x.strip() for x in fields)))
            r["cohort"] = r.get("cohort") or DEFAULT_COHORT
            r["block"] = r.get("block", "")
            for col in ("sample_id", "cohort", "condition"):
                req(ID_RE.match(r[col]), f"Samplesheet line {lineno}: invalid {col} '{r[col]}' "
                                         "(allowed: letters, digits, . _ -)")
            if r["block"]:
                req(ID_RE.match(r["block"]), f"Samplesheet line {lineno}: invalid block '{r['block']}'")
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
    """Per-cohort design: group sizes, resolved thresholds, test/descriptive mode, blocking,
    and that the condition-blind pre-filter is looser than the final filter."""
    f, st = cfg["filters"], cfg["statistics"]
    n = {}
    for r in rows:
        if r["role"] in ("case", "control"):
            n.setdefault(r["cohort"], {"case": 0, "control": 0})[r["role"]] += 1

    mspg = f["min_samples_per_group"]
    usable = sorted(c for c, v in n.items() if v["case"] >= mspg and v["control"] >= mspg)
    ri = cfg["comparison"]["replicate_in"]
    if ri == "all":
        cohorts = usable
        req(cohorts, f"No cohort has >= {mspg} case AND >= {mspg} control samples")
        dropped = sorted(set(n) - set(usable))
        if dropped:
            print(f"WARNING: cohorts {dropped} lack case or control samples and are not compared "
                  "(their samples still count toward candidate selection).", file=sys.stderr)
    else:
        cohorts = sorted(ri)
        for c in cohorts:
            req(c in n, f"replicate_in cohort '{c}' not found among case/control samples")
            req(c in usable, f"replicate_in cohort '{c}' has {n[c]} samples; need >= {mspg} of each")

    design = {"cohorts": {}, "external": {}}
    for c in cohorts:
        nc, nn = n[c]["case"], n[c]["control"]
        need = resolve_min(f["min_case_samples_detected"], nc)
        allow = resolve_max(f["max_control_samples_detected"], nn)
        req(need <= nc, f"cohort '{c}': min_case_samples_detected={f['min_case_samples_detected']} "
                        f"exceeds its {nc} case sample(s)")
        grp = [r for r in rows if r["cohort"] == c and r["role"] in ("case", "control")]
        blocks = [r["block"] for r in grp]
        blocked = any(blocks)
        if blocked:
            req(all(blocks), f"cohort '{c}': 'block' must be set for all or none of its case/control samples")
            for b in sorted(set(blocks)):
                roles = {r["role"] for r in grp if r["block"] == b}
                if len(roles) < 2:
                    print(f"WARNING: cohort '{c}', block '{b}' contains only {roles.pop()} samples; "
                          "they cannot inform the case-control comparison in the test.", file=sys.stderr)
        rank, ncol = design_rank([r["role"] for r in grp], blocks if blocked else None)
        df = len(grp) - rank
        testable = nc >= 2 and nn >= 2 and rank == ncol and df >= 1
        if st["test"] == "none":
            mode = "descriptive"
        elif testable:
            mode = "test"
        elif st["test"] == "edger":
            raise ConfigError(
                f"cohort '{c}' cannot be tested (case={nc}, control={nn}, "
                f"{'condition confounded with block, ' if rank < ncol else ''}residual df={df}); "
                "statistics.test: edger requires >= 2 samples per group and >= 1 residual df")
        else:
            mode = "descriptive"
            why = ("condition is confounded with block" if rank < ncol else
                   "no residual degrees of freedom" if nc >= 2 and nn >= 2 else
                   "a group has a single sample")
            print(f"WARNING: cohort '{c}' runs in DESCRIPTIVE mode ({why}): no statistical test, "
                  "results are exploratory.", file=sys.stderr)
        design["cohorts"][c] = {"n_case": nc, "n_control": nn, "need_case": need,
                                "allow_control": allow, "mode": mode, "blocked": blocked,
                                "residual_df": df}

    for r in rows:
        if r["role"] == "external":
            design["external"].setdefault(r["condition"], {"n": 0})["n"] += 1
    for cond, v in design["external"].items():
        v["allow"] = resolve_max(f["max_external_samples_detected"], v["n"])

    # A k-mer passing the final filter is detected in >= sum_c need_case_c case samples;
    # the condition-blind pre-filter (over all case+control samples) must not demand more.
    implied = sum(d["need_case"] for d in design["cohorts"].values())
    ms = cfg["candidates"]["min_samples"]
    if ms == "auto":
        ms = implied
    req(ms <= implied,
        f"candidates.min_samples={ms} is stricter than the final filter implies (>= {implied} "
        "case samples detected); k-mers that would pass could be lost. Lower it or use 'auto'.")
    design["candidates_min_samples"] = ms
    return n, cohorts, design


def write_outputs(cfg, rows, out):
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "resolved_config.json"), "w") as fh:
        json.dump(cfg, fh, indent=2, sort_keys=True)

    # 'block' is last so positional readers of the first seven columns stay valid
    cols = ["index", "sample_id", "cohort", "condition", "role", "file_type", "path", "block"]
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
        "N_TESTED_GROUPS": sum(r["role"] in ("case", "control") for r in rows),
        "K": cfg["kmer"]["k"],
        "CANONICAL": str(cfg["kmer"]["canonical"]).lower(),
        "KMC_MEMORY_GB": cfg["kmer"]["kmc_memory_gb"],
        "MIRTRACE_SPECIES": cfg["qc"]["mirtrace_species"],
        "MIRTRACE_ADAPTER": cfg["qc"]["adapter"] or "",
        "MIRTRACE_PHRED": cfg["qc"]["phred_offset"] or "",
        "CAND_MIN_COUNT": cfg["candidates"]["min_count"],
        "CAND_MIN_SAMPLES": cfg["design"]["candidates_min_samples"],
        "NORM_METHOD": cfg["normalization"]["method"],
        "KEEP_INTERMEDIATES": str(ex["keep_intermediates"]).lower(),
        "EXECUTOR": ex["executor"],
        "ENVIRONMENT": ex["environment"],
        "CONTAINER": ex["container"] if ex["environment"] == "container" else "",
        "MODULES": " ".join(ex["modules"]),
        "PATH_PREPEND": ":".join(ex["path_prepend"]),
        "R_LIBS_DIR": ex["r_libs"] or "",
        "NEEDS_EDGER": str(any(d["mode"] == "test" for d in cfg["design"]["cohorts"].values())).lower(),
        "BIND": ",".join(ex["bind"]),
        "SLURM_ACCOUNT": ex["slurm"]["account"] or "",
        "SAMTOOLS": t["samtools"], "PIGZ": t["pigz"], "MIRTRACE": t["mirtrace"],
        "KMC": t["kmc"], "KMC_TOOLS": t["kmc_tools"], "RSCRIPT": t["rscript"],
        "MERGETAGS": t["mergetags"],
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
        if a.container is not None:            # before validation: it can supply a missing path
            cfg["execution"]["container"] = a.container or None
        validate_config(cfg)
        base = os.path.dirname(cfg_path)
        cfg["samplesheet"] = resolve(cfg["samplesheet"], base)
        cfg["outdir"] = resolve(cfg["outdir"], base)
        ex = cfg["execution"]
        if ex["container"]:
            ex["container"] = resolve(ex["container"], base)
        ex["path_prepend"] = [resolve(d, base) for d in ex["path_prepend"]]
        if ex["r_libs"]:
            ex["r_libs"] = resolve(ex["r_libs"], base)
        if not a.no_file_check:
            if ex["environment"] == "container":
                req(os.path.exists(ex["container"]), f"container not found: {ex['container']}")
            for d in ex["path_prepend"] + ([ex["r_libs"]] if ex["r_libs"] else []):
                req(os.path.isdir(d), f"directory not found: {d}")
        cfg["config_file"] = cfg_path
        unsafe = re.compile(r"[\s'\"\\$`]")
        for label, val in ([("outdir", cfg["outdir"])] + [(f"tools.{k}", v) for k, v in cfg["tools"].items()]
                           + [("execution.path_prepend", d) for d in cfg["execution"]["path_prepend"]]
                           + [("execution.r_libs", cfg["execution"]["r_libs"] or "")]
                           + [("execution.container", cfg["execution"]["container"] or "")]):
            req(not unsafe.search(str(val)), f"{label} must not contain whitespace, quotes, '$', '`' or '\\': {val}")

        rows = assign_roles(cfg, read_samplesheet(cfg["samplesheet"], not a.no_file_check))
        counts, cohorts, design = check_design(cfg, rows)
        cfg["replicate_cohorts"] = cohorts
        cfg["design"] = design
        write_outputs(cfg, rows, os.path.abspath(a.out))
    except ConfigError as e:
        sys.exit(f"CONFIG ERROR: {e}")

    ext = {}
    for r in rows:
        if r["role"] == "external":
            ext[r["condition"]] = ext.get(r["condition"], 0) + 1
    print(f"Project {cfg['project']}: {len(rows)} samples", file=sys.stderr)
    for c in sorted(counts):
        d = cfg["design"]["cohorts"].get(c)
        tag = (f"{d['mode']}{', blocked' if d['blocked'] else ''}; needs >= {d['need_case']} case, "
               f"<= {d['allow_control']} control detected") if d else "candidates only"
        print(f"  {c:<12} case={counts[c]['case']:<4} control={counts[c]['control']:<4} [{tag}]", file=sys.stderr)
    for cond, k in sorted(ext.items()):
        print(f"  external '{cond}': {k} (<= {cfg['design']['external'][cond]['allow']} may be detected)", file=sys.stderr)
    print(f"  candidates: >= {cfg['candidates']['min_count']} reads in >= "
          f"{cfg['design']['candidates_min_samples']} case/control samples (condition-blind)", file=sys.stderr)


if __name__ == "__main__":
    main()

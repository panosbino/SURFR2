#!/usr/bin/env python3
"""Synthetic small-RNA dataset with planted ground truth for SURFR2 regression testing.
Writes data/, samples.tsv and truth.tsv into the current directory. Needs samtools."""
import random, gzip, os, subprocess
random.seed(42)
os.makedirs("data", exist_ok=True)
rc = lambda s: s[::-1].translate(str.maketrans("ACGT","TGCA"))
rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))
base = [(rnd(22), random.lognormvariate(0, 1.5)) for _ in range(60)]   # shared "miRNAs"
ONCOMIR = rnd(22)                       # hugely up in cases -> composition effect
S = {k: rnd(25) for k in ["S1_true", "S2_in_controls", "S3_cohortA_only", "S4_in_external", "S5_one_control"]}
open("truth.tsv","w").write("".join(f"{k}\t{v}\n" for k,v in S.items()))
design = [("A","cancer",8),("A","adjacent_normal",6),("B","cancer",7),("B","adjacent_normal",6),("SRA","noncancer_sra",6)]
rows = ["sample_id\tcohort\tcondition\tfile_type\tpath"]
for coh, cond, n in design:
    for i in range(n):
        sid = f"{coh}_{cond}_{i}"
        depth = random.randint(20000, 100000)        # 5x depth range
        pool = [(s, w) for s, w in base] + [(ONCOMIR, 40.0 if cond=="cancer" else 2.0)]
        tot = sum(w for _, w in pool)
        reads = []
        for s, w in pool: reads += [s] * int(depth * w / tot)
        per1e5 = lambda r: int(depth * r / 1e5)
        case = cond == "cancer"
        if case and i < int(0.75 * n): reads += [S["S1_true"]] * per1e5(60)
        if case or cond == "adjacent_normal": reads += [S["S2_in_controls"]] * per1e5(60)
        if case and coh == "A": reads += [S["S3_cohortA_only"]] * per1e5(60)
        if case or cond == "noncancer_sra": reads += [S["S4_in_external"]] * per1e5(60)
        if case or (cond == "adjacent_normal" and i == 0): reads += [S["S5_one_control"]] * per1e5(60 if case else 20)
        random.shuffle(reads)
        fq = f"data/{sid}.fastq.gz"
        with gzip.open(fq, "wt") as f:
            for j, r in enumerate(reads): f.write(f"@{sid}_{j}\n{r}\n+\n{'I'*len(r)}\n")
        if coh == "B" and case:
            # aligned-style BAM, half reads on '-' strand (SEQ stored reverse-complemented)
            sam = f"data/{sid}.sam"
            with open(sam, "w") as f:
                f.write("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000000\n")
                for j, r in enumerate(reads):
                    if j % 2: f.write(f"{sid}_{j}\t16\tchr1\t100\t60\t{len(r)}M\t*\t0\t0\t{rc(r)}\t{'I'*len(r)}\n")
                    else:     f.write(f"{sid}_{j}\t0\tchr1\t100\t60\t{len(r)}M\t*\t0\t0\t{r}\t{'I'*len(r)}\n")
                # a secondary alignment that must NOT be counted
                f.write(f"{sid}_0\t256\tchr1\t500\t0\t{len(reads[0])}M\t*\t0\t0\t{reads[0]}\t{'I'*len(reads[0])}\n")
            subprocess.run(["samtools","view","-b","-o",f"data/{sid}.bam",sam], check=True); os.remove(sam); os.remove(fq)
            rows.append(f"{sid}\t{coh}\t{cond}\tbam\tdata/{sid}.bam")
        else:
            rows.append(f"{sid}\t{coh}\t{cond}\tfastq\tdata/{sid}.fastq.gz")
        open(f"data/{sid}.nreads","w").write(str(len(reads)))
open("samples.tsv","w").write("\n".join(rows)+"\n")

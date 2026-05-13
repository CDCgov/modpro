
#!usr/bin/env python3

"""
Usage notes:
snakemake \
  --snakefile Snakefile.openfold \
  --cores 16 \
  --printshellcmds \
  --reason \
  --rerun-incomplete \
  --keep-going \
  --latency-wait 60

conda env

```python
conda:
    "envs/openfold.yaml"
```
retries
```python
retries: 3
```

gpu usage
```bash
nvidia-smi | tee -a {log}
```

run metadata
```bash
date | tee -a {log}
hostname | tee -a {log}
env | sort | tee -a {log}
```
"""

import pandas as pd
from pathlib import Path

configfile: "config/config.yaml"

manifest = pd.read_csv(config["manifest"], sep="\t")
VARIANTS = manifest["variant_hash"].tolist()

FASTA_MAP = {
    row.variant_hash: row.fasta_path
    for _, row in manifest.iterrows()
}

Path("logs").mkdir(exist_ok=True)
Path("benchmarks").mkdir(exist_ok=True)

rule all:
    input:
        expand("work/archive/{variant}.tar.gz", variant=VARIANTS)

rule validate_input:
    input:
        fasta=lambda wc: FASTA_MAP[wc.variant]
    output:
        validated="work/validated/{variant}.validated.fasta"
    log:
        "logs/validate/{variant}.log"
    benchmark:
        "benchmarks/validate/{variant}.txt"
    threads: 2
    shell:
        r"""
        set -euo pipefail

        mkdir -p work/validated logs/validate benchmarks/validate

        echo "[VALIDATE] Starting validation for {wildcards.variant}" | tee {log}

        python scripts/validate_input.py \
            --input {input.fasta} \
            --output {output.validated} \
            2>&1 | tee -a {log}

        echo "[VALIDATE] Finished validation" | tee -a {log}
        """

# MSA generation

rule generate_msa:
    input:
        fasta="work/validated/{variant}.validated.fasta"
    output:
        msa="work/msa/{variant}.a3m"
    log:
        "logs/msa/{variant}.log"
    benchmark:
        "benchmarks/msa/{variant}.txt"
    threads:
        config["msa"]["threads"]
    resources:
        mem_mb=config["resources"]["msa_mem_mb"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p work/msa logs/msa benchmarks/msa

        echo "[MSA] Starting MSA generation for {wildcards.variant}" | tee {log}

        python scripts/generate_msa.py \
            --input {input.fasta} \
            --db {config[msa][db]} \
            --threads {threads} \
            --output {output.msa} \
            2>&1 | tee -a {log}

        echo "[MSA] Completed MSA generation" | tee -a {log}
        """

# Extract metadata

rule extract_hash_chain:
    input:
        msa="work/msa/{variant}.a3m"
    output:
        metadata="work/hash_chain/{variant}.json"
    log:
        "logs/hash_chain/{variant}.log"
    benchmark:
        "benchmarks/hash_chain/{variant}.txt"
    shell:
        r"""
        set -euo pipefail

        mkdir -p work/hash_chain logs/hash_chain benchmarks/hash_chain

        echo "[HASH_CHAIN] Extracting metadata for {wildcards.variant}" | tee {log}

        python scripts/extract_hash_chain.py \
            --msa {input.msa} \
            --output {output.metadata} \
            2>&1 | tee -a {log}

        echo "[HASH_CHAIN] Extraction complete" | tee -a {log}
        """

# run openfold3 inference

rule infer_structure:
    input:
        msa="work/msa/{variant}.a3m",
        metadata="work/hash_chain/{variant}.json"

    output:
        pdb="work/predicted/{variant}.pdb"

    log:
        "logs/infer/{variant}.log"

    benchmark:
        "benchmarks/infer/{variant}.txt"

    threads: 16

    resources:
        gpu=config["openfold"]["gpu"],
        mem_mb=config["resources"]["infer_mem_mb"]

    shell:
        r"""
        set -euo pipefail

        mkdir -p work/predicted logs/infer benchmarks/infer

        echo "[INFER] Starting structure prediction for {wildcards.variant}" | tee {log}

        python scripts/infer_structure.py \
            --msa {input.msa} \
            --metadata {input.metadata} \
            --checkpoint {config[openfold][checkpoint]} \
            --output {output.pdb} \
            2>&1 | tee -a {log}

        echo "[INFER] Structure inference complete" | tee -a {log}
        """

# QC reports

rule generate_report:
    input:
        pdb="work/predicted/{variant}.pdb"

    output:
        report="work/reports/{variant}.qc.json"

    log:
        "logs/report/{variant}.log"

    benchmark:
        "benchmarks/report/{variant}.txt"

    shell:
        r"""
        set -euo pipefail

        mkdir -p work/reports logs/report benchmarks/report

        echo "[REPORT] Generating QC report for {wildcards.variant}" | tee {log}

        python scripts/generate_report.py \
            --pdb {input.pdb} \
            --output {output.report} \
            2>&1 | tee -a {log}

        echo "[REPORT] QC generation complete" | tee -a {log}
        """

# archive and upload hdfs

rule archive_results:
    input:
        pdb="work/predicted/{variant}.pdb",
        report="work/reports/{variant}.qc.json"

    output:
        archive="work/archive/{variant}.tar.gz"

    log:
        "logs/archive/{variant}.log"

    benchmark:
        "benchmarks/archive/{variant}.txt"

    shell:
        r"""
        set -euo pipefail

        mkdir -p work/archive logs/archive benchmarks/archive

        echo "[ARCHIVE] Archiving outputs for {wildcards.variant}" | tee {log}

        tar -czf {output.archive} \
            {input.pdb} \
            {input.report} \
            2>&1 | tee -a {log}

        echo "[ARCHIVE] Archive complete" | tee -a {log}


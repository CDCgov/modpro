
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

MANIFEST_DIR = "manifests"

OPENFOLD3_ROOT = "/scicomp/home-pure/qxa4/openfold3"
OPENFOLD3_SCRIPTS = f"{OPENFOLD3_ROOT}/openfold-3/scripts"
OPENFOLD3_MSA_DIR = f"{OPENFOLD3_SCRIPTS}/snakemake_msa"
OPENFOLD3_MSA_SNAKEFILE = f"{OPENFOLD3_MSA_DIR}/Snakefile"

OPENFOLD_ENV = "/scicomp/groups/OID/NCIRD/ID/VSDB/GAT/shared_conda_envs/of3-aln-env"
OPENFOLD_DB_PATH = "/scicomp/groups/OID/NCIRD/ID/VSDB/GAT/of3_dbs"

RUNNER_YAML = "config/runner.yaml"

JACKHMMER_THREADS = 16
HHBLITS_THREADS = 32
RUNNER_YAML = "config/runner.yaml"

HA_CHAIN_IDS = ["A", "B", "C"]   # trimer
NA_CHAIN_IDS = ["A", "B", "C", "D"]   # tetramer
QC_THRESHOLD = 70.0

#one rule to rule them all

rule all:
    input:
        []


# manifest to fasta
rule normalize_manifest_to_fasta:
    input:
        manifest = f"{MANIFEST_DIR}/{{manifest}}.csv"
    output:
        fasta = "results/{manifest}/seq/seq.fasta",
        seq_name_json = "results/{manifest}/config/seq_name.json"
    run:
        Path(f"results/{wildcards.manifest}/seq").mkdir(parents=True, exist_ok=True)
        Path(f"results/{wildcards.manifest}/config").mkdir(parents=True, exist_ok=True)
        lines = [
            line.strip()
            for line in open(input.manifest)
            if line.strip()
        ]
        records = []
        i = 0
        while i < len(lines):
            header = lines[i]
            if header.startswith(">"):
                name = header[1:].strip()
            else:
                name = header.strip()
            if i + 1 >= len(lines):
                raise ValueError(
                    f"Manifest {input.manifest} has hash/header without sequence: {name}"
                )
            seq = lines[i + 1].strip()
            records.append((name, seq))
            i += 2
        if not records:
            raise ValueError(f"No sequence records found in {input.manifest}")
        with open(output.fasta, "w") as out:
            for name, seq in records:
                out.write(f">{name}\n")
                out.write(f"{seq}\n")
        with open(output.seq_name_json, "w") as out:
            json.dump(
                [
                    {
                        "name": name,
                        "sequence_length": len(seq)
                    }
                    for name, seq in records
                ],
                out,
                indent=2
            )
 
# msa configs
rule create_openfold3_msa_gen_json:
    input:
        fasta = "results/{manifest}/seq/seq.fasta"
    output:
        msa_gen_json = f"{OPENFOLD3_MSA_DIR}/{{manifest}}_msa_gen.json"
    run:
        Path(OPENFOLD3_MSA_DIR).mkdir(parents=True, exist_ok=True)
        cfg = {
            "input_fasta": str(Path(input.fasta).resolve()),
            "openfold_env": OPENFOLD_ENV,
            "databases": ["uniref90", "uniprot", "mgnify", "bfd"],
            "base_database_path": OPENFOLD_DB_PATH,
            "output_directory": str(Path(f"results/{wildcards.manifest}/msa").resolve()),
            "jackhmmer_output_format": "sto",
            "jackhmmer_threads": JACKHMMER_THREADS,
            "hhblits_threads": HHBLITS_THREADS,
            "tmpdir": "/tmp",
            "run_template_search": False
        }
        with open(output.msa_gen_json, "w") as out:
            json.dump(cfg, out, indent=4)

rule create_openfold3_prediction_config:
    input:
        fasta = "results/{manifest}/seq/seq.fasta"
    output:
        config_json = "results/{manifest}/config/openfold3_input.json"
    run:
        Path(f"results/{wildcards.manifest}/config").mkdir(parents=True, exist_ok=True)
        records = []
        name = None
        seq_parts = []
        for line in open(input.fasta):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_parts)))
                name = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)
        if name is not None:
            records.append((name, "".join(seq_parts)))
        if not records:
            raise ValueError(f"No FASTA records found in {input.fasta}")
        of3_records = []
        for name, seq in records:
            msa_dir = str(Path(f"results/{wildcards.manifest}/msa").resolve())
            of3_records.append({
                "molecule_type": "protein",
                "chain_ids": "A",
                "sequence": seq,
                "use_msas": True,
                "use_main_msas": True,
                "use_paired_msas": False,
                "main_msa_file_paths": msa_dir,
                "paired_msa_file_paths": msa_dir,
                "template_alignment_file_path": msa_dir
            })
        with open(output.config_json, "w") as out:
            if len(of3_records) == 1:
                json.dump(of3_records[0], out, indent=2)
            else:
                json.dump(of3_records, out, indent=2)

# run msa

rule run_openfold3_msas:
    input:
        msa_gen_json = f"{OPENFOLD3_MSA_DIR}/{{manifest}}_msa_gen.json"
    output:
        done = "results/{manifest}/msa/msa.done"
    shell:
        """
        mkdir -p results/{wildcards.manifest}/msa

        cd {OPENFOLD3_MSA_DIR}

        snakemake \
          --configfile {input.msa_gen_json} \
          --cores {HHBLITS_THREADS} \
          --rerun-incomplete

        touch {output.done}
        """

# inference
rule run_openfold3_prediction:
    input:
        query_json = "results/{manifest}/openfold3/query.json"
    output:
        done = "results/{manifest}/openfold3/openfold_output/prediction.done"
    params:
        output_dir = "results/{manifest}/openfold3/openfold_output",
        runner_yaml = RUNNER_YAML
    shell:
        """
        mkdir -p {params.output_dir}

        run_openfold predict \
          --query_json {input.query_json} \
          --use_msa_server=False \
          --output_dir {params.output_dir} \
          --runner_yaml {params.runner_yaml}

        touch {output.done}
        """
# qc

rule check_qc:
    input:
        plddt = "results/{manifest}/reports/plddt.txt"
    output:
        flag = "results/{manifest}/qc/passed.flag"
    params:
        threshold = 70.0
    run:
        Path(f"results/{wildcards.manifest}/qc").mkdir(parents=True, exist_ok=True)

        score = float(open(input.plddt).read().strip())

        if score < params.threshold:
            raise ValueError(
                f"QC threshold not met for {wildcards.manifest}: {score}"
            )
        shell("touch {output.flag}")

# ha1 extraction

rule extract_ha1:
    input:
        pdb = "results/{manifest}/prediction/model.pdb",
        qc = "results/{manifest}/qc/passed.flag"
    output:
        ha1 = "results/{manifest}/ha1/ha1_only.pdb"
    shell:
        """
        mkdir -p results/{wildcards.manifest}/ha1
        python scripts/extract_ha1.py {input.pdb} > {output.ha1}
        """

rule archive_full_pdb:
    input:
        pdb = "results/{manifest}/prediction/model.pdb"
    output:
        archived = "results/{manifest}/archive/full_length.pdb"
    shell:
        """
        mkdir -p results/{wildcards.manifest}/archive
        cp {input.pdb} {output.archived}
        """

rule rosetta_score_ha1:
    input:
        ha1 = "results/{manifest}/ha1/ha1_only.pdb"
    output:
        score = "results/{manifest}/ha1/score.sc"
    shell:
        """
        rosetta_score.linuxgccrelease \
          -s {input.ha1} \
          -out:file:scorefile {output.score}
        """

rule rmsd_maxsub_ha1:
    input:
        ha1 = "results/{manifest}/ha1/ha1_only.pdb"
    output:
        metrics = "results/{manifest}/ha1/rmsd_maxsub.txt"
    shell:
        """
        python scripts/calc_rmsd_maxsub.py {input.ha1} > {output.metrics}
        """

rule rosetta_relax:
    input:
        ha1 = "results/{manifest}/ha1/ha1_only.pdb"
    output:
        relaxed = "results/{manifest}/ha1/relaxed.pdb"
    shell:
        """
        rosetta_relax.linuxgccrelease \
          -s {input.ha1} \
          -out:path:pdb results/{wildcards.manifest}/ha1

        mv results/{wildcards.manifest}/ha1/*relaxed*.pdb {output.relaxed}
        """

rule select_best_ha1:
    input:
        original = "results/{manifest}/ha1/ha1_only.pdb",
        relaxed = "results/{manifest}/ha1/relaxed.pdb",
        score = "results/{manifest}/ha1/score.sc",
        metrics = "results/{manifest}/ha1/rmsd_maxsub.txt",
        archived = "results/{manifest}/archive/full_length.pdb"
    output:
        selected = "results/{manifest}/final/ha1_selected.pdb"
    shell:
        """
        mkdir -p results/{wildcards.manifest}/final

        python scripts/select_best.py \
          --original {input.original} \
          --relaxed {input.relaxed} \
          --score {input.score} \
          --metrics {input.metrics} \
          --out {output.selected}
        """

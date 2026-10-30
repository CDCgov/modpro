#!/usr/bin/env python3
"""
Build the protein_modeling.logging table (schema version 1) for one run.

One row per variant in the run manifest. Each step column is True when the
files that step produces exist and are non-empty, False otherwise. `sia_dock`
is not implemented yet and is written as null (nullable boolean); it is not
part of the status decision until it is implemented.

status = "SUCCESS" when every implemented step column is True, else "FAIL".

Schema (v1):
    variant_hash         string
    run_date             timestamp
    msa_done             boolean
    structure_done       boolean
    rosetta_relax_done   boolean
    rosetta_score        boolean
    glyc_distance        boolean
    structural_distance  boolean
    epitope_calc         boolean
    sia_dock             boolean (null until implemented)
    status               string

Usage (normally invoked by Snakefile-thursdays):

    python pipeline_status_table.py \
        --run-id 20261027_b6edee \
        --variants-json results/20261027_b6edee/inputs/variants.json \
        --step 'msa_done=results/{run_id}/openfold3/{variant_hash}/msa/msa.done' \
        --step 'structure_done=results/{run_id}/openfold3/{variant_hash}/prediction/..._model.pdb' \
        --step 'rosetta_score=path/one|path/two' \
        ... \
        --output results/logging_qc/20261027_b6edee/protein_modeling_logging_v1.parquet \
        --tsv    results/logging_qc/20261027_b6edee/protein_modeling_logging_v1.tsv \
        --staging-dir /path/to/staging

Templates may use {run_id} and {variant_hash}. Several templates for one step
are separated by "|" and must all exist for the step to count as done.
Always exits 0 when the table was written; the status is data, not an error.
"""
import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

SCHEMA_VERSION = "1"
TABLE_NAME = "protein_modeling.logging"

STEP_COLUMNS = [
    "msa_done",
    "structure_done",
    "rosetta_relax_done",
    "rosetta_score",
    "glyc_distance",
    "structural_distance",
    "epitope_calc",
    "sia_dock",
]
NOT_IMPLEMENTED = {"sia_dock"}
IMPLEMENTED_STEPS = [c for c in STEP_COLUMNS if c not in NOT_IMPLEMENTED]
COLUMNS = ["variant_hash", "run_date"] + STEP_COLUMNS + ["status"]


def run_date_from_id(run_id):
    text = run_id[:8]
    if len(text) != 8 or not text.isdigit():
        raise ValueError(f"run_id must start with YYYYMMDD: {run_id}")
    return datetime.strptime(text, "%Y%m%d")


def parse_steps(step_args):
    steps = {}
    for item in step_args:
        if "=" not in item:
            raise ValueError(f"--step must be NAME=TEMPLATE[|TEMPLATE...]: {item}")
        name, templates = item.split("=", 1)
        name = name.strip()
        if name not in STEP_COLUMNS:
            raise ValueError(f"Unknown step {name!r}; expected one of {STEP_COLUMNS}")
        steps[name] = [t for t in templates.split("|") if t.strip()]
    return steps


def file_ok(path):
    """True when the file exists and carries evidence the step ran.

    A ".done" marker is created with `touch` and is legitimately ZERO bytes,
    so for markers existence is the whole test. Everything else (parquet,
    pdb, tsv, Rosetta .sc) must be non-empty, which is what catches a step
    that created its output file and then died before writing to it.

    Requiring non-empty for markers was failing msa_done and rosetta_score
    for every variant on 2026-10-27 while the pipeline was in fact fine.
    """
    p = Path(path)
    if not p.is_file():
        return False
    if p.suffix == ".done":
        return True
    return p.stat().st_size > 0


def step_done(templates, run_id, variant_hash):
    if not templates:
        return False
    return all(
        file_ok(t.format(run_id=run_id, variant_hash=variant_hash))
        for t in templates
    )


def build_rows(run_id, variants, steps):
    run_date = run_date_from_id(run_id)
    rows = []
    for v in variants:
        row = {"variant_hash": v, "run_date": run_date}
        for col in STEP_COLUMNS:
            if col in NOT_IMPLEMENTED:
                row[col] = None
            else:
                row[col] = step_done(steps.get(col, []), run_id, v)
        row["status"] = "SUCCESS" if all(row[c] for c in IMPLEMENTED_STEPS) else "FAIL"
        rows.append(row)
    return rows


def to_frame(rows):
    df = pd.DataFrame(rows, columns=COLUMNS)
    df["variant_hash"] = df["variant_hash"].astype("string")
    df["run_date"] = pd.to_datetime(df["run_date"])
    for col in IMPLEMENTED_STEPS:
        df[col] = df[col].astype("bool")
    for col in NOT_IMPLEMENTED:
        df[col] = pd.array([pd.NA] * len(df), dtype="boolean")
    df["status"] = df["status"].astype("string")
    return df


def write_parquet(df, path):
    table = pa.Table.from_pandas(df, preserve_index=False)
    meta = dict(table.schema.metadata or {})
    meta[b"table"] = TABLE_NAME.encode()
    meta[b"schema_version"] = SCHEMA_VERSION.encode()
    meta[b"generated_at"] = datetime.now().isoformat().encode()
    table = table.replace_schema_metadata(meta)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, path)


def summary(df, run_id):
    n = len(df)
    n_ok = int((df["status"] == "SUCCESS").sum()) if n else 0
    lines = [f"{TABLE_NAME} v{SCHEMA_VERSION} run_id={run_id}: {n_ok}/{n} variants SUCCESS"]
    for _, r in df[df["status"] != "SUCCESS"].iterrows():
        missing = [c for c in IMPLEMENTED_STEPS if not bool(r[c])]
        lines.append(f"  {r['variant_hash']}  FAIL  missing: {', '.join(missing)}")
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--variants-json", required=True,
                    help="variants.json written by the load_run_manifest checkpoint")
    ap.add_argument("--step", action="append", default=[],
                    help="NAME=TEMPLATE[|TEMPLATE...]; repeatable")
    ap.add_argument("--output", required=True, help="parquet path")
    ap.add_argument("--tsv", help="optional TSV copy (used by thursday_complete)")
    ap.add_argument("--staging-dir", help="optional directory to copy the parquet into")
    args = ap.parse_args(argv)

    with open(args.variants_json) as fh:
        variants = [r["variant_hash"] for r in json.load(fh)]
    steps = parse_steps(args.step)
    missing_defs = [c for c in IMPLEMENTED_STEPS if c not in steps]
    if missing_defs:
        print(f"WARNING: no --step given for {missing_defs}; they will be False", file=sys.stderr)

    df = to_frame(build_rows(args.run_id, variants, steps))
    write_parquet(df, args.output)
    if args.tsv:
        Path(args.tsv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.tsv, sep="\t", index=False, date_format="%Y-%m-%d")
    if args.staging_dir:
        staging = Path(args.staging_dir)
        staging.mkdir(parents=True, exist_ok=True)
        (staging / Path(args.output).name).write_bytes(Path(args.output).read_bytes())

    print(summary(df, args.run_id))
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())


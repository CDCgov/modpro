#!/usr/bin/env python3

import argparse
import csv
import smtplib
import socket
import sys
from pathlib import Path
from email.message import EmailMessage
import yaml

CONFIG_PATH = Path(__file__).resolve().parents[1] / "config.yaml"

with open(CONFIG_PATH) as handle:
    config = yaml.safe_load(handle)

EMAIL_TO = config["email"]["to"]
EMAIL_FROM = config["email"]["from"]

ROSETTA_READY_DIR = Path(
    config["rosetta"]["outputs"]["rosetta_ready_dir"]
)

HA1_BOUNDS_BY_SUBTYPE = {
    "H1": (18, 344),
    "H3": (17, 345),
    "Bvic": (16, 362),
    "H5": None,
}


def normalize_subtype(value):
    if value is None:
        return "UNKNOWN"

    s = (
        value.strip()
        .lower()
        .replace("_", "")
        .replace("-", "")
        .replace(" ", "")
    )

    if s.startswith("h1"):
        return "H1"

    if s.startswith("h3"):
        return "H3"

    if s.startswith("h5"):
        return "H5"

    if s in {
        "b",
        "bvic",
        "bvictoria",
        "b/vic",
        "b/victoria",
        "hvic",
        "hvictoria",
        "h/vic",
        "h/victoria",
        "victoria",
        "vic",
    }:
        return "Bvic"

    return "UNKNOWN"


def send_failure_email(subject, body):
    msg = EmailMessage()
    msg["To"] = EMAIL_TO
    msg["From"] = EMAIL_FROM
    msg["Subject"] = subject
    msg.set_content(body)

    with smtplib.SMTP("localhost") as smtp:
        smtp.send_message(msg)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract chain A and trim HA2/signal peptide "
            "using subtype-specific HA1 residue bounds."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input OpenFold/AlphaFold PDB"
    )

    parser.add_argument(
        "--metadata_path",
        required=True,
        help="Run metadata file"
    )

    parser.add_argument(
        "--variant",
        required=True,
        help="Variant hash"
    )

    parser.add_argument(
        "--run",
        required=True,
        help="Run ID"
    )

    parser.add_argument(
        "--output",
        required=False,
        help="Output HA1 PDB"
    )

    return parser.parse_args()


def residue_number(line):
    try:
        return int(line[22:26].strip())
    except ValueError:
        return None


def chain_id(line):
    return line[21].strip()


def read_metadata(metadata_path, variant_hash):
    metadata_path = Path(metadata_path)

    if not metadata_path.exists():
        raise FileNotFoundError(
            f"Metadata file does not exist: {metadata_path}"
        )

    with metadata_path.open(newline="") as handle:
        reader = csv.DictReader(handle)

        required = {"variant_hash", "protein", "subtype"}

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                f"Metadata file missing required columns "
                f"{sorted(missing)}. "
                f"Found columns: {reader.fieldnames}"
            )

        for row in reader:
            if row.get("variant_hash", "").strip() == variant_hash:

                protein = (
                    row.get("protein", "")
                    .strip()
                    .upper()
                )

                raw_subtype = (
                    row.get("subtype", "")
                    .strip()
                )

                subtype = normalize_subtype(raw_subtype)

                return {
                    "protein": protein,
                    "raw_subtype": raw_subtype,
                    "subtype": subtype,
                }

    raise ValueError(
        f"Variant {variant_hash} not found in "
        f"metadata file: {metadata_path}"
    )


def trim_ha_chain_a(input_pdb, output_pdb, subtype):

    if subtype not in HA1_BOUNDS_BY_SUBTYPE:
        raise ValueError(
            f"Unsupported normalized subtype: {subtype}"
        )

    bounds = HA1_BOUNDS_BY_SUBTYPE[subtype]

    if bounds is None:
        raise ValueError(
            f"HA1 bounds not defined for subtype {subtype}"
        )

    ha1_start, ha1_end = bounds

    input_pdb = Path(input_pdb)
    output_pdb = Path(output_pdb)

    output_pdb.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    if not input_pdb.exists():
        raise FileNotFoundError(
            f"Input PDB does not exist: {input_pdb}"
        )

    kept_atoms = 0
    skipped_wrong_chain = 0
    skipped_outside_bounds = 0

    with input_pdb.open() as inp, output_pdb.open("w") as out:

        out.write(f"REMARK source_pdb={input_pdb}\n")
        out.write("REMARK chain=A only\n")
        out.write(f"REMARK subtype={subtype}\n")
        out.write(
            f"REMARK HA1_bounds={ha1_start}..{ha1_end}\n"
        )
        out.write(
            "REMARK residue numbering preserved\n"
        )

        for line in inp:

            record = line[:6].strip()

            if record not in {"ATOM", "HETATM"}:
                continue

            if chain_id(line) != "A":
                skipped_wrong_chain += 1
                continue

            resnum = residue_number(line)

            if resnum is None:
                continue

            if not (ha1_start <= resnum <= ha1_end):
                skipped_outside_bounds += 1
                continue

            out.write(line)
            kept_atoms += 1

        out.write("TER\n")
        out.write("END\n")

    if kept_atoms == 0:
        raise RuntimeError(
            f"No atoms written to {output_pdb}. "
            f"Check chain IDs and residue numbering."
        )

    return {
        "kept_atoms": kept_atoms,
        "skipped_wrong_chain": skipped_wrong_chain,
        "skipped_outside_bounds": skipped_outside_bounds,
        "ha1_start": ha1_start,
        "ha1_end": ha1_end,
    }


def main():

    args = parse_args()

    input_pdb = Path(args.input)
    metadata_path = Path(args.metadata_path)

    output_pdb = (
        Path(args.output)
        if args.output
        else (
            ROSETTA_READY_DIR
            / args.run
            / f"{args.variant}.pdb"
        )
    )

    try:

        meta = read_metadata(
            metadata_path,
            args.variant
        )

        print(f"host={socket.gethostname()}")
        print(f"run={args.run}")
        print(f"variant={args.variant}")
        print(f"input_pdb={input_pdb}")
        print(f"metadata={metadata_path}")
        print(f"protein={meta['protein']}")
        print(f"raw_subtype={meta['raw_subtype']}")
        print(f"normalized_subtype={meta['subtype']}")
        print(f"output_pdb={output_pdb}")

        if meta["protein"] != "HA":
            raise RuntimeError(
                f"Variant {args.variant} "
                f"has protein={meta['protein']}; "
                "this script only handles HA."
            )

        stats = trim_ha_chain_a(
            input_pdb=input_pdb,
            output_pdb=output_pdb,
            subtype=meta["subtype"],
        )

        print("trim_complete=true")

        for key, value in stats.items():
            print(f"{key}={value}")

    except Exception as exc:

        subject = (
            f"[modpro] HA1 trim failed "
            f"for {args.variant}"
        )

        body = f"""
HA1 trim failed.

Host:
{socket.gethostname()}

Run:
{args.run}

Variant:
{args.variant}

Input:
{input_pdb}

Metadata:
{metadata_path}

Output:
{output_pdb}

Error:
{repr(exc)}
"""

        try:
            send_failure_email(subject, body)
        except Exception as email_exc:
            print(
                f"Failed to send failure email: "
                f"{email_exc}",
                file=sys.stderr,
            )

        raise


if __name__ == "__main__":
    main()

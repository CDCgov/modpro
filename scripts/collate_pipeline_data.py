#!/usr/bin/env python3
"""
Collates pipeline data into a parquet file for SQL upload.
Extracts metrics from various pipeline output files.

Fixes applied 2026-10-27:

1. CRASH: main() called collate_variant_data() with four positional arguments
   while the function accepted three ("takes 3 positional arguments but 4
   were given"). results_dir is now an accepted, documented parameter. It is
   still unused; the reference scores all live under logging_qc.

2. run_date was datetime.now(), i.e. the moment of collation, not the date of
   the run. It is now derived from the first eight characters of run_id, the
   same rule the Snakefiles use, so collated_data.run_date matches
   protein_modeling.logging.run_date. Falls back to now() if run_id is not
   dated.

3. The reference .sc file was picked with an unsorted glob, so which
   reference the rms/maxsub came from varied with filesystem order.

4. GRAIN CHANGE: one row per variant-reference pair.
   compare_to_reference scores each variant against every reference PDB
   (about six of them), but the script kept only one and discarded the rest.
   It now emits one row per reference, sorted by reference name within each
   variant, so a run of 11 variants against 6 references gives 66 rows
   instead of 11.

   Columns are unchanged. variant_hash, run_date, relaxed_score,
   unrelaxed_score and avg_plddt are per-variant and repeat across that
   variant's rows; reference, rms and maxsub vary.

   Downstream consumers that assumed one row per variant need updating: the
   natural key is now (variant_hash, run_date, reference). To recover the old
   shape, group by variant_hash and pick a reference, or aggregate rms.
"""

import argparse
import json
import re
from pathlib import Path
from datetime import datetime
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import shutil


def run_date_from_id(run_id):
    """YYYY-MM-DD datetime from the first eight characters of run_id.

    Returns None when run_id is not date-prefixed, so the caller can fall
    back rather than crash.
    """
    date_text = str(run_id)[:8]
    if len(date_text) == 8 and date_text.isdigit():
        try:
            return datetime.strptime(date_text, "%Y%m%d")
        except ValueError:
            return None
    return None


def extract_value_from_log(log_file, key_pattern):
    """Extract a value from a log file using a regex pattern."""
    try:
        with open(log_file, 'r') as f:
            content = f.read()
            match = re.search(key_pattern, content)
            if match:
                return float(match.group(1))
    except Exception as e:
        print(f"Warning: Could not extract from {log_file}: {e}")
    return None


def extract_score_jd2_value(score_file, field_name):
    """Extract a named value from a Rosetta score_jd2.sc file."""
    try:
        with open(score_file, 'r') as f:
            lines = f.readlines()
            # Find the SCORE header line
            header_line = None
            score_line = None
            for line in lines:
                if line.startswith('SCORE:'):
                    if 'total_score' in line:
                        header_line = line
                    else:
                        score_line = line
                        break

            if header_line and score_line:
                # Parse header to find column index
                header_parts = header_line.split()
                try:
                    col_idx = header_parts.index(field_name)
                    score_parts = score_line.split()
                    return float(score_parts[col_idx])
                except (ValueError, IndexError) as e:
                    print(f"Warning: Could not find {field_name} in {score_file}: {e}")
    except Exception as e:
        print(f"Warning: Could not parse {score_file}: {e}")
    return None


def extract_plddt_avg(plddt_log):
    """Extract average pLDDT from plddt.log file."""
    try:
        with open(plddt_log, 'r') as f:
            content = f.read()
            # Look for avg=XX.X pattern
            match = re.search(r'avg=([0-9.]+)', content)
            if match:
                return float(match.group(1))
    except Exception as e:
        print(f"Warning: Could not extract pLDDT avg from {plddt_log}: {e}")
    return None


def collate_variant_data(run_id, variant_hash, logging_qc_dir, results_dir=None):
    """
    Collate all data for a single variant, one row per reference.

    compare_to_reference scores the variant against every reference PDB and
    writes one <ref_name>.score_jd2.sc per reference. This returns one row
    for each of them, so no comparison is discarded.

    Args:
        run_id: The run identifier
        variant_hash: The variant hash
        logging_qc_dir: Path to logging_qc directory
        results_dir: Path to the results directory. Accepted because the
            caller passes it; not currently read, since every metric below
            comes from logging_qc. Keeping the parameter is what stops the
            "takes 3 positional arguments but 4 were given" crash.

    Returns:
        List of row dicts, one per reference. Empty list when the variant has
        no usable reference data, which the caller reports as skipped.
    """
    logging_qc_path = Path(logging_qc_dir) / run_id

    # compare_to_reference writes one <ref_name>.score_jd2.sc per reference
    # into the directory <variant_hash>.reference_compare.
    reference_dir = logging_qc_path / f"{variant_hash}.reference_compare"

    if reference_dir.is_dir():
        # Sorted so row order is stable across runs and machines.
        score_files = sorted(reference_dir.glob("*.sc"))
        if not score_files:
            print(f"Warning: No .sc files found in reference_compare directory {reference_dir}")
            return []
    elif reference_dir.exists():
        # Older layout: a single file rather than a directory.
        score_files = [reference_dir]
    else:
        print(f"Warning: No reference_compare file/directory found at {reference_dir}")
        return []

    references_data = []
    for score_file in score_files:
        try:
            with open(score_file, 'r') as f:
                lines = f.readlines()
        except Exception as e:
            print(f"Warning: Could not read reference score file {score_file}: {e}")
            continue
        # Require at least one data row beyond the header.
        if not [l for l in lines[1:] if len(l.strip().split()) >= 2]:
            print(f"Warning: No data rows in {score_file}")
            continue
        references_data.append({
            'reference': score_file.stem.replace('.score_jd2', ''),
            'rms': extract_score_jd2_value(str(score_file), 'rms'),
            'maxsub': extract_score_jd2_value(str(score_file), 'maxsub'),
        })

    if not references_data:
        print(f"Warning: No valid reference data extracted for {variant_hash}")
        return []

    # ---- per-variant metrics, read once and repeated on every row ----

    # Find best_energy decision TSV to get relaxed_score and original_score
    decision_file = logging_qc_path / f"{variant_hash}_best_energy_ha1.tsv"
    relaxed_score = None
    original_score = None

    if decision_file.exists():
        try:
            with open(decision_file, 'r') as f:
                lines = f.readlines()
                if len(lines) >= 2:
                    # Parse header
                    header = lines[0].strip().split('\t')
                    # Parse data row
                    data_row = lines[1].strip().split('\t')

                    # Find column indices
                    try:
                        relaxed_score_idx = header.index('relaxed_score')
                        original_score_idx = header.index('original_score')

                        if relaxed_score_idx < len(data_row):
                            relaxed_score = float(data_row[relaxed_score_idx])
                        if original_score_idx < len(data_row):
                            original_score = float(data_row[original_score_idx])
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Could not parse decision file columns: {e}")
                else:
                    print(f"Warning: Decision file has fewer than 2 lines: {decision_file}")
        except Exception as e:
            print(f"Warning: Could not read decision file {decision_file}: {e}")
    else:
        print(f"Warning: Decision file not found at {decision_file}")

    # Find plddt log to get average pLDDT
    plddt_log = logging_qc_path / f"{variant_hash}.plddt.log"
    avg_plddt = None

    if plddt_log.exists():
        avg = extract_plddt_avg(plddt_log)
        if avg is not None:
            avg_plddt = avg
    else:
        print(f"Warning: pLDDT log not found at {plddt_log}")

    # run_date is the date of the RUN, taken from run_id, not the time this
    # script happens to run. Keeps collated_data aligned with
    # protein_modeling.logging.
    run_date = run_date_from_id(run_id)
    if run_date is None:
        run_date = datetime.now()
        print(f"Warning: run_id {run_id} is not date-prefixed; using current time as run_date")

    return [
        {
            'variant_hash': variant_hash,
            'run_date': run_date,
            'reference': ref['reference'],
            'rms': ref['rms'],
            'maxsub': ref['maxsub'],
            'relaxed_score': relaxed_score,
            'unrelaxed_score': original_score,
            'avg_plddt': avg_plddt,
        }
        for ref in references_data
    ]


def main():
    parser = argparse.ArgumentParser(
        description='Collate pipeline data into parquet file for SQL upload'
    )
    parser.add_argument('--run-id', required=True, help='Run ID')
    parser.add_argument('--variant-hashes', required=True, help='Comma-separated variant hashes')
    parser.add_argument('--logging-qc-dir', required=True, help='Path to logging_qc directory')
    parser.add_argument('--results-dir', required=True, help='Path to results directory')
    parser.add_argument('--staging-dir', required=True, help='Path to staging directory')
    parser.add_argument('--output', required=True, help='Output parquet file path (will copy to staging)')

    args = parser.parse_args()

    variant_hashes = [v.strip() for v in args.variant_hashes.split(',') if v.strip()]
    all_data = []
    skipped = []

    for variant_hash in variant_hashes:
        rows = collate_variant_data(
            args.run_id,
            variant_hash,
            args.logging_qc_dir,
            args.results_dir
        )
        if rows:
            all_data.extend(rows)
        else:
            skipped.append(variant_hash)
            print(f"Skipping variant {variant_hash} due to missing data")

    if not all_data:
        print(
            f"Error: No valid data collected from any of {len(variant_hashes)} variants. "
            f"Looked under {Path(args.logging_qc_dir) / args.run_id}"
        )
        return 1

    # Create DataFrame and save to parquet without header
    df = pd.DataFrame(all_data)

    # Ensure column order
    column_order = [
        'variant_hash', 'run_date', 'reference', 'rms', 'maxsub',
        'relaxed_score', 'unrelaxed_score', 'avg_plddt'
    ]
    df = df[column_order]

    # Write to parquet file
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write parquet file in Impala-compatible format (v1.0)
    table = pa.Table.from_pandas(df, preserve_index=False)
    pq.write_table(
        table,
        args.output,
        version='1.0',
        coerce_timestamps='ms'
    )

    # Copy to staging directory
    staging_dir = Path(args.staging_dir)
    staging_dir.mkdir(parents=True, exist_ok=True)
    staging_path = staging_dir / output_path.name
    shutil.copy(args.output, staging_path)

    print(f"Successfully wrote {len(df)} records to {args.output}")
    print(
        f"  {df['variant_hash'].nunique()} variants x "
        f"{df['reference'].nunique()} references"
    )
    if skipped:
        print(f"Skipped {len(skipped)} variant(s) with missing data: {', '.join(skipped)}")
    print(f"Copied to staging: {staging_path}")

    return 0


if __name__ == '__main__':
    exit(main())


#!/usr/bin/env python3

from pathlib import Path
import subprocess
import time
import hashlib

MANIFEST_DIR = Path("./manifests")
DONE_DIR = MANIFEST_DIR / ".processed"
DONE_DIR.mkdir(exist_ok=True)

POLL_SECONDS = 10

def file_is_stable(path: Path, wait_seconds: int = 5) -> bool:
    size1 = path.stat().st_size
    time.sleep(wait_seconds)
    size2 = path.stat().st_size
    return size1 == size2

def manifest_id(path: Path) -> str:
    return path.stem

def run_snakemake(path: Path):
    mid = path.stem
    target = f"results/{mid}/final/ha1_selected.pdb"

    cmd = [
        "snakemake",
        "--snakefile", str(WORKFLOW_DIR / "Snakefile"),
        "--directory", str(WORKFLOW_DIR),
        target,
        "--cores", "8",
        "--rerun-incomplete"
    ]
    print(f"Running Snakemake for manifest: {path}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("Snakemake FAILED")
        print(result.stderr)
        fatal_errors = [
            "no Snakefile found",
            "WorkflowError",
            "SyntaxError"
        ]
        if any(err in result.stderr for err in fatal_errors):
            print("Fatal error detected. Exiting monitor.")
            raise SystemExit(1)
        raise subprocess.CalledProcessError(
            result.returncode, cmd, result.stdout, result.stderr
        )
    print("Snakemake completed successfully")
    marker = DONE_DIR / f"{path.name}.done"
    marker.write_text("done\n")

def main():
    while True:
        for path in MANIFEST_DIR.glob("*.csv"):
            marker = DONE_DIR / f"{path.name}.done"
            if marker.exists():
                continue
            if not file_is_stable(path):
                continue
            try:
                run_snakemake(path)
            except subprocess.CalledProcessError as e:
                print(f"Snakemake failed for {path}: {e}")
        time.sleep(POLL_SECONDS)

if __name__ == "__main__":
    main()

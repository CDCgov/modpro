#!/usr/bin/env python3

from pathlib import Path
import subprocess
import time
import smtplib
from email.message import EmailMessage
import socket
import sys

WORKFLOW_DIR = Path("/scicomp/groups/modpro").resolve()
SNAKEFILE = WORKFLOW_DIR / "Snakefile"
MANIFEST_DIR = WORKFLOW_DIR / "manifests"
DONE_DIR = MANIFEST_DIR / ".processed"
DONE_DIR.mkdir(exist_ok=True)
<<<<<<< HEAD
WORKFLOW_DIR = Path("/scicomp/groups/modpro").resolve()

=======
EMAIL_TO = "qxa4@cdc.gov"
EMAIL_FROM = "qxa4@cdc.gov"
>>>>>>> 5c3193d (update monitor script)
POLL_SECONDS = 10
CORES = "8"

def send_failure_email(subject: str, body: str):
    msg = EmailMessage()
    msg["To"] = EMAIL_TO
    msg["From"] = EMAIL_FROM
    msg["Subject"] = subject
    msg.set_content(body)

    # Most HPC/Linux systems expose local sendmail via localhost.
    # If this fails, the error is printed before the monitor exits.
    with smtplib.SMTP("localhost") as smtp:
        smtp.send_message(msg)

def file_is_stable(path: Path, wait_seconds: int = 5) -> bool:
    size1 = path.stat().st_size
    time.sleep(wait_seconds)
    size2 = path.stat().st_size
    return size1 == size2

def run_snakemake(path: Path):
    mid = path.stem
    target = f"results/{mid}/final/ha1_selected.pdb"
    cmd = [
        "snakemake",
        "--snakefile", str(SNAKEFILE),
        "--directory", str(WORKFLOW_DIR),
        target,
        "--cores", CORES,
        "--rerun-incomplete",
    ]
    print(f"Running Snakemake for manifest: {path}", flush=True)
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        failed_marker = DONE_DIR / f"{path.name}.failed"
        failed_marker.write_text(
            f"FAILED\n\nCommand:\n{' '.join(cmd)}\n\n"
            f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}\n"
        )
        subject = f"[modpro] Snakemake failed for {path.name}"
        body = f"""
Snakemake failed and the manifest monitor is exiting.

Host:
{socket.gethostname()}

Workflow directory:
{WORKFLOW_DIR}

Manifest:
{path}

Target:
{target}

Exit code:
{result.returncode}

Command:
{' '.join(cmd)}

STDOUT:
{result.stdout}

STDERR:
{result.stderr}

A failed marker was written to:
{failed_marker}
"""
        print("Snakemake FAILED. Sending failure email and exiting.", flush=True)
        print(result.stdout, flush=True)
        print(result.stderr, flush=True)
        try:
            send_failure_email(subject, body)
            print(f"Failure email sent to {EMAIL_TO}", flush=True)
        except Exception as email_error:
            print(f"Failed to send email: {email_error}", flush=True)
        sys.exit(result.returncode)
    print("Snakemake completed successfully", flush=True)
    done_marker = DONE_DIR / f"{path.name}.done"
    done_marker.write_text("done\n")

def main():
    if not SNAKEFILE.exists():
        msg = f"Snakefile not found: {SNAKEFILE}"
        try:
            send_failure_email("[modpro] Manifest monitor failed to start", msg)
        except Exception as email_error:
            print(f"Failed to send startup failure email: {email_error}", flush=True)
        raise SystemExit(msg)
    while True:
        for path in sorted(MANIFEST_DIR.glob("*.csv")):
            done_marker = DONE_DIR / f"{path.name}.done"
            failed_marker = DONE_DIR / f"{path.name}.failed"
            if done_marker.exists() or failed_marker.exists():
                continue
            if not file_is_stable(path):
                continue
            run_snakemake(path)
        time.sleep(POLL_SECONDS)

if __name__ == "__main__":
    main()

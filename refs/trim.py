#!/usr/bin/env python3

from pathlib import Path
import shutil

#REF_DIR = Path("/scicomp/groups/modpro/refs")
REF_DIR = Path.cwd()
HA1_BOUNDS_BY_PREFIX = {
#    "Hvic": (16, 362),   # seasonal-style B HA
    "B": (16, 362),   # BRISBANE60-style B HA
    "H1": (18, 344),   # CALI07
    "H3": (17, 345),   # HK4801
}

def residue_number(line):
    try:
        return int(line[22:26].strip())
    except ValueError:
        return None

def chain_id(line):
    return line[21].strip()

def trim_pdb_in_place(path: Path):
    prefix = path.name[:2]
    if path.name.startswith("H1"):
        start, end = HA1_BOUNDS_BY_PREFIX["H1"]
    elif path.name.startswith("H3"):
        start, end = HA1_BOUNDS_BY_PREFIX["H3"]
    elif path.name.startswith("B_"):
        start, end = HA1_BOUNDS_BY_PREFIX["B"]
    else:
        raise ValueError(f"Unsupported reference prefix for {path.name}")
    backup = path.with_suffix(path.suffix + ".bak")
    shutil.copy2(path, backup)
    kept = 0
    with backup.open() as inp, path.open("w") as out:
        out.write(f"REMARK trimmed HA1 from original file: {backup.name}\n")
        out.write(f"REMARK kept chain A only\n")
        out.write(f"REMARK kept original numbering; residues {start}..{end}\n")
        for line in inp:
            record = line[0:6].strip()
            if record not in {"ATOM", "HETATM"}:
                continue
            if chain_id(line) != "A":
                continue
            resnum = residue_number(line)
            if resnum is None:
                continue
            if start <= resnum <= end:
                out.write(line)
                kept += 1
        out.write("TER\n")
        out.write("END\n")
    if kept == 0:
        shutil.copy2(backup, path)
        raise RuntimeError(f"No atoms kept for {path}; restored from backup")
    print(
        f"Trimmed {path.name}: "
        f"chain A only, residues {start}..{end}; "
        f"backup={backup.name}"
    )

def main():
    pdbs = sorted(REF_DIR.glob("*_HA_*.pdb"))
    if not pdbs:
        raise RuntimeError(f"No reference PDBs found in {REF_DIR}")
    for path in pdbs:
        trim_pdb_in_place(path)

if __name__ == "__main__":
    main()

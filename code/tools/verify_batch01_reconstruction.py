#!/usr/bin/env python3
"""Restore and verify original HAR Batch01 inputs from the GitHub repository.

All 19 originals can be reconstructed byte-for-byte using the canonical paths
in provenance/batches/B01_2026-11-21/BATCH01_FILE_LEDGER.csv.

Usage:
    python code/tools/verify_batch01_reconstruction.py
    python code/tools/verify_batch01_reconstruction.py --output /tmp/har-batch01

Only Python standard library; no credentials or network are needed once the repo
has been cloned. External GEO and published chromatin files retain source identity.
"""
import argparse
import csv
import hashlib
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LEDGER = ROOT / "provenance/batches/B01_2026-11-21/BATCH01_FILE_LEDGER.csv"


def restore_one(row, output):
    relative_path = row["canonical_github_path"]
    count = int(row["part_count"])
    location = ROOT / relative_path
    destination = output / row["original_filename"]
    destination.parent.mkdir(parents=True, exist_ok=True)
    files = [location] if count == 1 else [
        location / f"part-{i:05d}" for i in range(1, count + 1)
    ]
    if any(not f.is_file() for f in files):
        missing = [str(f) for f in files if not f.is_file()]
        raise FileNotFoundError(f"{row['original_filename']}: missing {missing}")
    sha = hashlib.sha256()
    size = 0
    with destination.open("wb") as dst:
        for f in files:
            with f.open("rb") as src:
                while True:
                    chunk = src.read(1024 * 1024)
                    if not chunk:
                        break
                    dst.write(chunk)
                    sha.update(chunk)
                    size += len(chunk)
    if size != int(row["bytes"]) or sha.hexdigest() != row["sha256"]:
        destination.unlink(missing_ok=True)
        raise ValueError(f"Integrity FAILED: {row['original_filename']}")
    return destination, size


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "reconstructed_batch01")
    args = parser.parse_args()
    if not LEDGER.is_file():
        parser.error("Repository ledger not found. Run inside the full GitHub repository.")
    with LEDGER.open(encoding="utf-8", newline="") as fh:
        records = list(csv.DictReader(fh))
    if len(records) != 19:
        raise ValueError(f"Expected 19 original Batch01 files, got {len(records)}")
    restored = []
    for row in records:
        target, size = restore_one(row, args.output)
        restored.append(target)
        print(f"PASS {row['index']:>2}/19  {size:>9} bytes  {target.name}")
    print(f"Verified {len(restored)}/19; all SHA-256 checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Run FragGeneScanRs with a compatibility fallback across CLI styles.

Inputs:
  1) FASTA with reads/fragments
  2) output FAA path
  3) model (default: complete)

Usage:
  python scripts/run_fraggenescanrs.py reads.fa genes.faa [model]
"""

import shutil
import subprocess
import sys
from pathlib import Path


def run_command(cmd, stdin_path=None, stdout_path=None):
    stdin_handle = open(stdin_path, "rb") if stdin_path else None
    stdout_handle = open(stdout_path, "wb") if stdout_path else None
    try:
        proc = subprocess.run(cmd, stdin=stdin_handle, stdout=stdout_handle, stderr=subprocess.PIPE)
        return proc.returncode, proc.stderr.decode(errors="replace")
    finally:
        if stdin_handle is not None:
            stdin_handle.close()
        if stdout_handle is not None:
            stdout_handle.close()


def main():
    if len(sys.argv) < 3:
        print("Usage: run_fraggenescanrs.py <reads.fa> <genes.faa> [model]", file=sys.stderr)
        return 2

    reads_fa = Path(sys.argv[1])
    out_faa = Path(sys.argv[2])
    model = sys.argv[3] if len(sys.argv) > 3 else "complete"

    if not reads_fa.exists():
        print(f"Error: input FASTA not found: {reads_fa}", file=sys.stderr)
        return 1

    exe = shutil.which("FragGeneScanRs") or shutil.which("fraggenescanrs")
    if exe is None:
        print("Error: FragGeneScanRs executable not found in PATH", file=sys.stderr)
        return 1

    out_faa.parent.mkdir(parents=True, exist_ok=True)

    # Style 1 (documented backward-compatible mode): stdin -> stdout
    code, err = run_command([exe, "-t", model], stdin_path=str(reads_fa), stdout_path=str(out_faa))
    if code == 0 and out_faa.exists() and out_faa.stat().st_size > 0:
        return 0

    # Style 2 (subcommand style in some builds): predict
    tmp_err = err
    code, err = run_command([exe, "predict", "-t", model, "-s", str(reads_fa), "-o", str(out_faa)])
    if code == 0 and out_faa.exists() and out_faa.stat().st_size > 0:
        return 0

    print("Error: FragGeneScanRs failed in both CLI modes.", file=sys.stderr)
    print("Mode1 stderr:", file=sys.stderr)
    print(tmp_err.strip(), file=sys.stderr)
    print("Mode2 stderr:", file=sys.stderr)
    print(err.strip(), file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

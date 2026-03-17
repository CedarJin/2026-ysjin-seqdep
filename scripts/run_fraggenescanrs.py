#!/usr/bin/env python3
"""
Run FragGeneScanRs with compatibility fallbacks across CLI styles.

Usage:
  python scripts/run_fraggenescanrs.py reads.fa genes.faa [model] [--meta-out metadata.tsv] [--threads N]
"""

import argparse
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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run FragGeneScanRs and write predicted proteins (and optional metadata)."
    )
    parser.add_argument("reads_fa", help="Input reads/fragments in FASTA format")
    parser.add_argument("genes_faa", help="Output predicted proteins (FAA)")
    parser.add_argument(
        "model",
        nargs="?",
        default="complete",
        help="FragGeneScanRs training model (default: complete)",
    )
    parser.add_argument(
        "--meta-out",
        default=None,
        help="Optional output path for FragGeneScanRs metadata (-m)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads for FragGeneScanRs (-p)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    reads_fa = Path(args.reads_fa)
    out_faa = Path(args.genes_faa)
    model = args.model
    meta_out = Path(args.meta_out) if args.meta_out else None
    threads = max(1, int(args.threads))

    if not reads_fa.exists():
        print(f"Error: input FASTA not found: {reads_fa}", file=sys.stderr)
        return 1

    exe = shutil.which("FragGeneScanRs") or shutil.which("fraggenescanrs")
    if exe is None:
        print("Error: FragGeneScanRs executable not found in PATH", file=sys.stderr)
        return 1

    out_faa.parent.mkdir(parents=True, exist_ok=True)
    if meta_out is not None:
        meta_out.parent.mkdir(parents=True, exist_ok=True)

    # When metadata is requested, use explicit file-output style to guarantee .meta creation.
    if meta_out is not None:
        code, err = run_command(
            [
                exe,
                "-s",
                str(reads_fa),
                "-t",
                model,
                "-p",
                str(threads),
                "-a",
                str(out_faa),
                "-m",
                str(meta_out),
            ]
        )
        if (
            code == 0
            and out_faa.exists()
            and out_faa.stat().st_size > 0
            and meta_out.exists()
            and meta_out.stat().st_size > 0
        ):
            return 0
        print("Error: FragGeneScanRs failed while writing metadata.", file=sys.stderr)
        print(err.strip(), file=sys.stderr)
        return 1

    # Style 1 (legacy compatibility mode): stdin -> stdout (proteins only)
    code, err = run_command([exe, "-t", model], stdin_path=str(reads_fa), stdout_path=str(out_faa))
    if code == 0 and out_faa.exists() and out_faa.stat().st_size > 0:
        return 0

    # Style 2 (explicit file-output mode)
    tmp_err = err
    code, err = run_command(
        [exe, "-s", str(reads_fa), "-t", model, "-p", str(threads), "-a", str(out_faa)]
    )
    if code == 0 and out_faa.exists() and out_faa.stat().st_size > 0:
        return 0

    # Style 3 (subcommand mode in some builds)
    style2_err = err
    code, err = run_command([exe, "predict", "-t", model, "-s", str(reads_fa), "-o", str(out_faa)])
    if code == 0 and out_faa.exists() and out_faa.stat().st_size > 0:
        return 0

    print("Error: FragGeneScanRs failed in all tested CLI modes.", file=sys.stderr)
    print("Mode1 stderr:", file=sys.stderr)
    print(tmp_err.strip(), file=sys.stderr)
    print("Mode2 stderr:", file=sys.stderr)
    print(style2_err.strip(), file=sys.stderr)
    print("Mode3 stderr:", file=sys.stderr)
    print(err.strip(), file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

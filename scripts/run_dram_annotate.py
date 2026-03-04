#!/usr/bin/env python3
"""
Run DRAM.py annotate, streaming its log. When DRAM prints messages about
missing dbCAN/glycan subfamily descriptions, we:

- let DRAM/DatabaseHandler handle subfamily→family fallback internally, and
- emit a short NOTE in the log, but do NOT fail the run.

Usage:
  python scripts/run_dram_annotate.py [DRAM arguments...]
Example:
  python scripts/run_dram_annotate.py -i 'bins/*.fa' -o out --threads 16 ...
"""

import re
import subprocess
import sys

# Pattern that indicates glycan/CAZy subfamily description lookup failed
DBCAN_DESC_MISSING = re.compile(
    r"No descriptions were found.*dbcan_description",
    re.IGNORECASE,
)


def main():
    argv = sys.argv[1:]
    if not argv:
        print("Usage: run_dram_annotate.py <DRAM.py annotate arguments...>", file=sys.stderr)
        sys.exit(1)

    cmd = ["DRAM.py", "annotate"] + argv
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    for line in proc.stdout:
        print(line, end="")
        if DBCAN_DESC_MISSING.search(line):
            # This indicates that DRAM could not find a dbCAN description for
            # some IDs. DatabaseHandler now first tries a subfamily→family
            # fallback; if that also fails we simply continue without a
            # description. Here we only surface an informational NOTE.
            sys.stdout.flush()
            print(
                "\n[run_dram_annotate] NOTE: DRAM reported missing dbCAN/glycan subfamily "
                "description for some IDs. Subfamily→family fallback (if available) will "
                "be used; otherwise annotation continues without a description.",
                file=sys.stderr,
            )
    ret = proc.wait()
    sys.exit(ret)


if __name__ == "__main__":
    main()

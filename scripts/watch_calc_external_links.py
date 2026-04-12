#!/usr/bin/env python3
"""Watch comprehensive test CSVs and rebuild Calc link workbooks on change."""

from __future__ import annotations

import argparse
import subprocess
import time
from pathlib import Path


def snapshot(paths: list[Path]) -> dict[str, tuple[int, int]]:
    state: dict[str, tuple[int, int]] = {}
    for path in paths:
        stat = path.stat()
        state[str(path)] = (stat.st_mtime_ns, stat.st_size)
    return state


def run_refresh(script: Path, csv_glob: str) -> int:
    proc = subprocess.run([str(script), csv_glob], check=False)
    return proc.returncode


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "csv_glob",
        nargs="?",
        default="package_tests/comprehensive_tests_results_nc_1_*.csv",
    )
    parser.add_argument("--interval", type=float, default=5.0)
    parser.add_argument("--skip-initial-refresh", action="store_true")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    refresh_script = root / "scripts" / "refresh_calc_external_links.sh"

    csv_paths = sorted(root.glob(args.csv_glob))
    if not csv_paths:
        raise SystemExit(f"No CSV files matched: {args.csv_glob}")

    print(f"Watching {len(csv_paths)} CSV files.")
    print(f"Polling interval: {args.interval:.1f}s")
    print(f"Refresh command: {refresh_script} {args.csv_glob}")

    current = snapshot(csv_paths)
    if not args.skip_initial_refresh:
        print("Running initial refresh...")
        code = run_refresh(refresh_script, args.csv_glob)
        if code != 0:
            print(f"Initial refresh failed with exit code {code}.")

    while True:
        time.sleep(args.interval)
        latest_paths = sorted(root.glob(args.csv_glob))
        if not latest_paths:
            continue
        latest = snapshot(latest_paths)
        if latest != current:
            print("Change detected. Refreshing Calc link workbooks...")
            code = run_refresh(refresh_script, args.csv_glob)
            if code == 0:
                current = latest
                print("Refresh complete.")
            else:
                print(f"Refresh failed with exit code {code}.")


if __name__ == "__main__":
    raise SystemExit(main())

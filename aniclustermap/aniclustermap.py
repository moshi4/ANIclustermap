#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import platform
import subprocess as sp
from pathlib import Path

__version__ = "0.1.0"


def main():
    """ANIclustermap main function for entrypoint"""
    # Get argument values
    args = get_args()

    indir: Path = args.indir
    outdir: Path = args.outdir
    thread_num: int = args.thread_num

    run(indir, outdir, thread_num)


def run(indir: Path, outdir: Path, thread_num: int = 1) -> None:
    """Run ANIclustermap workflow"""
    add_bin_path()


def add_bin_path() -> None:
    """Add executable binary (fastANI) path to PATH"""
    os_name = platform.system()  # 'Windows' or 'Darwin' or 'Linux'
    bin_path = Path(__file__).parent / "bin" / os_name
    sep = ";" if os_name == "Windows" else ":"
    env_path = f"{bin_path}{sep}{os.environ['PATH']}"
    os.environ["PATH"] = env_path


def get_args() -> argparse.Namespace:
    """Get arguments

    Returns:
        argparse.Namespace: Argument values
    """
    description = "Draw ANI(Average Nucleotide Identity) clustermap"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-i",
        "--indir",
        required=True,
        type=Path,
        help="Input genome sequences directory",
        metavar="",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        type=Path,
        help="Output directory",
        metavar="",
    )
    cpu_num = os.cpu_count()
    default_thread_num = 1 if cpu_num is None or cpu_num == 1 else cpu_num - 1
    parser.add_argument(
        "-t",
        "--thread_num",
        type=int,
        help=f"fastANI thread numter parameter (Default: {default_thread_num})",
        default=default_thread_num,
        metavar="",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"v{__version__}",
        help="Print version information",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()

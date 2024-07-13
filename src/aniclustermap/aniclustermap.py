from __future__ import annotations

import argparse
import csv
import os
import platform
import re
import shlex
import subprocess as sp
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hc
import seaborn as sns
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, is_color_like
from scipy.cluster.hierarchy import ClusterNode
from seaborn.matrix import ClusterGrid

from aniclustermap.logger import get_logger

__version__ = "1.4.0"

logger = get_logger("aniclustermap")

RUN_MODE = Literal["fastani", "skani"]


def main():
    """ANIclustermap main function for entrypoint"""
    # Get argument values
    args = get_args()
    run(**args.__dict__)


def run(
    indir: Path,
    outdir: Path,
    mode: RUN_MODE = "fastani",
    thread_num: int = 1,
    overwrite: bool = False,
    fig_width: int = 10,
    fig_height: int = 10,
    dendrogram_ratio: float = 0.15,
    cmap_colors: list[str] | None = None,
    cmap_gamma: float = 1.0,
    cmap_ranges: list[float] | None = None,
    cbar_pos: tuple[float, float, float, float] = (0.02, 0.8, 0.05, 0.18),
    annotation: bool = False,
    annotation_fmt: str = ".3g",
    skani_c_param: int = 30,
) -> None:
    """Run ANIclustermap workflow"""
    outdir.mkdir(exist_ok=True)
    workdir = outdir / "work"
    workdir.mkdir(exist_ok=True)

    ani_result_file = workdir / f"{mode}_result"
    ani_matrix_file = Path(str(ani_result_file) + ".matrix")
    ani_matrix_tsv_file = workdir / f"{mode}_matrix.tsv"
    if not ani_matrix_tsv_file.exists() or overwrite:
        # Run ANI calculation
        fasta_list_file = workdir / "genome_fasta_file_list.txt"
        genome_num = write_genome_fasta_list(indir, fasta_list_file)
        if genome_num <= 1:
            logger.error("ERROR: Number of input genome fasta file is less than 1.")
            exit(1)

        logger.info(f"# Step1: Run {mode} between all-vs-all {genome_num} genomes.")
        add_bin_path()
        if mode == "fastani":
            run_fastani(fasta_list_file, ani_result_file, thread_num)
        else:
            run_skani(fasta_list_file, ani_matrix_file, thread_num, skani_c_param)
        ani_df = parse_ani_matrix(ani_matrix_file)
        ani_df.to_csv(ani_matrix_tsv_file, sep="\t", index=False)
    else:
        logger.info(f"# Step1: Previous {mode} matrix result found. Skip {mode} run.")
        ani_df = pd.read_csv(ani_matrix_tsv_file, sep="\t", index_col=False)
        ani_df = ani_df.set_index(ani_df.columns)

    # Hierarchical clustering ANI matrix
    logger.info(f"# Step2: Clustering {mode} matrix by scipy's UPGMA method.")
    linkage = hc.linkage(ani_df, method="average")

    # Output dendrogram tree as newick format tree
    tree = hc.to_tree(linkage)
    if isinstance(tree, ClusterNode):
        dendrogram_newick_file = outdir / "ANIclustermap_dendrogram.nwk"
        with open(dendrogram_newick_file, "w") as f:
            leaf_names = list(map(str, ani_df.columns))
            f.write(dendrogram2newick(tree, tree.dist, leaf_names))
    else:
        raise ValueError("Invalid hierarchy cluster detected!!")

    # Draw ANI clustermap
    logger.info("# Step3: Using clustered matrix, draw ANI clustermap by seaborn.\n")
    cmap_colors = ["lime", "yellow", "red"] if cmap_colors is None else cmap_colors
    if cmap_ranges is None:
        mycmap = LinearSegmentedColormap.from_list(
            "mycmap", colors=cmap_colors, gamma=cmap_gamma
        )
        opts = {}
    else:
        mycmap = LinearSegmentedColormap.from_list(
            "mycmap", colors=cmap_colors, gamma=cmap_gamma, N=len(cmap_ranges) - 1
        )
        opts = {"norm": BoundaryNorm(cmap_ranges, len(cmap_ranges) - 1)}
    mycmap.set_under("lightgrey")

    min_ani = min(filter(lambda v: v != 0, np.array(ani_df).flatten()))
    g: ClusterGrid = sns.clustermap(
        data=ani_df,
        col_linkage=linkage,
        row_linkage=linkage,
        figsize=(fig_width, fig_height),
        annot=annotation,
        fmt=annotation_fmt,
        cmap=mycmap,
        dendrogram_ratio=dendrogram_ratio,
        xticklabels=False,
        yticklabels=True,
        vmin=min_ani,
        vmax=100,
        cbar=True,
        cbar_pos=cbar_pos,
        cbar_kws={
            "label": "ANI (%)",
            "orientation": "vertical",
            "spacing": "proportional",
            # "extend": "min",
            # "extendfrac": 0.1,
        },
        tree_kws={"linewidths": 1.5},
        **opts,  # type: ignore
    )
    # Get clustered ani matrix dataframe
    clustered_ani_df = get_clustered_matrix(ani_df, g)
    clustered_ani_matrix_tsv_file = outdir / "ANIclustermap_matrix.tsv"
    clustered_ani_df.to_csv(clustered_ani_matrix_tsv_file, sep="\t", index=False)

    # Output ANI clustermap figure
    ani_clustermap_png_file = outdir / "ANIclustermap.png"
    ani_clustermap_svg_file = outdir / "ANIclustermap.svg"
    g.savefig(ani_clustermap_png_file)
    g.savefig(ani_clustermap_svg_file)


def write_genome_fasta_list(
    target_dir: Path,
    list_outfile: Path,
    exts: list[str] | None = None,
) -> int:
    """Write genome fasta file list for ANI comparison

    Parameters
    ----------
    target_dir : Path
        Target genome fasta directory
    list_outfile : Path
        List of genome fasta file path
    exts : list[str] | None, optional
        Genome fasta target extension list
        (Default: `.fa`,`.fna`,`.fna.gz`,`.fasta`)

    Returns
    -------
    count : int
        Number of file
    """
    # Get target file path list
    file_path_list = []
    exts = [".fa", ".fna", ".fna.gz", ".fasta"] if exts is None else exts
    for ext in exts:
        file_path_list.extend(target_dir.glob(f"*{ext}"))

    # Write file path list
    contents = "\n".join([str(f) for f in file_path_list])
    with open(list_outfile, "w") as f:
        f.write(contents)

    return len(file_path_list)


def run_fastani(
    genome_fasta_list_file: Path,
    fastani_result_file: Path,
    thread_num: int,
) -> None:
    """Run fastANI

    Parameters
    ----------
    genome_fasta_list_file : Path
        Genome fasta file list
    fastani_result_file : Path
        fastANI result output file
    thread_num : int
        Thread number for fastANI run
    """
    cmd = f"fastANI --ql {genome_fasta_list_file} --rl {genome_fasta_list_file} -o {fastani_result_file} -t {thread_num} --matrix"  # noqa: E501
    run_cmd(cmd)
    if fastani_result_file.exists():
        fastani_result_file.unlink()


def run_skani(
    genome_fasta_list_file: Path,
    result_file: str | Path,
    thread_num: int,
    c_param: int = 30,
) -> None:
    """Run skani

    Parameters
    ----------
    genome_fasta_list_file : Path
        Genome fasta file list
    result_file : str | Path
        skani result output file
    thread_num : int
        Thread number for skani run
    c_param : int, optional
        Compression factor parameter
    """
    cmd = f"skani triangle -l {genome_fasta_list_file} -o {result_file} -t {thread_num} -c {c_param}"  # noqa: E501
    run_cmd(cmd)


def run_cmd(cmd: str) -> None:
    """Run command"""
    cmd_res = sp.run(shlex.split(cmd), capture_output=True, text=True)
    if cmd_res.returncode != 0:
        logger.error("Failed to run command below!!")
        logger.error(f"$ {cmd}")
        stdout_lines = cmd_res.stdout.splitlines()
        if len(stdout_lines) > 0:
            logger.error("STDOUT:")
            for line in stdout_lines:
                logger.error(f"> {line}")
        stderr_lines = cmd_res.stderr.splitlines()
        if len(stderr_lines) > 0:
            logger.error("STDERR:")
            for line in stderr_lines:
                logger.error(f"> {line}")
        exit(1)


def add_bin_path() -> None:
    """Add executable binary (fastANI) path to PATH"""
    os_name = platform.system()  # 'Windows' or 'Darwin' or 'Linux'
    bin_path = Path(__file__).parent / "bin" / os_name
    sep = ";" if os_name == "Windows" else ":"
    env_path = f"{os.environ['PATH']}{sep}{bin_path}"
    os.environ["PATH"] = env_path


def parse_ani_matrix(matrix_file: Path) -> pd.DataFrame:
    """Parse ANI matrix as Dataframe

    Parameters
    ----------
    matrix_file : Path
        All-vs-All ANI matrix file

    Returns
    -------
    df : pd.DataFrame
        Dataframe of ANI matrix
    """
    names: list[str] = []
    ani_values_list: list[list[float]] = []
    with open(matrix_file) as f:
        reader = csv.reader(f, delimiter="\t")
        genome_num = int(next(reader)[0].rstrip("\n"))
        for row in reader:
            name = Path(row[0]).with_suffix("").name
            name = re.sub("\\.fna$", "", name)
            names.append(name)
            ani_values = list(map(lambda d: 0.0 if d == "NA" else float(d), row[1:]))
            ani_values.extend([0] * (genome_num - len(ani_values)))
            ani_values_list.append(ani_values)

    df = pd.DataFrame(data=ani_values_list, columns=names, index=names)
    for i in range(genome_num):
        df.iat[i, i] = 100
    for i, name in enumerate(names):
        for j, d in enumerate(df[name][i:]):
            df.iat[i, i + j] = d
    return df


def dendrogram2newick(
    node: ClusterNode, parent_dist: float, leaf_names: list[str], newick: str = ""
) -> str:
    """Convert scipy dendrogram tree to newick format tree

    Parameters
    ----------
    node : ClusterNode
        Tree node
    parent_dist : float
        Parent distance
    leaf_names : list[str]
        Leaf names
    newick : str, optional
        Newick format string (Used in recursion)

    Returns
    -------
    newick : str
        Newick format tree
    """
    if node.is_leaf():
        return f"{leaf_names[node.id]}:{(parent_dist - node.dist):.2f}{newick}"
    else:
        if len(newick) > 0:
            newick = f"):{(parent_dist - node.dist):.2f}{newick}"
        else:
            newick = ");"
        if node.left is None or node.right is None:
            raise ValueError
        newick = dendrogram2newick(node.left, node.dist, leaf_names, newick)
        newick = dendrogram2newick(node.right, node.dist, leaf_names, f",{newick}")
        newick = f"({newick}"
        return newick


def get_clustered_matrix(original_df: pd.DataFrame, g: ClusterGrid) -> pd.DataFrame:
    """Get clustered ANI matrix

    Parameters
    ----------
    original_df : pd.DataFrame
        Original dataframe before clustering
    g : ClusterGrid
        Cluster grid (`clustermap` return value)

    Returns
    -------
    df : pd.DataFrame
        Clustered matrix dataframe
    """
    clustered_row_index = original_df.index[g.dendrogram_row.reordered_ind]
    clustered_col_index = original_df.columns[g.dendrogram_col.reordered_ind]
    return original_df.loc[clustered_row_index, clustered_col_index]  # type: ignore


def get_args() -> argparse.Namespace:
    """Get arguments

    Returns
    -------
    args : argparse.Namespace
        Argument values
    """
    description = "Draw ANI(Average Nucleotide Identity) clustermap"
    parser = argparse.ArgumentParser(
        usage="ANIclustermap -i [Genome fasta directory] -o [output directory]",
        description=description,
        add_help=False,
    )

    parser.add_argument(
        "-i",
        "--indir",
        required=True,
        type=Path,
        help="Input genome fasta directory (*.fa|*.fna[.gz]|*.fasta)",
        metavar="I",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        type=Path,
        help="Output directory",
        metavar="O",
    )
    default_mode = "fastani"
    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        help="ANI calculation mode ('fastani'[default]|'skani')",
        default=default_mode,
        choices=["fastani", "skani"],
        metavar="",
    )
    cpu_num = os.cpu_count()
    default_thread_num = 1 if cpu_num is None or cpu_num == 1 else cpu_num - 1
    parser.add_argument(
        "-t",
        "--thread_num",
        type=int,
        help=f"Thread number parameter (Default: {default_thread_num})",
        default=default_thread_num,
        metavar="",
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite previous ANI calculation result (Default: OFF)",
        action="store_true",
    )
    default_fig_width = 10
    parser.add_argument(
        "--fig_width",
        type=int,
        help=f"Figure width (Default: {default_fig_width})",
        default=default_fig_width,
        metavar="",
    )
    default_fig_height = 10
    parser.add_argument(
        "--fig_height",
        type=int,
        help=f"Figure height (Default: {default_fig_height})",
        default=default_fig_height,
        metavar="",
    )
    default_dendrogram_ratio = 0.15
    parser.add_argument(
        "--dendrogram_ratio",
        type=float,
        help=f"Dendrogram ratio to figsize (Default: {default_dendrogram_ratio})",
        default=default_dendrogram_ratio,
        metavar="",
    )
    default_cmap_colors = "lime,yellow,red"
    parser.add_argument(
        "--cmap_colors",
        type=str,
        help=f"cmap interpolation colors parameter (Default: '{default_cmap_colors}')",
        default=default_cmap_colors,
        metavar="",
    )
    default_cmap_gamma = 1.0
    parser.add_argument(
        "--cmap_gamma",
        type=float,
        help=f"cmap gamma parameter (Default: {default_cmap_gamma})",
        default=default_cmap_gamma,
        metavar="",
    )
    parser.add_argument(
        "--cmap_ranges",
        type=str,
        help="Range values (e.g. 80,90,95,100) for discrete cmap (Default: None)",
        default=None,
        metavar="",
    )
    default_cbar_pos = (0.02, 0.8, 0.05, 0.18)
    parser.add_argument(
        "--cbar_pos",
        type=float,
        nargs=4,
        help=f"Colorbar position (Default: {default_cbar_pos})",
        default=default_cbar_pos,
        metavar="",
    )
    parser.add_argument(
        "--annotation",
        help="Show ANI value annotation (Default: OFF)",
        action="store_true",
    )
    default_annotation_fmt = ".3g"
    parser.add_argument(
        "--annotation_fmt",
        type=str,
        help=f"Annotation value format (Default: '{default_annotation_fmt}')",
        default=default_annotation_fmt,
        metavar="",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"v{__version__}",
        help="Print version information",
    )
    parser.add_argument(
        "-h",
        "--help",
        help="Show this help message and exit",
        action="help",
    )

    default_skani_c_param = 30
    parser.add_argument(
        "--skani_c_param",
        type=int,
        help=argparse.SUPPRESS,
        default=default_skani_c_param,
    )

    args = parser.parse_args()

    # Validate cmap color string
    cmap_colors = args.cmap_colors.split(",")
    for cmap_color in cmap_colors:
        if not is_color_like(cmap_color):
            parser.error(
                f"--cmap_colors: '{cmap_color}' is not valid color like string!!"
            )
    if len(cmap_colors) <= 1:
        parser.error("--cmap_colors: Multiple colors are expected.")

    # Validate range values for discrete cmap
    if args.cmap_ranges is not None:
        try:
            cmap_ranges = [float(v) for v in args.cmap_ranges.split(",")]
        except ValueError:
            parser.error("--cmap_ranges: Contains Non-float values.")
        if len(cmap_ranges) <= 1:
            parser.error("--cmap_ranges: Multiple range values are expected.")
        if cmap_ranges != sorted(cmap_ranges):
            parser.error("--cmap_ranges: Specify range values in ascending order.")
        for cmap_range in cmap_ranges:
            if not 70 <= cmap_range <= 100:
                parser.error("--cmap_ranges: Range values must be 70 <= value <= 100.")
        if max(cmap_ranges) != 100:
            parser.error("--cmap_ranges: Max range value must be 100.")

    args.cmap_colors = args.cmap_colors.split(",")
    if args.cmap_ranges is not None:
        args.cmap_ranges = [float(v) for v in args.cmap_ranges.split(",")]

    return args


if __name__ == "__main__":
    main()

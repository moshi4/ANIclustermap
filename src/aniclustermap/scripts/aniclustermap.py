import logging
import platform
import sys
from enum import Enum
from functools import partial
from pathlib import Path
from typing import Annotated, Optional

import pandas as pd
import scipy.cluster.hierarchy as hc
import typer
from scipy.spatial.distance import squareform
from typer import Option, Typer

from aniclustermap import __version__, const
from aniclustermap.ani.matrix import get_clustered_matrix
from aniclustermap.ani.tools import FastAni, SkAni
from aniclustermap.cluster import clustermap
from aniclustermap.logger import init_logger
from aniclustermap.tree import to_newick_tree
from aniclustermap.utils import exit_handler, logging_timeit

Option = partial(Option, metavar="")

app = Typer(add_completion=False)


def version_callback(v: bool):
    """Callback function for print version"""
    if v:
        print(f"v{__version__}")
        raise typer.Exit()


class AniCalcMode(str, Enum):
    fastani = "fastani"
    skani = "skani"


@app.command(
    no_args_is_help=True,
    epilog=None,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@logging_timeit
@exit_handler
def cli(
    indir: Annotated[
        Path,
        Option(
            "-i",
            "--indir",
            help="Input genome fasta directory (*.fa|*.fna[.gz]|*.fasta)",
            show_default=False,
            exists=True,
        ),
    ],
    outdir: Annotated[
        Path,
        Option("-o", "--outdir", help="Output directory", show_default=False),
    ],
    mode: Annotated[
        AniCalcMode,
        Option("--mode", help="ANI calculation tool (fastani|skani)"),
    ] = AniCalcMode.fastani,
    thread_num: Annotated[
        int,
        Option("-t", "--thread_num", help="Thread number parameter"),
    ] = const.DEFAULT_CPU,
    overwrite: Annotated[
        bool,
        Option("--overwrite", help="Overwrite previous ANI calculation result"),
    ] = False,
    fig_width: Annotated[
        int,
        Option("--fig_width", help="Figure width"),
    ] = 10,
    fig_height: Annotated[
        int,
        Option("--fig_height", help="Figure height"),
    ] = 10,
    dendrogram_ratio: Annotated[
        float,
        Option("--dendrogram_ratio", help="Dendrogram ratio to figsize"),
    ] = 0.15,
    cmap_colors: Annotated[  # type: ignore
        str,
        Option(
            "--cmap_colors",
            help="cmap interpolation colors parameter",
        ),
    ] = "lime,yellow,red",
    cmap_gamma: Annotated[
        float,
        Option("--cmap_gamma", help="cmap gamma parameter"),
    ] = 1.0,
    cmap_ranges: Annotated[  # type: ignore
        Optional[str],
        Option(
            "--cmap_ranges",
            help="Range values (e.g. 80,90,95,100) for discrete cmap",
            show_default=False,
        ),
    ] = None,
    cbar_pos: Annotated[
        tuple[float, float, float, float],
        Option("--cbar_pos", help="Colorbar position"),
    ] = (0.02, 0.85, 0.04, 0.15),
    annotation: Annotated[
        bool,
        Option("--annotation", help="Show ANI value annotation"),
    ] = False,
    annotation_fmt: Annotated[
        str,
        Option("--annotation_fmt", help="Annotation value format"),
    ] = ".3g",
    quiet: Annotated[
        bool,
        Option("--quiet", help="No print log on screen"),
    ] = False,
    debug: Annotated[
        bool,
        Option("--debug", help="Print debug log", hidden=True),
    ] = False,
    _: Annotated[
        bool,
        Option(
            "-v",
            "--version",
            help="Print version information",
            callback=version_callback,
            is_eager=True,
        ),
    ] = False,
) -> None:
    """Draw ANI(Average Nucleotide Identity) clustermap"""
    args = locals()
    outdir.mkdir(exist_ok=True)

    # Setup logger
    log_file = outdir / "aniclustermap.log"
    init_logger(quiet=quiet, verbose=debug, log_file=log_file)
    logger = logging.getLogger(__name__)

    logger.info(f"Run ANIclustermap v{__version__}")
    logger.info(f"$ {Path(sys.argv[0]).name} {' '.join(sys.argv[1:])}")
    logger.info(f"Operating System: {sys.platform}")
    logger.info(f"Python Version: v{platform.python_version()}")
    for name, value in args.items():
        if name not in ("quiet", "debug", "_"):
            logger.info(f"Parameter: {name}={value}")

    # Calculate ANI matrix by fastANI or skani
    tmpdir = None
    if debug:
        tmpdir = outdir / "tmp"
        tmpdir.mkdir(exist_ok=True)
    ani_matrix_tsv_file = outdir / f"{mode.value}_matrix.tsv"
    if not ani_matrix_tsv_file.exists() or overwrite:
        logger.info(f"Run {mode.value} between all-vs-all genomes")
        AniCalcTool = dict(fastani=FastAni, skani=SkAni)[mode.value]
        ani_matrix_df = AniCalcTool(indir, tmpdir).run(thread_num=thread_num)
        ani_matrix_df.to_csv(ani_matrix_tsv_file, sep="\t", index=False)
        logger.info("Write all-vs-all genomes ANI matrix result")
        logger.info(f"=> {ani_matrix_tsv_file}")
    else:
        logger.info(
            f"Previous {mode.value} matrix result found (={ani_matrix_tsv_file})"
        )
        logger.info(f"Skip {mode.value} run")
        ani_matrix_df = pd.read_csv(ani_matrix_tsv_file, sep="\t", index_col=False)
        ani_matrix_df = ani_matrix_df.set_index(ani_matrix_df.columns)

    # Hierarchical clustering ANI matrix
    logger.info(f"Clustering {mode.value} ANI matrix by scipy UPGMA method")
    linkage = hc.linkage(squareform(100 - ani_matrix_df), method="average")
    dendrogram_newick_file = outdir / "ANIclustermap_dendrogram.nwk"
    to_newick_tree(ani_matrix_df, linkage, dendrogram_newick_file)
    logger.info("Write newick format cluster dendrogram")
    logger.info(f"=> {dendrogram_newick_file}")

    # Draw ANI clustermap
    g = clustermap(
        ani_matrix_df,
        linkage=linkage,
        cmap_colors=cmap_colors,
        cmap_ranges=cmap_ranges,
        cmap_gamma=cmap_gamma,
        fig_width=fig_width,
        fig_height=fig_height,
        dendrogram_ratio=dendrogram_ratio,
        cbar_pos=cbar_pos,
        annotation=annotation,
        annotation_fmt=annotation_fmt,
    )

    # Get clustered ani matrix dataframe
    clustered_ani_df = get_clustered_matrix(ani_matrix_df, g)
    clustered_ani_matrix_tsv_file = outdir / "ANIclustermap_matrix.tsv"
    clustered_ani_df.to_csv(clustered_ani_matrix_tsv_file, sep="\t", index=False)
    logger.info("Write clustered ANI matrix")
    logger.info(f"=> {clustered_ani_matrix_tsv_file}")

    # Output ANI clustermap figure
    logger.info("Using clustered ANI matrix, draw ANI clustermap by seaborn")
    aniclustermap_png_file = outdir / "ANIclustermap.png"
    aniclustermap_svg_file = aniclustermap_png_file.with_suffix(".svg")
    g.savefig(aniclustermap_png_file)
    logger.info(f"=> {aniclustermap_png_file}")
    g.savefig(aniclustermap_svg_file)
    logger.info(f"=> {aniclustermap_svg_file}")


if __name__ == "__main__":
    app()

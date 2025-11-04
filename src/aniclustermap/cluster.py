from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hc
import seaborn as sns
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import LinearSegmentedColormap as LSC
from scipy.spatial.distance import squareform
from seaborn.matrix import ClusterGrid


def clustermap(
    ani_matrix_df: str | Path | pd.DataFrame,
    *,
    linkage: np.ndarray | None = None,
    cmap_colors: str | list[str] | None = None,  # type: ignore
    cmap_ranges: str | list[float] | None = None,  # type: ignore
    cmap_gamma: float = 1.0,
    fig_width: int = 10,
    fig_height: int = 10,
    dendrogram_ratio: float = 0.15,
    cbar_pos: tuple[float, float, float, float] = (0.02, 0.8, 0.05, 0.18),
    annotation: bool = False,
    annotation_fmt: str = ".3g",
) -> ClusterGrid:
    """Create clustermap of ANI matrix

    Parameters
    ----------
    ani_matrix_df : str | Path | pd.DataFrame
        ANI matrix file or dataframe
    linkage : np.ndarray | None, optional
        Scipy clustering result. If None, cluster in this method.
    cmap_colors : str | list[str] | None, optional
        cmap interpolation colors parameter. By default, `[lime, yellow, red]`
    cmap_ranges : str | list[str] | None, optional
        Range values (e.g. 80,90,95,100) for discrete cmap
    fig_width : int, optional
        Figure width
    fig_height : int, optional
        Figure height
    dendrogram_ratio : float, optional
        Dendrogram ratio
    cbar_pos : tuple[float, float, float, float], optional
        Colorbar position
    annotation : bool, optional
        If True, show ANI value annotation
    annotation_fmt : str, optional
        Annotation value format

    Returns
    -------
    g : ClusterGrid
        Cluster grid
    """
    if isinstance(ani_matrix_df, (str, Path)):
        ani_matrix_df = pd.read_csv(ani_matrix_df, sep="\t", encoding="utf-8")

    if linkage is None:
        linkage = hc.linkage(squareform(100 - ani_matrix_df), method="average")

    if cmap_colors is None:
        cmap_colors = ["lime", "yellow", "red"]
    if isinstance(cmap_colors, str):
        cmap_colors: list[str] = cmap_colors.split(",")

    if isinstance(cmap_ranges, str):
        cmap_ranges: list[float] = list(map(float, cmap_ranges.split(",")))
    if cmap_ranges is None:
        mycmap = LSC.from_list("mycmap", cmap_colors, gamma=cmap_gamma)
        opts = dict()
    else:
        N = len(cmap_ranges) - 1
        mycmap = LSC.from_list("mycmap", cmap_colors, gamma=cmap_gamma, N=N)
        opts = dict(norm=BoundaryNorm(cmap_ranges, N))
    mycmap.set_under("lightgrey")

    min_ani = min(filter(lambda v: v != 0, np.array(ani_matrix_df).flatten()))
    g: ClusterGrid = sns.clustermap(
        data=ani_matrix_df,
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
        },
        tree_kws={"linewidths": 1.5},
        **opts,  # type: ignore
    )
    return g

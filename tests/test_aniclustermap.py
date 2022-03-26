import os
from pathlib import Path

from aniclustermap import aniclustermap


def test_aniclustermap(genome_fasta_dir: Path, tmp_path: Path):
    """Check only no errors on runtime"""
    cpu_num = os.cpu_count()
    thread_num = cpu_num - 1 if cpu_num is not None and cpu_num != 1 else 1

    aniclustermap.run(
        indir=genome_fasta_dir,
        outdir=tmp_path,
        thread_num=thread_num,
        fig_width=10,
        fig_height=10,
        dendrogram_ratio=0.15,
        cmap_colors=["blue", "red"],
        cmap_gamma=1.0,
        annotation=False,
    )

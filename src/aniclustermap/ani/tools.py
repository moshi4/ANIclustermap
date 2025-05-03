from __future__ import annotations

import logging
import shlex
import subprocess as sp
import tempfile
from abc import ABC, abstractmethod
from pathlib import Path

import pandas as pd

from aniclustermap import const
from aniclustermap.ani.matrix import parse_ani_matrix

logger = logging.getLogger(__name__)


class AniCalcTool(ABC):
    """ANI Calculation Tool Abstract Base Class"""

    def __init__(
        self,
        genome_fasta_dir: str | Path,
        outdir: str | Path | None = None,
    ):
        """
        Parameters
        ----------
        genome_fasta_dir : str | Path
            Genome fasta directory
        outdir : str | Path | None, optional
            Output directory. If None, temporary directory is used.
        """
        self._genome_fasta_dir = Path(genome_fasta_dir)
        self._outdir = None if outdir is None else Path(outdir)

    @abstractmethod
    def run(self) -> pd.DataFrame:
        """Run ANI calculation tool"""
        raise NotImplementedError

    @abstractmethod
    def get_tool_name(self) -> str:
        """Get tool name"""
        raise NotImplementedError

    def _run_cmd(
        self,
        cmd: str,
        stdout_file: str | Path | None = None,
    ) -> None:
        """Run command

        Parameters
        ----------
        cmd : str
            Command to run
        stdout_file : str | Path | None, optional
            Write stdout result if file is set
        """
        logger.info(f"$ {cmd}")
        cmd_args = shlex.split(cmd)
        try:
            cmd_res = sp.run(cmd_args, capture_output=True, text=True, check=True)
            # Write stdout result if stdout_file is set
            if stdout_file:
                logger.info(f"> Save cmd stdout results to '{stdout_file}'")
                with open(stdout_file, "w", encoding="utf-8") as f:
                    f.write(cmd_res.stdout)
        except sp.CalledProcessError as e:
            returncode, stdout, stderr = e.returncode, str(e.stdout), str(e.stderr)
            logger.error(f"Failed to run command below ({returncode=})")
            logger.error(f"$ {cmd}")
            stdout_lines = stdout.splitlines()
            if len(stdout_lines) > 0:
                logger.error("STDOUT:")
                for line in stdout_lines:
                    logger.error(f"> {line}")
            stderr_lines = stderr.splitlines()
            if len(stderr_lines) > 0:
                logger.error("STDERR:")
                for line in stderr_lines:
                    logger.error(f"> {line}")
            raise
        except FileNotFoundError:
            name = self.get_tool_name()
            logger.error(f"{name} is not installed? Please check installation.")
            raise

    def _write_fasta_list(self, outfile: Path) -> int:
        """Write fasta file list for calculate ANI

        Parameters
        ----------
        outfile : Path
            Output fasta list file

        Returns
        -------
        count : int
            Number of target fasta file
        """
        # Search target fasta files from target directory
        target_exts = (".fa", ".fna", ".fna.gz", ".fasta")
        target_files = []
        for ext in target_exts:
            target_files.extend(self._genome_fasta_dir.glob(f"*{ext}"))

        if len(target_files) <= 1:
            raise ValueError("Number of input genome fasta file is less than 1.")

        # Write fasta list
        with open(outfile, "w") as f:
            f.write("\n".join(map(str, target_files)))

        genome_num = len(target_files)
        logger.info(f"Write list of {genome_num} genome fasta file path")
        logger.info(f"=> {outfile}")

        return genome_num


class FastAni(AniCalcTool):
    """FastANI ANI Calculation Class"""

    def get_tool_name(self) -> str:
        """Get tool name"""
        return "fastANI"

    def run(self, *, thread_num: int | None = None) -> pd.DataFrame:
        """Run fastANI

        Parameters
        ----------
        thread_num : int | None, optional
            Thread number. If None, `MaxThread-1` is set.

        Returns
        -------
        ani_matrix_df : pd.DataFrame
            ANI matrix dataframe
        """
        thread_num = const.DEFAULT_CPU if thread_num is None else thread_num
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) if self._outdir is None else self._outdir
            genome_fasta_list_file = outdir / "genome_fasta_list.txt"
            self._write_fasta_list(genome_fasta_list_file)
            ani_result_file = outdir / "fastani_result"
            cmd = f"{self.get_tool_name()} --ql {genome_fasta_list_file} --rl {genome_fasta_list_file} -o {ani_result_file} -t {thread_num} --matrix"  # noqa: E501
            self._run_cmd(cmd)
            ani_matrix_file = Path(f"{ani_result_file}.matrix")
            return parse_ani_matrix(ani_matrix_file)


class SkAni(AniCalcTool):
    """SkAni ANI Calculation Class"""

    def get_tool_name(self) -> str:
        """Get tool name"""
        return "skani"

    def run(self, *, thread_num: int | None = None) -> pd.DataFrame:
        """Run skani

        Parameters
        ----------
        thread_num : int | None, optional
            Thread number. If None, `MaxThread-1` is set.

        Returns
        -------
        ani_matrix_df : pd.DataFrame
            ANI matrix dataframe
        """
        thread_num = const.DEFAULT_CPU if thread_num is None else thread_num
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) if self._outdir is None else self._outdir
            genome_fasta_list_file = outdir / "genome_fasta_list.txt"
            self._write_fasta_list(genome_fasta_list_file)
            ani_matrix_file = outdir / "skani_result.matrix"
            cmd = f"{self.get_tool_name()} triangle -l {genome_fasta_list_file} -o {ani_matrix_file} -t {thread_num} -c 30"  # noqa: E501
            self._run_cmd(cmd)
            return parse_ani_matrix(ani_matrix_file)

# ANIclustermap: Draw Clustermap of All-vs-All ANI

![Python3](https://img.shields.io/badge/Language-Python3-steelblue)
![OS](https://img.shields.io/badge/OS-Windows_|_Mac_|_Linux-steelblue)
![License](https://img.shields.io/badge/License-MIT-steelblue)

## Installation

ANIclustermap is implemented in Python3.

Install PyPI stable version with pip:

    pip install aniclustermap

Install latest development version with pip:

    pip install git+https://github.com/moshi4/ANIclustermap.git

COGclassifier uses `fastANI` for All-vs-All ANI calculation.  
RPS-BLAST(v1.3.3) is bundled in ANIclustermap package ([src/aniclustermap/bin](https://github.com/moshi4/ANIclustermap/tree/main/src/aniclustermap/bin)).  

## Workflow

Description of ANIclustermap's automated workflow.

### 1. Calculate All-vs-All ANI by fastANI

### 2. Draw clustermap by Seaborn

## Usage

### Basic Command

    ANIclustermap -i [Genome fasta directory] -o [output directory]

### Options

    -h, --help            show this help message and exit
    -i , --indir          Input genome fasta directory
    -o , --outdir         Output directory
    -t , --thread_num     fastANI thread number parameter (Default: MaxThread - 1)
    -v, --version         Print version information

### Example Command

    COGclassifier -i ./example/input/ecoli.faa -o ./ecoli_cog_classifier

### Example API

```python
from aniclustermap import aniclustermap

genome_fasta_dir = "./example/input/"
outdir = "./ani_clustermap_result"
aniclustermap.run(genome_fasta_dir, outdir)
```

## Output Contents

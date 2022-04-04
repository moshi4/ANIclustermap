#!/usr/bin/bash
OUTDIR=output
mkdir -p $OUTDIR

ANIclustermap -i ./input/minimal_dataset -o ${OUTDIR}/01_minimal_dataset \
              --fig_width 8 --fig_height 5

ANIclustermap -i ./input/minimal_dataset -o ${OUTDIR}/02_minimal_dataset_annotation \
              --fig_width 8 --fig_height 5 --annotation

ANIclustermap -i ./input/small_dataset -o ${OUTDIR}/03_small_dataset \
              --fig_width 15

ANIclustermap -i ./input/small_dataset -o ${OUTDIR}/04_small_dataset_annotation \
              --fig_width 15 --annotation

ANIclustermap -i ./input/normal_dataset -o ${OUTDIR}/05_normal_dataset \
              --fig_width 15

ANIclustermap -i ./input/normal_dataset -o ${OUTDIR}/06_normal_dataset_annotation \
              --fig_width 20 --annotation

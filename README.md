# `HEAL: Homoeologous Exchange Automated Labeller`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)

A Snakemake workflow for automated labelling of homoeologous exchanges (HE).


## Wa

Assumes all read files have the same read length (same sequencing experiment)
Input desired bin size, else defaults to 10kb 
Use config file to set healr parameters. Else, defaults to default behaviour. 
We assume data is either paired or unpaired.
Assumes all short read files have same length of reads

If you wish to preserve reads that are softclipped in the bismark alignment, potentially to recover exact recombination breakpoints, you can do so by passing the config argument SOFT_CLIP=true. This will keep soft clipped reads in the alignment step (For which aligner??? Not both I think). In our experience this has not really worked but it could be done. For DNA seq data, softclipped reads are always allowed in bwa mem2 alignments. 
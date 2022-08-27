#!/bin/bash
set -o xtrace

# First, generate subsampled reads (if we need to)
# This validates md5sums and then runs seqtk to (stably) select 5M reads/paired end (max) = total 10M paired reads
snakemake --cores 100 --use-conda --printshellcmds --rerun-incomplete --keep-going --rerun-triggers mtime --until seqtk_subsample

# Next, run the subsampled analyses
# This runs the entire analysis pipeline, exclusively on _subsampled data
snakemake --cores 100 --use-conda --printshellcmds --rerun-incomplete --keep-going --rerun-triggers mtime --config subsampled=True

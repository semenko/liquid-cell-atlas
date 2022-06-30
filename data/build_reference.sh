#!/bin/bash
# Author: Nick Semenkovich <semenko@alum.mit.edu>
set -x

SCRIPT=$(realpath "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

#### Download iGenomes Reference
# NOTE: Doesn't support TLS :(
if [ ! -f "$SCRIPTPATH/reference/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa" ]
then
    wget --no-config --no-clobber --directory-prefix="$SCRIPTPATH/reference" http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/NCBI/GRCh38/Homo_sapiens_NCBI_GRCh38.tar.gz
    tar --keep-old-files -xvf "$SCRIPTPATH/reference/Homo_sapiens_NCBI_GRCh38.tar.gz" --directory "$SCRIPTPATH/reference"
fi


#### Build Bismark Reference

./bismark_genome_preparation --verbose ../liquid-cell-atlas/data/reference/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/ --genomic_composition --parallel 10
# ./bismark_genome_preparation --verbose ../liquid-cell-atlas/data/reference/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/
#### Download Biscuit Reference
# NOTE: We download the reference here, but you can also make this de-novo using their perl code, via:
# perl tools/build_biscuit_QC_assets.pl -v --outdir reference/biscuit --ref reference/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa

if [ ! -f "$SCRIPTPATH/reference/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.faT" ]
then
    wget --no-config --no-clobber --directory-prefix="$SCRIPTPATH/reference" https://github.com/huishenlab/biscuit/releases/download/v1.0.2.20220113/hg38_biscuit_qc_assets.zip
    unzip -n "$SCRIPTPATH/reference/hg38_biscuit_qc_assets.zip" -d "$SCRIPTPATH/reference"
fi

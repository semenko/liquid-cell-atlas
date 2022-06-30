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

### Download Biscuit Reference


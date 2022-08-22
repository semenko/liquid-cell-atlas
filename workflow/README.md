# Raw Sequencing Processing

A standardized WGBS pipeline for our internally generated WGBS bisulfite & EM-seq data. This goes from .fastq to methylation calls (via bwa-meth) and extensive QC and plotting, using a Snakemake pipeline.

## Background & Trimming Approach

I strongly suggest reading work from Felix Krueger (author of Bismark) as background. In particular:
- TrimGalore's [RRBS guide](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/RRBS_Guide.pdf)
- The Babraham [WGBS/RRBS tutorials](https://www.bioinformatics.babraham.ac.uk/training.html#bsseq)

For similar pipelines and inspiration, see:
- NEB's [EM-seq pipeline](https://github.com/nebiolabs/EM-seq/)
- Felix Krueger's [Nextflow WGBS Pipeline](https://github.com/FelixKrueger/nextflow_pipelines/blob/master/nf_bisulfite_WGBS)
- The Snakepipes [WGBS pipeline](https://snakepipes.readthedocs.io/en/latest/content/workflows/WGBS.html)

## Reference Genome

I chose GRCh38, with these specifics:
- No patches
- Includes the hs38d1 decoy
- Includes Alt chromosomes
- Applies the [U2AF1 masking file](https://genomeref.blogspot.com/2021/07/one-of-these-things-doest-belong.html)
- Applies the [Encode DAC exclusion](https://www.encodeproject.org/annotations/ENCSR636HFF/)

You can see a good explanation of the rationale for some of these components at [this NCBI explainer](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/README_analysis_sets.txt).

## Requirements

All software requirements are specified in `env.yaml` except for:
- NEB's [mark-nonconverted-reads.py package](https://github.com/nebiolabs/mark-nonconverted-reads)
- [wgbs_tools](https://github.com/nloyfer/wgbs_tools)


## Trimming Parameters
 
 --nextseq INT?: see https://sequencing.qcfail.com/articles/illumina-2-colour-chemistry-can-overcall-high-confidence-g-bases/
 This sets --nextseq-trim=3'CUTOFF within Cutadapt and ignores G base quality -- mutually exclusive with -q

Very stringend emseq clipping? https://github.com/FelixKrueger/Bismark/issues/419
--clip_R1 10 --clip_R2 10 ?
"""
For paired-end BS-Seq, it is recommended to remove the first few bp because the end-repair reaction may introduce a bias towards low methylation. Please refer to the M-bias plot section in the Bismark User Guide for some examples.
"""

`--clip_R1 5 --clip_R2 3 --three_prime_clip_R1 1 --three_prime_clip_R2 1`
--dovetail ?


## Pipeline Graph

Here's a high-level overview of the Snakemake pipeline, generated via `snakemake --rulegraph | dot -Tpng > rules.png`

<p align="center">
<img src="https://user-images.githubusercontent.com/167135/185484931-ccfa0549-6898-44e1-9be2-ee0cf25ee6b2.png" width="500">
</p>

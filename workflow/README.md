# Raw Sequencing Processing

A standardized WGBS pipeline for our internally generated WGBS bisulfite & EM-seq data. This goes from .fastq to methylation calls (via bwa-meth) and extensive QC and plotting, using a Snakemake pipeline.


## Output Data Structure

Raw data files from [data](../data) are processed and analyzed by this snakemake workflow. Within each project directory, the output is (roughly) structured as:

    SAMPLE_01/                  # e.g. melanoma / crc / healthy, etc.
    │   SAMPLE.bam              # The final alignment file 
    |   SAMPLE.bam.bai          #   (and its index)
    |── biscuit_qc/             # biscuit QC.sh text files
    |── epibeds/                # epibed files (bgzip-compressed & tabix-indexed)
    ├── fastp/                  # fastp statistics & logs
    ├── fastqc/                 # fastqc graphs 
    ├── goleft/                 # goleft coverage plots
    ├── logs/                   # runlogs from each pipeline component
    ├── methyldackel/           # mbias plots
    ├── raw/
    │   ├── ...fastq.gz         # Raw reads
    |   ├── ...md5.txt          # Checksums and validation
    ├── samtools/               # samtools statistics
    SAMPLE_02/
    ...
    ...
    multiqc/                    # A project-level multiqc stats across all data

Note each project also has a `_subsampled` directory with identical structure, which is the result of the pipeline run on only 10M reads/sample.


## Reference Genome

I chose GRCh38, with these specifics:
- No patches
- Includes the hs38d1 decoy
- Includes Alt chromosomes
- Applies the [U2AF1 masking file](https://genomeref.blogspot.com/2021/07/one-of-these-things-doest-belong.html)
- Applies the [Encode DAC exclusion](https://www.encodeproject.org/annotations/ENCSR636HFF/)

You can see a good explanation of the rationale for some of these components at [this NCBI explainer](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/README_analysis_sets.txt).

## Requirements

All software requirements are now specified in `env.yaml`.

You'll need `snakemake` to run the pipeline.

## Trimming Approach

I chose a relatively conservative approach to trimming -- which is needed due to end-repair bias, adaptase bias, and more. 

For **EMseq**, I trim 10 bp everywhere, from offline discussions with NEB. See some of my notes here: https://github.com/FelixKrueger/Bismark/issues/509

For **BSseq**, I trim 15 bp 5' R2, and 10 bp everywhere else due to significant adaptase bias.

For all reads, I set `--trim_poly_g` (see https://sequencing.qcfail.com/articles/illumina-2-colour-chemistry-can-overcall-high-confidence-g-bases/) and set a `--length_required` (minimum read length) of 10 bp.

## No Quality Filtering

Notably I do NOT do quality filtering here (I set `--disable_quality_filtering`), and save this for downstream analyses as desired.

I experimented with more stringent quality filtering early on, and found it had little yield / performance benefit. 

## Background & Inspiration

I strongly suggest reading work from Felix Krueger (author of Bismark) as background. In particular:
- TrimGalore's [RRBS guide](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/RRBS_Guide.pdf)
- The Babraham [WGBS/RRBS tutorials](https://www.bioinformatics.babraham.ac.uk/training.html#bsseq)

For similar pipelines and inspiration, see:
- NEB's [EM-seq pipeline](https://github.com/nebiolabs/EM-seq/)
- Felix Krueger's [Nextflow WGBS Pipeline](https://github.com/FelixKrueger/nextflow_pipelines/blob/master/nf_bisulfite_WGBS)
- The Snakepipes [WGBS pipeline](https://snakepipes.readthedocs.io/en/latest/content/workflows/WGBS.html)


## Pipeline Graph

Here's a high-level overview of the Snakemake pipeline (generated via `snakemake --rulegraph | dot -Tpng > rules.png`)

<p align="center">
<img src="https://user-images.githubusercontent.com/167135/216741276-7113ab3c-b7fc-42f6-b917-77a9c6b68398.png" width="500">
</p>

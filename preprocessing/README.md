# Raw Sequencing Preprocessing

A standardized WGBS pre-processing pipeline for our internally generated Bisulfite & EM-seq data. This adapts [this Biscuit sample snakemake pipline](https://github.com/vari-bbc/Biscuit_Snakemake_Workflow/).

We take an opinionated stance on that pipeline and rip out lots of more general components (e.g. no fastq_screen, preseq, etc.)
and enables lots of QC & filters (e.g. adapter trimming) by default. 

## Running:

1. Run `download_and_extract_reference.sh`
    1. This downloads and extracts the Hg38 reference genomes via [Illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html)
1. Run ``
    1. This establises a conda environment for snakemake
    1. For simplicity, we do **not** include module/cluster support
1. Test your installation by running ``
1. Run on real data buy running ``


## Processing Steps
  + [run once] Generate asset files used during QC related rules
  + [default off] Modify and index genome reference to including methylation controls (?)
  + Trim adapters and/or hard clip R2
  + Run FastQC
  + Alignment, duplicate tagging, indexing, flagstat of input data (biscuitBlaster v1 and v2)
  + Methylation information extraction (BED Format)
  + Merge C and G beta values in CpG dinucleotide context
  + [default off] SNP and Epiread extraction
  + MultiQC with BICUIT QC modules specifically for methyaltion data
  + Generate plots of the observed / expected coverage ratio for different genomic features
  + Generate percentage of covered CpGs and CpG island coverage figures
  + [default off] QC methylated and unmethylated controls


## License
GPLv3

NOTE: Preprocessing must be GPLv3 licensed, as it derives from the Biscuit sample WGBS pipeline noted above.

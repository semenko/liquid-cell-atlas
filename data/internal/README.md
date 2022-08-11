# Internal Experiments

A collection of our internally generated sequencing data. Raw files aren't stored on github (internally, they're symlinked / hardlinked from different volumes).
  
These are parsed by snakemake and the preprocessing/analytical pipelines.

  
| Folder | Experiment Description |
| --- | --- |
| epcam_spike | PBMCs from healthy donors, spiked with varying percentages of HT-29 DNA (an EPCAM+ CRC cell line), sequenced with EM-seq and BS-seq. |
| healthy_pbmc | PBMCs from healthy donors, sequenced with EM-seq and BS-seq. |
| melanoma | Melanoma patient samples as part of a collaboration with Yale, including patients undergoing immunotherapy. Sequence with EM-seq and BS-seq. |
| crc | Colorectal tumor samples, including tumor, cfDNA, and sorted populations. |


## WashU / Chaudhuri Lab Users

The raw data backing these files are stored in RIS in Active/lca_data (e.g. `/storage1/fs1/aadel/Active/lca_data`), which is periodically mirrored to our local compute infrastructure (currently on aclm350 in `/logo2/lca_data`).

## Notes

Within the melanoma data, a few of the original .fastq files from the first batch (batch B00) were misplaced. The raw reads were recovered from .bam files as interleaved fastq files (to preserve all flagged reads & orphans) using `samtools fastq <BAM> | gzip --fast > <OUT>.fastq.gz`

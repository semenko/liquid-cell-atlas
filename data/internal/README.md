# Internal Experiments

A collection of our internally generated sequencing data. Raw files aren't stored on github (internally, they're symlinked / hardlinked from different volumes).
  
These are parsed by snakemake and the preprocessing/analytical pipelines.

  
| Folder | Experiment Description |
| --- | --- |
| epcam_spike | PBMCs from healthy donors, spiked with varying percentages of HT-29 DNA (an EPCAM+ CRC cell line), sequenced with EM-seq and BS-seq. |
| healthy_pbmc | PBMCs from healthy donors, sequenced with EM-seq and BS-seq. |
| melanoma | Melanoma patient samples as part of a collaboration with Yale, including patients undergoing immunotherapy. Sequence with EM-seq and BS-seq. |
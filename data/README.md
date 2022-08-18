# Sequencing Data

A collection of our sequencing data. Raw files aren't stored on github (internally, they're symlinked / mounted from different volumes). Our published data can be download via *** and in the public domain. Third party data (e.g. ENCODE and BLUEPRINT) are covered by their own 3rd party licenses.

These directories are parsed by snakemake and the preprocessing/analytical pipelines.
  
| Folder | Experiment Description |
| --- | --- |
| _external | External data from Blueprint & Encode. |
| epcam_spike | PBMCs from healthy donors, spiked with varying percentages of HT-29 DNA (an EPCAM+ CRC cell line), sequenced with EM-seq and BS-seq. |
| healthy_pbmc | PBMCs from healthy donors, sequenced with EM-seq and BS-seq. |
| melanoma | Melanoma patient samples as part of a collaboration with Yale, including patients undergoing immunotherapy. Sequence with EM-seq and BS-seq. |
| crc | Colorectal tumor samples, including tumor, cfDNA, and sorted populations. |


## Sample Definition Files

Each internal experiment has an associated .csv with key experimental metadata and information. Most samples were sequenced by [MedGenome](https://research.medgenome.com/), though some were sequenced by the [McDonnell Genome Institute](https://www.genome.wustl.edu/). Of note, the "batch" in each sample is arbitrary -- and just a unique identifier for samples run on the same batch (e.g. batch B00 was not necessarily sequenced before or after batch B01). 

## WashU / Chaudhuri Lab Users

The raw data backing these files are stored in RIS in Active/lca_data (e.g. `/storage1/fs1/aadel/Active/lca_data`), which is periodically mirrored to our local compute infrastructure (currently on aclm350 in `/logo2/lca_data`).

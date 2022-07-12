# Raw Sequencing Preprocessing

A standardized WGBS pre-processing pipeline for our internally generated WGBS bisulfite & EM-seq data. This goes from .fastq to methylation calls via Bismark, using Snakemake.

## Background

I strongly suggest reading work from Felix Krueger (author of Bismark) to understand this approach. In particular:
- TrimGalore's [RRBS guide](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/RRBS_Guide.pdf)
- The Babraham [WGBS/RRBS tutorials](https://www.bioinformatics.babraham.ac.uk/training.html#bsseq)
- Felix Krueger's [Nextflow WGBS Pipeline](https://github.com/FelixKrueger/nextflow_pipelines/blob/master/nf_bisulfite_WGBS)

## Processing Steps:

- QC Analysis with FastQC
- Contamination testing with FastQ Screen
- Adapter & 3' trimming with TrimGalore
    - Note: I chose this over fastp & others as it was specifically designed for WGBS/RRBS (and mostly just wraps cutadapt). There is great documentation on its filtering/trimming decisions [in this guide](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/RRBS_Guide.pdf).
    - TrimGalore steps:
        - Trim low quality (Phred <20) basepairs from 3' end
        - Remove adapters (first 13bp of Illumina standard adapters)
        - Min read size 20bp (for each entry in a pair)
- Run Bismark


```
-- [Once] Generate Bismark index for GRCh38/hg38
-- FastQC Analysis
-- Trim Galore (default flags)
    |
    -- FastQ Screen (--bisulfite)
    -- FastQC (again)
    -- Bismark (default flags)
        |
        -- Bismark deduplication (? EmSEQ necessary?)
            |
            -- Bismark methylation extraction (--bedGraph --buffer 10G --parallel 4 --ignore_r2 2)       
-- bismark2report
-- bismark2summary
-- MultiQC
```

### Random thoughts

Trim galore:
`trim_galore --paired --quality 20 --phred33 --illumina --stringency 1 -e 0.1 --gzip --length 20 --output_dir WHAT --clip_R2 5? --cores 4 `
 
 --nextseq INT?: see https://sequencing.qcfail.com/articles/illumina-2-colour-chemistry-can-overcall-high-confidence-g-bases/
 This sets --nextseq-trim=3'CUTOFF within Cutadapt and ignores G base quality -- mutually exclusive with -q

Very stringend emseq clipping? https://github.com/FelixKrueger/Bismark/issues/419
--clip_R1 10 --clip_R2 10 ?
"""
For paired-end BS-Seq, it is recommended to remove the first few bp because the end-repair reaction may introduce a bias towards low methylation. Please refer to the M-bias plot section in the Bismark User Guide for some examples.
"""

`--clip_R1 5 --clip_R2 3 --three_prime_clip_R1 1 --three_prime_clip_R2 1`
--dovetail ?

## Running:

1. Run `download_and_extract_reference.sh`
    1. This downloads and extracts the Hg38 reference genomes via [Illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html)
1. Run ``
    1. This establises a conda environment for snakemake
    1. For simplicity, we do **not** include module/cluster support
1. Test your installation by running ``
1. Run on real data buy running ``



## License
MIT - derives some code from https://github.com/lyijin/bismsmark which is MIT licensed.

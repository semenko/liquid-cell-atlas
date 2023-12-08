#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##
# A command-line Python script to extract methylation data from an aligned .bam file.
# This assumes alignment with Biscuit, or with gem3-mapper and gemBS, though it should work with any aligner.
#
# TODO: Test with other aligners.
#
# Note this is designed for readability & ease of use (not speed) since it's a one-time use script.
#
#
# Copyright (c) 2023 Nick Semenkovich <semenko@alum.mit.edu>.
#   https://nick.semenkovich.com/
#
# This software is released under the MIT License:
#  <https://opensource.org/licenses/mit-license.php>


# Import modules
import argparse
import json
import os
import sys
import time

# Third party modules
import numpy as np
import pysam
import scipy.sparse
from tqdm import tqdm

# We use SeqIO to parse the reference .fa file
from Bio import SeqIO


## Globals
CHROMOSOMES = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
# dict: {chromosome: index}, e.g. {"chr1": 0, "chr2": 1, ...}
CHROMOSOMES_DICT = {ch: idx for idx, ch in enumerate(CHROMOSOMES)}

def get_cpg_sites_from_fasta(reference_fasta, verbose=False, skip_cache=False):
    """
    Generate a dict of *all* CpG sites across each chromosome in the reference genome.

    This is a dict of lists, where the key is the chromosome: e.g. "chr1"
    The value is a list of CpG sites: e.g. [0, 35, 190, 212, 1055, ...]

    We store this as a dict because it's easier to portably serialize to disk as JSON.

    Args:
        reference_fasta (str): Path to the reference genome .fa file.
        verbose (bool): Print verbose output.
        skip_cache (bool): Ignore any cached files (slow!).
    Returns:
        cpg_sites_dict (dict): A dict of CpG sites for each chromosome in the reference genome.
    """
    # TODO: Store hash/metadata / reference file info, etc.?
    cached_cpg_sites_json = os.path.splitext(reference_fasta)[0] + ".cpg_all_sites.json"

    if verbose:
        print(f"\nLoading all CpG sites for: {reference_fasta}")

    if os.path.exists(cached_cpg_sites_json) and not skip_cache:
        if verbose:
            print(f"\tLoading all CpG sites from cache: {cached_cpg_sites_json}")
        with open(cached_cpg_sites_json, "r", encoding="utf-8") as f:
            cpg_sites_dict = json.load(f)

        return cpg_sites_dict

    if verbose:
        print(f"\tNo cache of all CpG sites found (or --skip-cache=True), generating from: {reference_fasta}")

    # We want a list of all CpG sites in the reference genome, including the chromosome and position
    # We'll use this to define our embedding size
    # Each chromosome is identified with a line starting with ">" followed by the chromosome name

    # Store the CpG sites in a dict per chromosome
    cpg_sites_dict = {}

    # Iterate over sequences
    for seqrecord in tqdm(SeqIO.parse(reference_fasta, "fasta"), total=len(CHROMOSOMES), disable=not verbose):
        if seqrecord.id not in CHROMOSOMES:
            tqdm.write(f"\tSkipping chromosome {seqrecord.id}")
            continue
        sequence = seqrecord.seq

        # Find all CpG sites
        # The pos+1 is because we want to store the 1-based position, because .bed is wild and arguably 1-based maybe:
        # e.g. https://genome-blog.soe.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/
        # Regardless, $ samtools faidx GRCh38-DAC-U2AF1.fna chr1:0-0 throws an error, while
        # $ samtools faidx GRCh38-DAC-U2AF1.fna chr1:1-1 returns the correct base.

        cpg_indices = []
        search_str = "CG"
        pos = sequence.find(search_str)

        while pos != -1:
            cpg_indices.append(pos + 1)
            pos = sequence.find(search_str, pos + 1)

        cpg_sites_dict[seqrecord.id] = cpg_indices

    if verbose:
        print(f"\tSaving all cpg sites to cache: {cached_cpg_sites_json}")
    with open(cached_cpg_sites_json, "w", encoding="utf-8") as f:
        json.dump(cpg_sites_dict, f)

    return cpg_sites_dict


# TODO: Ponder this window size, as aligned reads might be larger by a bit... Is this useful?
def get_windowed_cpg_sites(reference_fasta, cpg_sites_dict, window_size, verbose=False, skip_cache=False):
    """
    Generate a dict of CpG sites for each chromosome in the reference genome.

    This is a dict of lists, where but each list contains a tuple of CpG ranges witin a window
    The key is the chromosome: e.g. "chr1"
    The value is a list of CpG sites: e.g. [(0, 35), (190, 212), (1055, ...)]

    We store this as a dict because it's easier to portably serialize to disk as JSON.

    NOTE: The window size is tricky -- it's set usually to the read size, but reads can be larger than the window size, 
    since they can map with deletions, etc. So we need to be careful here.

    Args:
        reference_fasta (str): Path to the reference genome .fa file.
        window_size (int): Size of the window to use.
        verbose (bool): Print verbose output.
        skip_cache (bool): Ignore any cached files (slow!).
    Returns:
        windowed_cpg_sites_dict (dict): A dict of CpG sites for each chromosome in the reference genome.
        windowed_cpg_sites_dict_reverse (dict): A dict of each per-window CpG site. # TODO: Clarify this & rename var?
    """

    # Check if we have a cached version of our windowed cpg sites
    # TODO: Update caching to use a hash of the reference genome, window size, etc.
    windowed_cpg_sites_cache = os.path.splitext(reference_fasta)[0] + ".cpg_windowed_sites.json"
    windowed_cpg_sites_reverse_cache = os.path.splitext(reference_fasta)[0] + ".cpg_windowed_sites_reverse.json"

    if os.path.exists(windowed_cpg_sites_cache) and os.path.exists(windowed_cpg_sites_reverse_cache) and not skip_cache:
        if verbose:
            print(f"\tLoading windowed CpG sites from caches:\n\t\t{windowed_cpg_sites_cache}\n\t\t{windowed_cpg_sites_reverse_cache}")
        with open(windowed_cpg_sites_cache, "r", encoding="utf-8") as f:
            windowed_cpg_sites_dict = json.load(f)
        with open(windowed_cpg_sites_reverse_cache, "r", encoding="utf-8") as f:
            # This wild object_hook is to convert the keys back to integers, since JSON only supports strings as keys
            windowed_cpg_sites_dict_reverse = json.load(f, object_hook=lambda d: {int(k) if k.isdigit() else k: v for k, v in d.items()})

        return windowed_cpg_sites_dict, windowed_cpg_sites_dict_reverse

    if verbose:
        print(f"\tNo cache found (or --skip-cache=True) at: {windowed_cpg_sites_cache} or {windowed_cpg_sites_reverse_cache}")

    # Let's generate ranges (windows) of all CpGs within READ_SIZE of each other
    # This is to subsequently query our .bam/.sam files for reads containing CpGs
    # TODO: Shrink read size a little to account for trimming?
    assert 10 < window_size < 500, "Read size is atypical, please double check (only validated for ~150 bp.)"

    if verbose:
        print(f"\n\tGenerating windowed CpG sites dict (window size = {window_size} bp.)\n")

    # This is a dict of lists, where but each list contains a tuple of CpG ranges witin a window
    # Key: chromosome, e.g. "chr1"
    # Value: a list of tuples, e.g. [(0,35), (190,212), (1055,)]
    windowed_cpg_sites_dict = {}

    # And a reverse dict of dicts where chrom->window_start->[cpgs]
    windowed_cpg_sites_dict_reverse = {}

    # Loop over all chromosomes
    for chrom, cpg_sites in tqdm(cpg_sites_dict.items(), disable=not verbose):
        windowed_cpg_sites_dict[chrom] = []
        windowed_cpg_sites_dict_reverse[chrom] = {}
        window_start = None
        window_end = None

        # Loop over all CpG sites
        cpg_sites_len = len(cpg_sites)
        temp_per_window_cpg_sites = []
        for i in range(cpg_sites_len):
            temp_per_window_cpg_sites.append(cpg_sites[i])

            if window_start is None:
                window_start = cpg_sites[i]

            # If we're at the end of the chromosome or the next CpG site is too far away
            if i + 1 == cpg_sites_len or (cpg_sites[i + 1] - cpg_sites[i]) > window_size:
                # We have a complete window
                window_end = cpg_sites[i]
                windowed_cpg_sites_dict[chrom].append((window_start, window_end))
                windowed_cpg_sites_dict_reverse[chrom][window_start] = temp_per_window_cpg_sites
                temp_per_window_cpg_sites = []
                window_start = None
                window_end = None

    # Save these to .json caches
    if verbose:
        print(f"\tSaving windowed CpG sites to caches:\n\t\t{windowed_cpg_sites_cache}\n\t\t{windowed_cpg_sites_reverse_cache}")

    with open(windowed_cpg_sites_cache, "w", encoding="utf-8") as f:
        json.dump(windowed_cpg_sites_dict, f)
    with open(windowed_cpg_sites_reverse_cache, "w", encoding="utf-8") as f:
        json.dump(windowed_cpg_sites_dict_reverse, f)

    return windowed_cpg_sites_dict, windowed_cpg_sites_dict_reverse


# TODO: Object orient this input / simplify the input?
# TODO: Ingest chr_to_cpg_to_embedding_dict instead?
def embedding_to_genomic_position(total_cpg_sites, cpg_sites_dict, cpgs_per_chr_cumsum, embedding_pos):
    """
    Given an embedding position, return the chromosome and position.

    Parameters
    ----------
    embedding_pos : int
        The embedding position.

    Returns
    -------
    tuple:
        The chromosome and position, e.g. ("chr1", 12345)
    """
    assert embedding_pos >= 0 and embedding_pos < total_cpg_sites

    # Find the index of the first element in cpgs_per_chr_cumsum that is greater than or equal to the embedding position
    chr_index = np.searchsorted(cpgs_per_chr_cumsum, embedding_pos + 1, side="left")

    # Now we know the chromosome, but we need to find the position within the chromosome
    # If this is the first chromosome, the position is just the embedding position
    if chr_index == 0:
        return CHROMOSOMES[chr_index], cpg_sites_dict[CHROMOSOMES[chr_index]][embedding_pos]
    # Otherwise, subtract the length of the previous chromosomes from the embedding position
    embedding_pos -= cpgs_per_chr_cumsum[chr_index - 1]
    return CHROMOSOMES[chr_index], cpg_sites_dict[CHROMOSOMES[chr_index]][embedding_pos]


# TODO: Object orient this input / simplify the input?
def genomic_position_to_embedding(chr_to_cpg_to_embedding_dict, cpgs_per_chr_cumsum, chrom, pos):
    """
    Given a genomic position, return the embedding position.

    Parameters
    ----------
    chrom : str
        The chromosome, e.g. "chr1"
    pos : int
        The position, e.g. 12345

    Returns
    -------
    embedding_pos : int
        The numerical CpG embedding position, e.g. 3493
    """
    # Find the index of the chromosome
    chr_index = CHROMOSOMES_DICT[chrom]
    # Find the index of the CpG site in the chromosome
    cpg_index = chr_to_cpg_to_embedding_dict[chrom][pos]
    # If this is the first chromosome, the embedding position is just the CpG index
    if chr_index == 0:
        return cpg_index
    # Otherwise, add the length of the previous chromosomes to the CpG index
    return cpg_index + cpgs_per_chr_cumsum[chr_index - 1]


def extract_methylation_data_from_bam(input_bam, total_cpg_sites, chr_to_cpg_to_embedding_dict, cpgs_per_chr_cumsum, windowed_cpg_sites_dict, windowed_cpg_sites_dict_reverse, quality_limit=0, verbose=False, debug=False):
    """
    Extract methylation data from a .bam file.

    Args:

    """
    input_bam_object = pysam.AlignmentFile(input_bam, "rb", require_index=True, threads=1)

    if verbose:
        print(f"\tTotal reads: {input_bam_object.mapped:,}\n")

    # Temporary storage for loading into a COO matrix
    # For a COO, we want three lists:
    # row: the read number (we'll store a read name -> ID dict perhaps?)
    # column: cpg #
    # data: methylation state (1 [methylated], 0 [unmethylated], and -1 [no data/snv/indel])
    next_coo_row_number = 0
    read_name_to_row_number = {}

    # TODO: Modularize/clean?
    coo_row = []
    coo_col = []
    coo_data = []

    # This is slow, but we only run it once and store the results for later
    for chrom, windowed_cpgs in tqdm(windowed_cpg_sites_dict.items(), disable=not verbose):
        for start_pos, stop_pos in windowed_cpgs:
            cpgs_within_window = windowed_cpg_sites_dict_reverse[chrom][start_pos]
            if debug:
                print(f"Window Position: {start_pos} - {stop_pos}")
                print(f"\tCovered CpGs: {cpgs_within_window}")
            # cpgs_within_window = [10469]
            # This returns an AlignmentSegment object from PySam

            # TODO: Think carefully here -- A read could span more than one segment, right?!?
            for aligned_segment in input_bam_object.fetch(contig=chrom, start=start_pos, end=stop_pos):
                if aligned_segment.mapping_quality < quality_limit:
                    continue
                if aligned_segment.is_duplicate:
                    continue
                if aligned_segment.is_qcfail:
                    continue
                if aligned_segment.is_secondary:
                    continue
                
                if debug:
                    print(aligned_segment.query_name)
                    # Validity tests
                    assert aligned_segment.is_mapped
                    assert aligned_segment.is_supplementary is False
                    assert aligned_segment.reference_start <= stop_pos
                    assert aligned_segment.reference_end >= start_pos

                    # Ensure alignment methylation tags exist
                    assert aligned_segment.has_tag("MD")  # Location of mismatches (methylation)
                    assert aligned_segment.has_tag("YD")  # Bisulfite conversion strand label (f: OT/CTOT C->T or r: OB/CTOB G->A)

                # TODO: We ignore paired/unpaired read status for now, and treat them the same
                # Should we treat paired reads / overlapping reads differently?

                bisulfite_parent_strand_is_reverse = None
                if aligned_segment.has_tag("YD"):  # Biscuit tag
                    yd_tag = aligned_segment.get_tag("YD")
                    if yd_tag == "f":  # Forward = C→T
                        # This read derives from OT/CTOT strand: C->T substitutions matter (C = methylated, T = unmethylated),
                        bisulfite_parent_strand_is_reverse = False
                    elif yd_tag == "r":  # Reverse = G→A
                        # This read derives from the OB/CTOB strand: G->A substitutions matter (G = methylated, A = unmethylated)
                        bisulfite_parent_strand_is_reverse = True
                elif aligned_segment.has_tag("XB"):  # gem3 / blueprint tag
                    xb_tag = aligned_segment.get_tag("XB")  # XB:C = Forward / Reference was CG
                    # TODO: Double check one or two of these gem3 tags manually.
                    if xb_tag == "C":
                        bisulfite_parent_strand_is_reverse = False
                    elif xb_tag == "G":  # XB:G = Reverse / Reference was GA
                        bisulfite_parent_strand_is_reverse = True

                # TODO: Think about this old note of mine -- is this correct?
                # We have paired-end reads; one half should (the "parent strand") has the methylation data.
                # The other half (the "daughter strand") was the complement created by PCR, which we don't care about.
                if bisulfite_parent_strand_is_reverse != aligned_segment.is_reverse:
                    # Skip if we're not on the bisulfite-converted parent strand.
                    if debug:
                        print("\tNot on methylated strand, ignoring.")
                    continue

                # get_aligned_pairs returns a list of tuples of (read_pos, ref_pos)
                # We filter this to only include the specific CpG sites from above
                this_segment_cpgs = [e for e in aligned_segment.get_aligned_pairs(matches_only=True) if e[1] + 1 in cpgs_within_window]

                # Ok we're on the same strand as the methylation (right?)
                # Let's compare the possible CpGs in this interval to the reference and note status
                #   A methylated C will be *unchanged* and read as C (pair G)
                #   An unmethylated C will be *changed* and read as T (pair A)
                for query_pos, ref_pos in this_segment_cpgs:
                    query_base = aligned_segment.query_sequence[query_pos]
                    # query_base2 = aligned_segment.get_forward_sequence()[query_pos] # raw off sequencer
                    # query_base3 = aligned_segment.query_alignment_sequence[query_pos] # this needs to be offset by the soft clip

                    # Store in our sparse array
                    coo_row.append(next_coo_row_number)
                    # instead of ref_pos, we want the index of the CpG in the list of CpGs
                    # coo_col.append(ref_pos)

                    # TODO: Object orient these inputs? -- lots of bad inheritence style here
                    coo_col.append(genomic_position_to_embedding(chr_to_cpg_to_embedding_dict, cpgs_per_chr_cumsum, chrom, ref_pos + 1))
                    # coo_data.append(1 if query_base == "C" else 0)
                    if query_base == "C":
                        # Methylated
                        coo_data.append(1)
                        if debug:
                            print(f"\t{query_pos} {ref_pos} C->{query_base} [Methylated]")
                    elif query_base == "T":
                        coo_data.append(0)
                        # Unmethylated
                        if debug:
                            print(f"\t{query_pos} {ref_pos} C->{query_base} [Unmethylated]")
                    else:
                        coo_data.append(-1)  # or just 0?
                        if debug:
                            print(f"\t{query_pos} {ref_pos} C->{query_base} [Unknown! SNV? Indel?]")

                assert aligned_segment.query_name not in read_name_to_row_number

                read_name_to_row_number[aligned_segment.query_name] = next_coo_row_number
                next_coo_row_number += 1
                if debug:
                    print("************************************************\n")

                # query_bp = aligned_segment.query_sequence[pileupread.query_position]
                # reference_bp = aligned_segment.get_reference_sequence()[aligned_segment.reference_start - pileupcolumn.reference_pos].upper()


    ## IIRC there's still a critical edge here, where sometimes we raise ValueError('row index exceeds matrix dimensions')

    # Validity checks -- TODO: move elsewhere? to loader?
    assert max(coo_data) <= 1
    assert min(coo_data) >= -1

    print("Debug info for coo_matrix dimensions:")
    print(f"\tcoo_row: {len(coo_row):,}")
    print(f"\tcoo row max: {max(coo_row):,}")
    print(f"\tcoo_col: {len(coo_col):,}")
    print(f"\tcoo col max: {max(coo_col):,}")
    print(f"\tcoo_data: {len(coo_data):,}")
    print(f"\tlen(read_name_to_row_number): {len(read_name_to_row_number):,}")
    print(f"\ttotal_cpg_sites: {total_cpg_sites:,}")

    # The size of the coo_matrix is:
    #   Number of rows = number of reads that pass our filters
    #   Number of columns = number of CpG sites

    return scipy.sparse.coo_matrix((coo_data, (coo_row, coo_col)), shape=(next_coo_row_number, total_cpg_sites))

    # return scipy.sparse.coo_matrix((coo_data, (coo_row, coo_col)), shape=(len(read_name_to_row_number) + 1, total_cpg_sites))


def main():
    """
    Our main function. Call argparse and whatnot.

    Options are:
    --input-path: Input .bam file OR directory to recursively process.
    --reference-fasta: Reference genome fasta file (critical to determine CpG sites).
    --verbose: Verbose output.
    --paranoid: Paranoid mode (extensive validity checking).
    # --output-file: Output file. Default is to use the extension .methylation.npz. (Only available for single .bam file input.)
    --overwrite: Overwrite output file if it exists.

    Returns: Nothing.
    """

    # Let's do our argparse stuff
    parser = argparse.ArgumentParser(description="Extract read-level methylation data from an aligned .bam file and store as a numpy array.")
    parser.add_argument("--input-path", help="Input .bam file OR directory to recursively process.", required=True)
    parser.add_argument("--reference-fasta", help="Reference genome fasta file (critical to determine CpG sites).", required=True)
    parser.add_argument("--quality-limit", help="Quality filter for aligned reads (default = 20)", default=20, type=int)
    parser.add_argument("--window-size", help="Window size (default = 150)", default=150, type=int)
    parser.add_argument("--verbose", help="Verbose output.", action="store_true")
    parser.add_argument("--skip-cache", help="De-novo generate CpG sites (slow).", action="store_true")
    parser.add_argument("--debug", help="Debug mode (extensive validity checking + debug messages).", action="store_true")
    # parser.add_argument("--output-file", help="Output file. Default is to use the extension .methylation.npz. (Only available for single .bam file input.)", default=None)
    parser.add_argument("--overwrite", help="Overwrite output file if it exists.", action="store_true")

    # Parse the arguments & set up some variables
    args = parser.parse_args()
    input_path = args.input_path
    reference_fasta = args.reference_fasta
    verbose = args.verbose
    quality_limit = args.quality_limit
    skip_cache = args.skip_cache
    debug = args.debug
    # output_file = args.output_file
    window_size = args.window_size

    time_start = time.time()
    # Print run information
    if verbose:
        print(f"Reference fasta: {reference_fasta}")
        print(f"Input path: {input_path}")

    # Check if input_path is a file or a directory
    if os.path.isfile(input_path):
        bams_to_process = [input_path]
    elif os.path.isdir(input_path):
        # Recursively find all .bam files in this path, and add them to a list
        bams_to_process = []
        for root, _, files in os.walk(input_path):
            for file in files:
                if file.endswith(".bam"):
                    bams_to_process.append(os.path.join(root, file))
    else:
        raise ValueError(f"Input path {input_path} is not a file or a directory.")

    if verbose:
        print(f"Found {len(bams_to_process)} .bam file(s) to process.")
       
    #if args.output_file is None:
    #    output_file = os.path.splitext(input_bam)[0] + ".methylation.npz"

    # Check input/output validity
    for bam_file in bams_to_process:
        assert os.access(bam_file, os.R_OK), f"Input file is not readable: {bam_file}"

        output_file = os.path.splitext(bam_file)[0] + ".methylation.npz"
        # Check if the output files exist or are writable
        if os.path.exists(output_file):
            if args.overwrite:
                if verbose:
                    print(f"\tOutput file exists and --overwrite specified. Will overwrite: {output_file}")
                assert os.access(output_file, os.W_OK), f"Output file is not writable: {output_file}"
            else:
                print(f"Exiting. An output file exists and --overwrite not specified: {output_file}")
                sys.exit(1)
        # Otherwise, check the path is writable
        else:
            assert os.access(os.path.dirname(os.path.abspath(output_file)), os.W_OK), f"Output file path is not writable: {output_file}"
    # We need to obtain all cpg sites in the reference genome
    # NOTE: This can take a while (~10 minutes for GRCh38 if not cached)
    cpg_sites_dict = get_cpg_sites_from_fasta(reference_fasta=reference_fasta, verbose=verbose, skip_cache=skip_cache)

    # How many CpG sites are there?
    total_cpg_sites = sum([len(v) for v in cpg_sites_dict.values()])
    if verbose:
        print(f"Total CpG sites: {total_cpg_sites:,}")

    assert total_cpg_sites > 28_000_000  # Validity check for hg38

    # Create a dictionary of chromosome -> CpG site -> embedding index for efficient lookup
    chr_to_cpg_to_embedding_dict = {ch: {cpg: idx for idx, cpg in enumerate(cpg_sites_dict[ch])} for ch in CHROMOSOMES}

    # Count the number of CpGs per chromosome
    cpgs_per_chr = {k: len(v) for k, v in cpg_sites_dict.items()}

    # Add up the number of CpGs per chromosome, e.g. chr1, then chr1+chr2, then chr1+chr2+chr3, etc
    cpgs_per_chr_cumsum = np.cumsum([cpgs_per_chr[k] for k in CHROMOSOMES])

    # TODO: Move these into to a formal test framework
    # TODO: Simplify the input framework (likely object orient the window / cpg dict?)
    # FYI embedding_to_genomic_position is unused currently :P
    ### Tests
    assert cpgs_per_chr_cumsum[-1] == total_cpg_sites
    assert embedding_to_genomic_position(total_cpg_sites, cpg_sites_dict, cpgs_per_chr_cumsum, 0) == ("chr1", cpg_sites_dict["chr1"][0])
    assert embedding_to_genomic_position(total_cpg_sites, cpg_sites_dict, cpgs_per_chr_cumsum, 1) == ("chr1", cpg_sites_dict["chr1"][1])
    # Edges
    assert embedding_to_genomic_position(total_cpg_sites, cpg_sites_dict, cpgs_per_chr_cumsum, cpgs_per_chr_cumsum[0]) == ("chr2", cpg_sites_dict["chr2"][0])
    assert embedding_to_genomic_position(total_cpg_sites, cpg_sites_dict, cpgs_per_chr_cumsum, cpgs_per_chr_cumsum[-1] - 1) == ("chrY", cpg_sites_dict["chrY"][-1])

    ### Tests
    assert genomic_position_to_embedding(chr_to_cpg_to_embedding_dict, cpgs_per_chr_cumsum, "chr1", cpg_sites_dict["chr1"][0]) == 0
    assert genomic_position_to_embedding(chr_to_cpg_to_embedding_dict, cpgs_per_chr_cumsum, "chr1", cpg_sites_dict["chr1"][1]) == 1
    # Edges
    assert genomic_position_to_embedding(chr_to_cpg_to_embedding_dict, cpgs_per_chr_cumsum, "chr2", cpg_sites_dict["chr2"][0]) == cpgs_per_chr_cumsum[0]
    assert genomic_position_to_embedding(chr_to_cpg_to_embedding_dict, cpgs_per_chr_cumsum, "chrY", cpg_sites_dict["chrY"][-1]) == cpgs_per_chr_cumsum[-1] - 1
    ########



    # Get our windowed_cpg_sites, hopefully cached!
    if verbose:
        print("\nLoading (or generating) windowed CpG sites for reference genome.")

    windowed_cpg_sites_dict, windowed_cpg_sites_dict_reverse = get_windowed_cpg_sites(reference_fasta=reference_fasta, cpg_sites_dict=cpg_sites_dict, window_size=window_size, verbose=verbose, skip_cache=skip_cache)

    if verbose:
        print(f"\nTime elapsed: {time.time() - time_start:.2f} seconds")

    #################################################
    # Operate over the input BAM files
    #################################################

    for i, input_bam in enumerate(bams_to_process):
        time_bam = time.time()
        output_file = os.path.splitext(input_bam)[0] + ".methylation.npz"
        # Extract methylation data as a COO sparse matrix
        if verbose:
            print("\n" + "=" * 80)
            print(f"Processing BAM file {i+1} of {len(bams_to_process)}")
            print(f"\nExtracting methylation data from: {input_bam}")

        methylation_data_coo = extract_methylation_data_from_bam(input_bam=input_bam, total_cpg_sites=total_cpg_sites, chr_to_cpg_to_embedding_dict=chr_to_cpg_to_embedding_dict,
                                                                 cpgs_per_chr_cumsum=cpgs_per_chr_cumsum, windowed_cpg_sites_dict=windowed_cpg_sites_dict, windowed_cpg_sites_dict_reverse=windowed_cpg_sites_dict_reverse,
                                                                 quality_limit=quality_limit, verbose=verbose, debug=debug)  # TODO: simplify these inputs!

        assert len(methylation_data_coo.toarray()[0]) == total_cpg_sites

        if verbose:
            print(f"\nWriting methylation data to: {output_file}")

        # Save the matrix, which is an ndarray of shape (n_reads, n_cpgs), to a file
        scipy.sparse.save_npz(output_file, methylation_data_coo, compressed=True)

        # Report performance time
        if verbose:
            print(f"\nTime for this bam: {time.time() - time_bam:.2f} seconds")
            print(f"\nTotal time elapsed: {time.time() - time_start:.2f} seconds")

    if verbose:
        print("\nRun complete.")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##
# A command-line Python script to extract methylation data from an aligned .bam file.
# This assumes alignment with Biscuit, or with gem3-mapper and gemBS, though it should work with any aligner.
# TODO: Test with other aligners.
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
from scipy.sparse import coo_matrix
from tqdm import tqdm

# We use SeqIO to parse the reference .fa file
from Bio import SeqIO

def get_cpg_sites_from_fasta(reference_fasta, chromosomes, verbose=False, skip_cache=False):
    """
    Generate a dict of *all* CpG sites across each chromosome in the reference genome.

    This is a dict of lists, where the key is the chromosome: e.g. "chr1"
    The value is a list of CpG sites: e.g. [0, 35, 190, 212, 1055, ...]

    We store this as a dict because it's easier to portably serialize to disk as JSON.

    Args:
        reference_fasta (str): Path to the reference genome .fa file.
        chromosomes (list): List of chromosomes to include.
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
        with open(cached_cpg_sites_json, "r", encoding='utf-8') as f:
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
    for seqrecord in tqdm(SeqIO.parse(reference_fasta, "fasta"), total=len(chromosomes), disable=not verbose):
        if seqrecord.id not in chromosomes:
            tqdm.write(f"\tSkipping chromosome {seqrecord.id}")
            continue
        sequence = seqrecord.seq

        # Find all CpG sites
        # The i+1 is because we want to store the 1-based position, because .bed is wild and arguably 1-based maybe:
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
    with open(cached_cpg_sites_json, "w", encoding='utf-8') as f:
        json.dump(cpg_sites_dict, f)

    return cpg_sites_dict


def get_windowed_cpg_sites(reference_fasta, cpg_sites_dict, window_size, verbose=False, skip_cache=False):
    """
    Generate a dict of CpG sites for each chromosome in the reference genome.

    This is a dict of lists, where but each list contains a tuple of CpG ranges witin a window
    The key is the chromosome: e.g. "chr1"
    The value is a list of CpG sites: e.g. [(0, 35), (190, 212), (1055, ...)]

    We store this as a dict because it's easier to portably serialize to disk as JSON.

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
    windowed_cpg_sites_cache = os.path.splitext(reference_fasta)[0] + '.cpg_windowed_sites.json'
    windowed_cpg_sites_reverse_cache = os.path.splitext(reference_fasta)[0] + '.cpg_windowed_sites_reverse.json'

    if os.path.exists(windowed_cpg_sites_cache) and os.path.exists(windowed_cpg_sites_reverse_cache) and not skip_cache:
        if verbose:
            print(f'\tLoading windowed CpG sites from caches:\n\t\t{windowed_cpg_sites_cache}\n\t\t{windowed_cpg_sites_reverse_cache}')
        with open(windowed_cpg_sites_cache, "r", encoding='utf-8') as f:
            windowed_cpg_sites_dict = json.load(f)
        with open(windowed_cpg_sites_reverse_cache, "r", encoding='utf-8') as f:
            # This wild object_hook is to convert the keys back to integers, since JSON only supports strings as keys
            windowed_cpg_sites_dict_reverse = json.load(f, object_hook=lambda d: {int(k) if k.isdigit() else k: v for k, v in d.items()})

        return windowed_cpg_sites_dict, windowed_cpg_sites_dict_reverse

    if verbose:
        print(f'\tNo cache found (or --skip-cache=True) at: {windowed_cpg_sites_cache} or {windowed_cpg_sites_reverse_cache}')

    # Let's generate ranges (windows) of all CpGs within READ_SIZE of each other
    # This is to subsequently query our .bam/.sam files for reads containing CpGs
    # TODO: Shrink read size a little to account for trimming?
    assert 10 < window_size < 500, "Read size is atypical, please double check (only validated for ~150 bp.)"

    if verbose:
        print(f'\n\tGenerating windowed CpG sites dict (window size = {window_size} bp.)\n')

    # This is a dict of lists, where but each list contains a tuple of CpG ranges witin a window
    # The key is the chromosome: e.g. "chr1"
    # The value is a list of tuples: e.g. [(0,35), (190,212), (1055,)]
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
            if i+1 == cpg_sites_len or (cpg_sites[i+1] - cpg_sites[i]) > window_size:
                # We have a complete window
                window_end = cpg_sites[i]
                windowed_cpg_sites_dict[chrom].append((window_start, window_end))
                windowed_cpg_sites_dict_reverse[chrom][window_start] = temp_per_window_cpg_sites
                temp_per_window_cpg_sites = []
                window_start = None
                window_end = None

    # Save these to .json caches
    if verbose:
        print(f'\tSaving windowed CpG sites to caches:\n\t\t{windowed_cpg_sites_cache}\n\t\t{windowed_cpg_sites_reverse_cache}')

    with open(windowed_cpg_sites_cache, "w", encoding='utf-8') as f:
        json.dump(windowed_cpg_sites_dict, f)
    with open(windowed_cpg_sites_reverse_cache, "w", encoding='utf-8') as f:
        json.dump(windowed_cpg_sites_dict_reverse, f)

    return windowed_cpg_sites_dict, windowed_cpg_sites_dict_reverse

# TODO: Object orient this input / simplify the input?
def embedding_to_genomic_position(chromosomes, total_cpg_sites, cpg_sites_dict, cpgs_per_chr_cumsum, embedding_pos):
    """
    Given an embedding position, return the chromosome and position.

    Parameters
    ----------
    embedding_pos : int
        The embedding position.

    Returns
    -------
    chr : str
        The chromosome.
    pos : int
        The position.
    """
    assert embedding_pos >= 0 and embedding_pos < total_cpg_sites

    # Find the index of the first element in cpgs_per_chr_cumsum that is greater than or equal to the embedding position
    chr_index = np.searchsorted(cpgs_per_chr_cumsum, embedding_pos+1, side="left")

    # Now we know the chromosome, but we need to find the position within the chromosome
    # If this is the first chromosome, the position is just the embedding position
    if chr_index == 0:
        return chromosomes[chr_index], cpg_sites_dict[chromosomes[chr_index]][embedding_pos]
    # Otherwise, subtract the length of the previous chromosomes from the embedding position
    embedding_pos -= cpgs_per_chr_cumsum[chr_index-1]
    return chromosomes[chr_index], cpg_sites_dict[chromosomes[chr_index]][embedding_pos]


# TODO: Object orient this input / simplify the input?
def genomic_position_to_embedding(chromosomes, cpg_sites_dict, cpgs_per_chr_cumsum, chrom, pos):
    """
    Given a genomic position, return the embedding position.

    Parameters
    ----------
    chrom : str
        The chromosome.
    pos : int
        The position.

    Returns
    -------
    embedding_pos : int
        The embedding position.
    """
    assert chrom in chromosomes
    
    # Find the index of the chromosome
    chr_index = chromosomes.index(chrom)
    # Find the index of the CpG site in the chromosome
    cpg_index = cpg_sites_dict[chrom].index(pos)
    # If this is the first chromosome, the embedding position is just the CpG index
    if chr_index == 0:
        return cpg_index
    # Otherwise, add the length of the previous chromosomes to the CpG index
    return cpg_index + cpgs_per_chr_cumsum[chr_index-1]


def extract_methylation_data_from_bam(input_bam, total_cpg_sites,
                                      chromosomes, cpg_sites_dict, cpgs_per_chr_cumsum,
                                      windowed_cpg_sites_dict, windowed_cpg_sites_dict_reverse,
                                      quality_limit=0, verbose=False, paranoid=False):
    """
    Extract methylation data from a .bam file.
    
    Args:

    """
    input_bam_object = pysam.AlignmentFile(input_bam, "rb", require_index=True, threads=1)
    
    if verbose:
        print(f'\tTotal reads: {input_bam_object.mapped}')

    # Temporary storage for loading into a COO matrix
    # For a COO, we want three lists:
    # row: the read number (we'll store a read name -> ID dict perhaps?)
    # column: cpg #
    # data: methylation state (1 [methylated], -1 [unmethylated], and 0 [no data])
    next_coo_row_number = 0
    read_name_to_row_number = {}

    # TODO: Modularize/clean?
    coo_row = []
    coo_col = []
    coo_data = []

    DEBUG = False

    # This is slow, but we only run it once and store the results for later
    for chrom, windowed_cpgs in tqdm(windowed_cpg_sites_dict.items(), disable=not verbose):
        for start_pos, stop_pos in windowed_cpgs:
            cpgs_within_window = windowed_cpg_sites_dict_reverse[chrom][start_pos]
            if DEBUG: print(f"\tPosition: {start_pos} - {stop_pos}")
            if DEBUG: print(f"\tCovered CpGs: {cpgs_within_window}")
            # cpgs_within_window = [10469]
            # This returns an AlignmentSegment object from PySam
            for aligned_segment in input_bam_object.fetch(contig=chrom, start=start_pos, end=stop_pos):
                if aligned_segment.mapping_quality < quality_limit:
                    # if DEBUG: print(f"Skipping read {aligned_segment.query_name} with MAPQ {aligned_segment.mapping_quality} < {QUALITY_THRESHOLD}")
                    continue

                if aligned_segment.is_duplicate:
                    continue
                if aligned_segment.is_qcfail:
                    continue
                if aligned_segment.is_secondary:
                    continue

                if DEBUG: print(aligned_segment)

                if paranoid:
                    # Validity tests
                    assert aligned_segment.is_mapped
                    assert aligned_segment.is_supplementary is False
                    assert aligned_segment.reference_start <= stop_pos
                    assert aligned_segment.reference_end >= start_pos
                    assert aligned_segment.query_alignment_sequence == aligned_segment.get_forward_sequence()

                    # Ensure alignment methylation tags exist
                    assert aligned_segment.has_tag("XB") # GEM3/Blueprint may add this (also in Biscuit, but different!)
                    assert aligned_segment.has_tag("MD") # Location of mismatches (methylation)
                    assert aligned_segment.has_tag("YD") # Bisulfite conversion strand label (f: OT/CTOT C->T or r: OB/CTOB G->A)

                # TODO: We ignore paired/unpaired read status for now, and treat them the same
                # Should we treat paired reads / overlapping reads differently?

                bisulfite_parent_strand_is_reverse = None
                if aligned_segment.has_tag("YD"): # Biscuit tag
                    yd_tag = aligned_segment.get_tag("YD")
                    if yd_tag == "f": # Forward = C→T
                        # This read derives from OT/CTOT strand: C->T substitutions matter (C = methylated, T = unmethylated),
                        bisulfite_parent_strand_is_reverse = False
                    elif yd_tag == "r": # Reverse = G→A
                        # This read derives from the OB/CTOB strand: G->A substitutions matter (G = methylated, A = unmethylated)
                        bisulfite_parent_strand_is_reverse = True
                elif aligned_segment.has_tag("XB"): # gem3 / blueprint tag
                    raise NotImplementedError("XB tag not validated yet")
                    xb_tag = aligned_segment.get_tag("XB") # XB:C = Forward / Reference was CG
                    # TODO: Double check one or two of these gem3 tags manually.
                    if xb_tag == "C":
                        bisulfite_parent_strand_is_reverse = False
                    elif xb_tag == "G": # XB:G = Reverse / Reference was GA
                        bisulfite_parent_strand_is_reverse = True

                # TODO: Think about this old note of mine -- is this correct?
                # We have paired-end reads; one half should (the "parent strand") has the methylation data.
                # The other half (the "daughter strand") was the complement created by PCR, which we don't care about.
                if bisulfite_parent_strand_is_reverse != aligned_segment.is_reverse:
                    # Skip if we're not on the bisulfite-converted parent strand.
                    if DEBUG: print(f"\t\t Not on methylated strand, skipping: {aligned_segment.query_name}")
                    continue

                # get_aligned_pairs returns a list of tuples of (read_pos, ref_pos)
                # We filter this to only include the specific CpG sites from above
                this_segment_cpgs = [e for e in aligned_segment.get_aligned_pairs(matches_only=True) if e[1]+1 in cpgs_within_window]

                # Ok we're on the same strand as the methylation (right?)
                # Let's compare the possible CpGs in this interval to the reference and note status
                #   A methylated C will be *unchanged* and read as C (pair G)
                #   An unmethylated C will be *changed* and read as T (pair A)
                if DEBUG: print(f"bisulfite_parent_strand_is_reverse: {bisulfite_parent_strand_is_reverse}")
                for query_pos, ref_pos in this_segment_cpgs:
                    query_base = aligned_segment.query_sequence[query_pos]
                    query_base2 = aligned_segment.get_forward_sequence()[query_pos]
                    query_base3 = aligned_segment.query_alignment_sequence[query_pos]
                    assert query_base == query_base2
                    assert query_base == query_base3
                    #ref_base = aligned_segment.get_reference_sequence()[ref_pos - aligned_segment.reference_start]
                    #assert ref_base.upper() == "C", f"Expected C at {ref_pos} but got {ref_base}"
                    # print(f"\t{query_pos} {ref_pos} {query_base} {query_base2} {query_base3} {ref_base}")
                    # print(f"\t{query_pos} {ref_pos} C->{query_base}")
                    # If query == reference, we're methylated.
                    # If query !== reference, we're unmethylated.
                    # NOTE: there's an edge where if query != reference AND query == random base, we're an SNV

                    # Store in our sparse array
                    coo_row.append(next_coo_row_number)
                    # instead of ref_pos, we want the index of the CpG in the list of CpGs
                    # coo_col.append(ref_pos)

                    # TODO: Object orient these inputs? -- lots of bad inheritence style here
                    coo_col.append(genomic_position_to_embedding(chromosomes, cpg_sites_dict, cpgs_per_chr_cumsum, chrom, ref_pos+1))
                    #coo_data.append(1 if query_base == "C" else 0)
                    if query_base == "C":
                        # Methylated
                        coo_data.append(1)
                        if DEBUG: print(f"\t\tMethylated")
                    elif query_base == "T":
                        coo_data.append(0)
                        # Unmethylated
                        if DEBUG: print(f"\t\tUnmethylated")
                    else:
                        coo_data.append(-1) # or just 0?
                        if DEBUG: print(f"\t\tUnknown!")
                        # raise ValueError(f"Unexpected query base {query_base} at {query_pos} {ref_pos}")

                read_name_to_row_number[aligned_segment.query_name] = next_coo_row_number
                next_coo_row_number += 1
                if DEBUG: print("************************************************\n")

                #query_bp = aligned_segment.query_sequence[pileupread.query_position]
                #reference_bp = aligned_segment.get_reference_sequence()[aligned_segment.reference_start - pileupcolumn.reference_pos].upper()

    return coo_matrix((coo_data, (coo_row, coo_col)), shape=(len(read_name_to_row_number), total_cpg_sites)).toarray()

def main():
    """
    Our main function. Call argparse and whatnot.

    Options are:
    --input-bam: Input .bam file to process.
    --reference-fasta: Reference genome fasta file (critical to determine CpG sites).
    --verbose: Verbose output.
    --paranoid: Paranoid mode (extensive validity checking).
    --output-file: Output file. Default is to use the input file name with a .methylation.npy extension.
    --overwrite: Overwrite output file if it exists.

    Returns: Nothing.
    """

    # Let's do our argparse stuff
    parser = argparse.ArgumentParser(description='Extract read-level methylation data from an aligned .bam file and store as a numpy array.')
    parser.add_argument('--input-bam', help='Input .bam file to process.', required=True)
    parser.add_argument('--reference-fasta', help='Reference genome fasta file (critical to determine CpG sites).', required=True)
    parser.add_argument('--quality-limit', help='Quality filter for aligned reads (default = 20)', default=20, type=int)
    parser.add_argument('--window-size', help='Window size (default = 150)', default=150, type=int)
    parser.add_argument('--verbose', help='Verbose output.', action='store_true')
    parser.add_argument('--skip-cache', help='De-novo generate CpG sites (slow).', action='store_true')
    parser.add_argument('--paranoid', help='Paranoid mode (extensive validity checking).', action='store_true')
    parser.add_argument('--output-file', help='Output file. Default is to use the input file name with a .methylation.npy extension.', default=None)
    parser.add_argument('--overwrite', help='Overwrite output file if it exists.', action='store_true')

    # Parse the arguments & set up some variables
    args = parser.parse_args()
    input_bam = args.input_bam
    reference_fasta = args.reference_fasta
    verbose = args.verbose
    quality_limit = args.quality_limit
    skip_cache = args.skip_cache
    paranoid = args.paranoid
    output_file = args.output_file
    window_size = args.window_size

    if args.output_file is None:
        output_file = os.path.splitext(input_bam)[0] + '.methylation.npy'

    time_start = time.time()
    # Print run information
    if verbose:
        print(f'Input bam: {input_bam}')
        print(f'Reference fasta: {reference_fasta}')
        print(f'Output file: {output_file}')

    # Check if the input file exists & readable
    assert os.access(input_bam, os.R_OK), f"Input file is not readable: {input_bam}"

    # Check if the output file exists
    if os.path.exists(output_file):
        if args.overwrite:
            if verbose:
                print('\tOutput file exists and --overwrite specified. Overwriting.')
            assert os.access(output_file, os.W_OK), f"Output file is not writable: {output_file}"
        else:
            print('\tOutput file exists and --overwrite not specified. Exiting.')
            sys.exit(1)

    chromosomes = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]

    # We need to obtain all cpg sites in the reference genome
    # NOTE: This can take a while (~10 minutes for GRCh38 if not cached)
    cpg_sites_dict = get_cpg_sites_from_fasta(reference_fasta=reference_fasta,
                                              chromosomes=chromosomes,
                                              verbose=verbose,
                                              skip_cache=skip_cache)

    # How many CpG sites are there?
    total_cpg_sites = sum([len(v) for v in cpg_sites_dict.values()])
    if verbose:
        print(f"Total CpG sites: {total_cpg_sites:,}")

    assert(total_cpg_sites>28_000_000) # Validity check for hg38

    # Count the number of CpGs per chromosome
    cpgs_per_chr = {k: len(v) for k, v in cpg_sites_dict.items()}

    # Add up the number of CpGs per chromosome, e.g. chr1, then chr1+chr2, then chr1+chr2+chr3, etc
    cpgs_per_chr_cumsum = np.cumsum([cpgs_per_chr[k] for k in chromosomes])


    # TODO: Move these into to a formal test framework
    # TODO: Simplify the input framework (object orient the window / cpg dict?)
    ### Tests
    assert cpgs_per_chr_cumsum[-1] == total_cpg_sites
    assert embedding_to_genomic_position(chromosomes, total_cpg_sites, cpg_sites_dict, cpgs_per_chr_cumsum, 0) == ("chr1", cpg_sites_dict["chr1"][0])
    assert embedding_to_genomic_position(chromosomes, total_cpg_sites, cpg_sites_dict, cpgs_per_chr_cumsum, 1) == ("chr1", cpg_sites_dict["chr1"][1])
    # Edges
    assert embedding_to_genomic_position(chromosomes, total_cpg_sites, cpg_sites_dict, cpgs_per_chr_cumsum, cpgs_per_chr_cumsum[0]) == ("chr2", cpg_sites_dict["chr2"][0])
    assert embedding_to_genomic_position(chromosomes, total_cpg_sites, cpg_sites_dict, cpgs_per_chr_cumsum, cpgs_per_chr_cumsum[-1]-1) == ("chrY", cpg_sites_dict["chrY"][-1])
    ### Tests
    assert genomic_position_to_embedding(chromosomes, cpg_sites_dict, cpgs_per_chr_cumsum, "chr1", cpg_sites_dict["chr1"][0]) == 0
    assert genomic_position_to_embedding(chromosomes, cpg_sites_dict, cpgs_per_chr_cumsum, "chr1", cpg_sites_dict["chr1"][1]) == 1
    # Edges
    assert genomic_position_to_embedding(chromosomes, cpg_sites_dict, cpgs_per_chr_cumsum, "chr2", cpg_sites_dict["chr2"][0]) == cpgs_per_chr_cumsum[0]
    assert genomic_position_to_embedding(chromosomes, cpg_sites_dict, cpgs_per_chr_cumsum, "chrY", cpg_sites_dict["chrY"][-1]) == cpgs_per_chr_cumsum[-1]-1
    ########

    # Get our windowed_cpg_sites, hopefully cached!
    if verbose:
        print('\nLoading (or generating) windowed CpG sites for reference genome.')

    windowed_cpg_sites_dict, windowed_cpg_sites_dict_reverse = get_windowed_cpg_sites(reference_fasta=reference_fasta,
                                                                                      cpg_sites_dict=cpg_sites_dict,
                                                                                      window_size=window_size,
                                                                                      verbose=verbose,
                                                                                      skip_cache=skip_cache)
    

    # Extract methylation data as a COO sparse matrix
    if verbose:
        print(f'\nExtracting methylation data from: {input_bam}')

    methylation_data_coo = extract_methylation_data_from_bam(input_bam=input_bam,
                                                             total_cpg_sites=total_cpg_sites,
                                                         chromosomes=chromosomes, cpg_sites_dict=cpg_sites_dict, cpgs_per_chr_cumsum=cpgs_per_chr_cumsum, # TODO: simplify these inputs
                                                         windowed_cpg_sites_dict=windowed_cpg_sites_dict,
                                                         windowed_cpg_sites_dict_reverse=windowed_cpg_sites_dict_reverse,
                                                         quality_limit=quality_limit,
                                                         verbose=verbose,
                                                         paranoid=paranoid)

    assert(len(methylation_data_coo[0]) == total_cpg_sites)
    
    if verbose:
        print(f'\nWriting methylation data to: {output_file}')

    # Save the matrix, which is an ndarray of shape (n_reads, n_cpgs), to a file
    np.save(output_file, methylation_data_coo, allow_pickle=True)

    # Report performance time
    if verbose:
        time_end = time.time()
        print(f'\nTime elapsed: {time_end - time_start:.2f} seconds')

if __name__ == '__main__':
    main()
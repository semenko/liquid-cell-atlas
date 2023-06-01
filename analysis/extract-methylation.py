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
#import numpy as np
#import pysam
# from scipy.sparse import coo_matrix
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
    cached_cpg_sites_json = os.path.splitext(reference_fasta)[0] + "cpg_all_sites.json"

    if verbose:
        print(f"\tLoading all CpG sites for: {reference_fasta}")

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


def get_windowed_cpg_sites(reference_fasta, window_size, chromosomes, verbose=False, skip_cache=False):
    """
    Generate a dict of CpG sites for each chromosome in the reference genome.

    This is a dict of lists, where but each list contains a tuple of CpG ranges witin a window
    The key is the chromosome: e.g. "chr1"
    The value is a list of CpG sites: e.g. [(0, 35), (190, 212), (1055, ...)]

    We store this as a dict because it's easier to portably serialize to disk as JSON.

    Args:
        reference_fasta (str): Path to the reference genome .fa file.
        window_size (int): Size of the window to use.
        chromosomes (list): List of chromosomes to include.
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
            windowed_cpg_sites_dict_reverse = json.load(f)

        return windowed_cpg_sites_dict, windowed_cpg_sites_dict_reverse

    if verbose:
        print(f'\tNo cache found (or --skip-cache=True) at: {windowed_cpg_sites_cache} or {windowed_cpg_sites_reverse_cache}')

    # Let's generate ranges (windows) of all CpGs within READ_SIZE of each other
    # This is to subsequently query our .bam/.sam files for reads containing CpGs
    # TODO: Shrink read size a little to account for trimming?
    assert 10 < window_size < 500, "Read size is atypical, please double check (only validated for ~150 bp.)"

    # We need to obtain all cpg sites in the reference genome
    # NOTE: This can take a while (~10 minutes for GRCh38 if not cached)
    cpg_sites_dict = get_cpg_sites_from_fasta(reference_fasta=reference_fasta,
                                              chromosomes=chromosomes,
                                              verbose=verbose,
                                              skip_cache=skip_cache)

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


def extract_methylation_data_from_bam(input_bam, windowed_cpg_sites_dict, windowed_cpg_sites_dict_reverse, quality_limit=0, verbose=False):

    return False

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
    parser.add_argument('--output-file', help='Output file. Default is to use the input file name with a .methylation.npy extension.', default=None)
    parser.add_argument('--overwrite', help='Overwrite output file if it exists.', action='store_true')

    # Parse the arguments & set up some variables
    args = parser.parse_args()
    input_bam = args.input_bam
    reference_fasta = args.reference_fasta
    verbose = args.verbose
    quality_limit = args.quality_limit
    skip_cache = args.skip_cache
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

    # Get our windowed_cpg_sites, hopefully cached!
    if verbose:
        print('\nLoading (or generating) windowed CpG sites for reference genome.')

    windowed_cpg_sites_dict, windowed_cpg_sites_dict_reverse = get_windowed_cpg_sites(reference_fasta=reference_fasta,
                                                                                      window_size=window_size,
                                                                                      chromosomes=chromosomes,
                                                                                      verbose=verbose,
                                                                                      skip_cache=skip_cache)
    

    # Extract methylation data
    if verbose:
        print(f'\nExtracting methylation data from: {input_bam}')

    methylation_data = extract_methylation_data_from_bam(input_bam=input_bam,
                                                         windowed_cpg_sites_dict=windowed_cpg_sites_dict,
                                                         windowed_cpg_sites_dict_reverse=windowed_cpg_sites_dict_reverse,
                                                         quality_limit=quality_limit,
                                                         verbose=verbose)
    
    # Write methylation data to a COO sparse matrix
    if verbose:
        print(f'\nWriting methylation data to: {output_file}')

    #write_methylation_data(methylation_data=methylation_data,
    #                       output_file=output_file,
    #                       verbose=verbose,
    #                       paranoid=paranoid)
    

    # Report performance time
    if verbose:
        time_end = time.time()
        print(f'\nTime elapsed: {time_end - time_start:.2f} seconds')

if __name__ == '__main__':
    main()
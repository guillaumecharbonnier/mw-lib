#!/usr/bin/env python

import sys
import os
import pybedtools
from bx.bbi.bigwig_file import BigWigFile
from gtftoolkit.utils import intervals, make_tmp_file
from tempfile import NamedTemporaryFile
import multiprocessing
import pandas as pd
from itertools import repeat
import numpy as np
import argparse


op = argparse.ArgumentParser(description="Create a count matrix from a set of "
                            "BigWig files. The resolution is set through the "
                            "'window_size' argument.", prog=sys.argv[0])
op.add_argument('-g', 
                '--genome',
                type=str,
                help='Genome (e.g mm9)', required=True)
op.add_argument('-w', '--window_size',
                type=int, 
                default=3000,
                help='The size of the windows where the mean of counts will be computed.')
op.add_argument('-t', '--threads',
                type=int, 
                default=1,
                help='Number of threads.')
op.add_argument('-b', '--bw_list', nargs='+', 
                help='A list of bigwig file separated by spaces.')
op.add_argument('-o', '--out_file',
                type=argparse.FileType('w'),
                help='Output file name.')  
op.add_argument('-v', '--verbose',
                action='store_true',
                help='Verbose mode.')  
op.add_argument('-m', '--min_row_sum',
                type=float,
                default=1,
                help="A row (region) should sum to more than 'min_row_sum' over "
                "the bigWigs to be included in the final matrix.")  
options = op.parse_args()



def big_wig_summary_worker((span,
                            bw_list,
                            region_bed_file_name,
                            nb_proc,
                            verbose
                            )):


    results = list()
    bw_label = [os.path.basename(p) for p in bw_list]
    bw_label = [os.path.splitext(os.path.basename(p))[0] for p in bw_list]
    
    if verbose:
        sys.stderr.write("Processing: " + region_bed_file_name)
    for big_wig, cpt in zip(bw_list, range(len(bw_list))):

        bigwig = BigWigFile(open(big_wig, "r"))
        if verbose:

            sys.stderr.write(
                "Computing coverage for file: " +
                big_wig +
                " [" +
                str(multiprocessing.current_process()) +
                "], " +
                str(span[1] - span[0]) + " chunks to process.\n")

        bed_windows = pybedtools.BedTool(region_bed_file_name)

        chr_cur = None
        
        # Loop through bed lines (features object)

        for i in bed_windows[slice(span[0], span[1])]:

            if chr_cur == None:
                chr_cur = i.chrom

            else:
                if i.chrom != chr_cur:
                    chr_cur = i.chrom

            # Note: bigWig is zero-based/half open as bed.
            bw_sum = bigwig.query(i.chrom,
                                  i.start,
                                  i.end,
                                  1)


            if bw_sum is not None:
                bw_sum = bw_sum[0]['mean']
                bw_sum = np.nan_to_num(bw_sum)
                bw_sum = np.round(bw_sum, 2)
            else:
                bw_sum=0.00

            results.append((i.chrom + ":" + str(i.start),
                            bw_label[cpt], 
                            float(bw_sum)))

            
    if verbose:

        sys.stderr.write(
            "Computing coverage for file: " +
            big_wig +
            " [" +
            str(multiprocessing.current_process()) +
            "]. Job done.\n")


    return results


def bw_coverage(
        out_file=None,
        bw_list=None,
        nb_proc=None,
        genome=None,
        verbose=False,
        min_row_sum=None,
        window_size=None):
    """
    Compute transcript coverage with one or several bigWig.
    -------------------------------------------------------
    Uses bx-python as interface to kent utilities.

    """
    
    window = pybedtools.BedTool()

    # Produce a set of windowed genomic regions 
    tmp_file = NamedTemporaryFile(delete=False)
    a = window.window_maker(genome=genome, w=window_size)
    a.saveas(tmp_file.name)

    tokens = intervals(range(len(pybedtools.BedTool(tmp_file.name))), nb_proc)

    pool = multiprocessing.Pool(nb_proc)
    coverage_list = pool.map_async(big_wig_summary_worker,
                                   zip(tokens,
                                       repeat(bw_list),
                                       repeat(tmp_file.name),
                                       repeat(nb_proc),
                                       repeat(verbose))).get(9999999)

    # Unlist the list of list
    coverage_list = [item for sublist in coverage_list for item in sublist]

    # Prepare a data.frame to collect the results
    # Loop through the bigWig file to compute coverage
    dataframe = pd.DataFrame(columns=None)

    if verbose:
        sys.stderr.write("Retrieving results.\n")

    for i in coverage_list:
        dataframe.ix[i[0], i[1]] = float(i[2])

    dataframe =dataframe.loc[(dataframe.sum(axis=1) > min_row_sum), ]
    dataframe.to_csv(out_file, sep="\t")


bw_coverage(
        out_file=options.out_file.name,
        bw_list=options.bw_list,
        nb_proc=options.threads,
        genome=options.genome,
        window_size=options.window_size,
        min_row_sum=options.min_row_sum,
        verbose=options.verbose)
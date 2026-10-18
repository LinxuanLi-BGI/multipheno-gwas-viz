#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse

import pandas as pd
import numpy as np

def reformat_chrom(chr_col:pd.Series)->pd.Series:
    """
    convert 'chrom' columns of SummaryStats object to int type.\n
    chrX and chrY is numbered after autosome,
    invaild value is converted to -1.

    Parameters:
    ----------
        chr_col : Series
            Series or array like, the data to convert.
    
    Returns:
    ----------
        Series
            converted chr columns
    """
    def reformat(x):
        try:
            # int type
            if isinstance(x, (int, np.integer, float, np.floating)):
                return x
            # string type
            elif type(x) is str:
                if x[0:3] == 'chr':
                    x = x[3:]
                if x == "X":
                    return 23
                if x == "Y":
                    return 24
                if x == 'M':
                    return 25
                else:
                    return int(x)
            # not supported type
            else:
                return -1
        except ValueError:
            return -1
    return np.frompyfunc(reformat,1,1)(chr_col).astype(int)

# get mapped gene from 
def mapped(chrom:int, pos:int, all_gene=False): 
    Chr_b= known_gene['#hg38.knownCanonical.chrom']==chrom
    this_chrom =  known_gene[Chr_b]
    this_chrom = this_chrom[this_chrom['hg38.kgXref.geneSymbol']!="Y_RNA"]
    Bg_b = this_chrom['hg38.knownCanonical.chromStart']<=pos
    Ed_b = this_chrom['hg38.knownCanonical.chromEnd']>=pos

    matched:pd.DataFrame = this_chrom[Bg_b&Ed_b].copy()
    if len(matched) == 0:
        return np.nan
    elif len(matched) == 1:
        return matched['hg38.kgXref.geneSymbol'].iloc[0]
    else:
        matched['center_distance'] = (
            (matched['hg38.knownCanonical.chromStart'] + matched['hg38.knownCanonical.chromEnd']) / 2 - pos
        ).abs()
        matched = matched.sort_values('center_distance')
        if all_gene:
            return '; '.join(matched['hg38.kgXref.geneSymbol'].astype(str).drop_duplicates())
        return matched['hg38.kgXref.geneSymbol'].iloc[0]

uf_mapped = np.frompyfunc(mapped, 3, 1)

def main(args):
    parser = argparse.ArgumentParser(description='Annotate SNPs with overlapping genes.')
    parser.add_argument('input_file')
    parser.add_argument('known_gene_file')
    parser.add_argument('output_dir')
    parser.add_argument('--all-gene', action='store_true',
                        help='Retain all overlapping genes, separated by "; ". '
                             'By default, select the overlapping gene with the nearest midpoint.')
    options = parser.parse_args(args[1:])

    global known_gene
    known_gene = pd.read_csv(options.known_gene_file, sep='\t')
    known_gene['#hg38.knownCanonical.chrom'] = reformat_chrom(known_gene['#hg38.knownCanonical.chrom'])

    ############################################################
    # 1. Load input
    ############################################################

    all_sig_path = options.input_file # ./output/step2.all_sig.csv
    output_path = options.output_dir # ./output/

    all_sig_df = pd.read_csv(all_sig_path, sep='\t') # chrom pos A1 A2 freq size beta se pvalue log10p

    ############################################################
    # 2. Extract unique SNP and map to gene label
    ############################################################

    all_sig_df['id'] = all_sig_df['chrom'].astype(str) + ':' + all_sig_df['pos'].astype(str)

    all_snp = pd.concat([
        all_sig_df.groupby('id')['chrom'].min(),
        all_sig_df.groupby('id')['pos'].min(),
        all_sig_df.groupby('id')['phenotype'].nunique().rename('N'),
        all_sig_df.groupby('id')['phenotype'].unique().map(lambda x:'; '.join(x)).rename('phenotypes')
    ],axis=1)

    all_snp['gene'] = uf_mapped(all_snp['chrom'], all_snp['pos'], options.all_gene)
    all_snp.to_csv(os.path.join(output_path, 'step3.all_unique_snps.csv'), sep='\t', index=False)

    ############################################################
    # 3. map SNP level gene annotation to sig df
    ############################################################

    all_sig_df['gene'] = all_sig_df['id'].map(all_snp['gene'])
    all_sig_df.to_csv(os.path.join(output_path, 'step3.all_sig_gene.csv'), sep='\t', index=False)

if __name__ == "__main__":
    main(sys.argv)

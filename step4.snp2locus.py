#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import pandas as pd
import numpy as np

def get_locusID(significants:pd.DataFrame, window:int = 500_000):

    locusID = 0
    last = None
    for key,this in significants.iterrows():
        if last is None:
            pass
        else:
            if last['chrom']!=this['chrom']:
                locusID += 1
            elif this['pos'] - last['pos'] > window:
                locusID += 1
            else:
                pass
        significants.loc[key,'locusID'] = locusID
        last = this
    significants['locusID'] = significants['locusID'].astype(int)
    return significants

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
def mapped(chrom:int, pos:int):
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
        return matched['hg38.kgXref.geneSymbol'].iloc[0]

def closed(chrom:int ,pos:int , enable = True, direction = True):
    mapped_gene = mapped(chrom, pos)
    if pd.isna(mapped_gene) and enable:
        Chr_b= known_gene['#hg38.knownCanonical.chrom']==chrom
        this_chrom = known_gene[Chr_b]
        this_chrom = this_chrom[this_chrom['hg38.kgXref.geneSymbol']!="Y_RNA"]
        if this_chrom.empty:
            return np.nan
        S_distance = (this_chrom['hg38.knownCanonical.chromStart'] - pos).abs()
        E_distance = (this_chrom['hg38.knownCanonical.chromEnd'] - pos).abs()

        tplt = "{}"
        if S_distance.min() < E_distance.min():
            K_distance:pd.Series = S_distance
            if direction:
                tplt = "<-{}"
        else:
            K_distance:pd.Series = E_distance
            if direction:
                tplt = "{}->"

        closest = K_distance.idxmin()
        return tplt.format(this_chrom.loc[closest,'hg38.kgXref.geneSymbol'])
    else:
        return mapped_gene
uf_closed = np.frompyfunc(closed, 4, 1)

def main(args):
    parser = argparse.ArgumentParser(description='Group SNPs into loci and annotate lead SNPs.')
    parser.add_argument('input_file')
    parser.add_argument('known_gene_file')
    parser.add_argument('locus_window', type=int)
    parser.add_argument('output_dir')
    parser.add_argument('--no-nearest-gene', action='store_true',
                        help='Leave intergenic lead SNPs unannotated instead of assigning '
                             'the gene with the nearest boundary.')
    parser.add_argument('--no-gene-direction', action='store_true',
                        help='Omit <- and -> arrows on nearest-gene labels for intergenic lead SNPs.')
    options = parser.parse_args(args[1:])

    global known_gene
    known_gene = pd.read_csv(options.known_gene_file, sep='\t')
    known_gene['#hg38.knownCanonical.chrom'] = reformat_chrom(known_gene['#hg38.knownCanonical.chrom'])

    ############################################################
    # 1. Load input
    ############################################################
    
    all_sig_path = options.input_file # ./output/step3.all_sig_gene.csv

    locus_window = options.locus_window

    output_path = options.output_dir

    all_sig_df = pd.read_csv(all_sig_path, sep='\t') # chrom pos phenotype category

    ############################################################
    # 2. Sort and mark locusID
    ############################################################

    all_sig_df = all_sig_df.sort_values(['phenotype','chrom','pos'])

    all_sig_phenotypes = all_sig_df['phenotype'].unique()
    all_sig_locus_list = []
    for P in all_sig_phenotypes:
        this = all_sig_df.query('phenotype == @P').copy()
        this = get_locusID(this, window=locus_window)
        all_sig_locus_list.append(this)
    all_sig_locus:pd.DataFrame = pd.concat(all_sig_locus_list, ignore_index=True)

    # UID: phenotype(0)
    all_sig_locus['UID'] = all_sig_locus['phenotype'] + '(' + all_sig_locus['locusID'].astype(str) + ')'

    all_sig_locus.to_csv(os.path.join(output_path, 'step4.all_sig_l.csv'), sep='\t', index=False)


    ############################################################
    # 3. locus level summary
    ############################################################
    Locus_UID = all_sig_locus.groupby(['phenotype','locusID'])['UID'].unique().map(lambda x:x[0])
    phenotype_cate = all_sig_locus.groupby(['phenotype','locusID'])['category'].unique().map(lambda x:x[0])
    # locus coordinates
    Locus_chrom = all_sig_locus.groupby(['phenotype','locusID'])['chrom'].min()
    Locus_start = all_sig_locus.groupby(['phenotype','locusID'])['pos'].min().rename('start')
    Locus_end = all_sig_locus.groupby(['phenotype','locusID'])['pos'].max().rename('end')

    # Flag mega-loci using their observed span, independently of the locus window.
    Locus_span = Locus_end - Locus_start
    for (phenotype, locusID), span in Locus_span[Locus_span > 500_000].items():
        locus_key = (phenotype, locusID)
        print(
            f"Warning: mega-locus detected: phenotype={phenotype}, locusID={locusID}, "
            f"chr{Locus_chrom.loc[locus_key]}:{Locus_start.loc[locus_key]}-{Locus_end.loc[locus_key]}, "
            f"span={span:,} bp (>500 kb). Consider progressively smaller locus windows "
            "(e.g., 200 kb, then 100 kb) and check whether locus assignments change.",
            file=sys.stderr,
        )

    # locus mapped gene: GeneA[4], GeneB[6] 
    Locus_N = all_sig_locus.groupby(['phenotype','locusID'])['pos'].count().rename('NSNP')
    # Count each retained overlapping gene separately in multi-gene annotations.
    locus_genes = all_sig_locus[['phenotype', 'locusID', 'gene']].copy()
    locus_genes['gene'] = locus_genes['gene'].map(
        lambda genes: genes.split('; ') if isinstance(genes, str) else genes
    )
    locus_genes = locus_genes.explode('gene')
    locus_gene_snp_count:pd.DataFrame = locus_genes.groupby(['phenotype','locusID'])['gene'].value_counts().rename('count').reset_index()
    locus_gene_snp_count['display'] = locus_gene_snp_count['gene'] + '[' + locus_gene_snp_count['count'].astype(str) + ']'
    Locus_all_gene_count = locus_gene_snp_count.groupby(['phenotype','locusID'])['display'].unique().map(', '.join).rename('all_snp_mapped_gene')

    # lead SNP
    Locus_topP = all_sig_locus.groupby(['phenotype','locusID'])['pvalue'].min()
    locus_top_logP = all_sig_locus.groupby(['phenotype','locusID'])['log10p'].max()
    Locus_toppos = pd.Series(
        data = all_sig_locus.loc[all_sig_locus.groupby(['phenotype','locusID'])['log10p'].idxmax().values,'pos'].values,
        index= Locus_N.index,
        name = 'top'
    )
    Locus_beta = pd.Series(
        data = all_sig_locus.loc[all_sig_locus.groupby(['phenotype','locusID'])['log10p'].idxmax().values, 'beta'].values,
        index= Locus_topP.index,
        name = 'beta'
    )
    Locus_A1 = pd.Series(
        data = all_sig_locus.loc[all_sig_locus.groupby(['phenotype','locusID'])['log10p'].idxmax().values, 'A1'].values,
        index= Locus_topP.index,
        name = 'A1'
    )
    # concat

    all_locus_summary:pd.DataFrame = pd.concat([
        Locus_UID,
        phenotype_cate,
        Locus_chrom,
        Locus_start,
        Locus_end, 
        Locus_N,
        Locus_all_gene_count,
        Locus_topP,
        locus_top_logP,
        Locus_toppos,
        Locus_beta,
        Locus_A1
    ] ,axis=1).reset_index()

    all_locus_summary['top_snp_nearest_gene'] = uf_closed(
        all_locus_summary['chrom'], all_locus_summary['top'],
        not options.no_nearest_gene, not options.no_gene_direction,
    )

    ############################################################
    # 4. cross-phenotype recurrence
    ############################################################

    for i, locus in all_locus_summary.iterrows():
        chrom = locus.chrom
        start = locus.start - locus_window
        end = locus.end + locus_window
        UID = locus.UID

        all_recurrence_in_window:pd.DataFrame = all_locus_summary.query('chrom == @chrom and top >= @start and top <= @end and UID != @UID')
        all_locus_summary.loc[i, 'Recurrence'] = ', '.join(all_recurrence_in_window.sort_values('log10p', ascending=False)['UID'])

    all_locus_summary.to_csv(os.path.join(output_path, 'step4.all_locus.csv'), sep='\t', index=False)

if __name__ == "__main__":
    main(sys.argv)

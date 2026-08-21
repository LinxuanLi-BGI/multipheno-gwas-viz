#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np

def format_trait_count(series: pd.Series) -> str:
    if series.empty:
        return ""
    parts = series.index + "[" + series.astype(str) + "]"
    return ", ".join(parts.tolist())

def main(args):
    ############################################################
    # 1. Load input
    ############################################################
    all_locus_path = args[1] # ./output/step5.all_locus_cata.csv

    gwas_catalog_path = args[2] # ./external/gwas_catalog_v1.0-associations_e116_r2026-07-19_full.zip

    replication_window = int(args[3]) # 500000

    output_path = args[4] # ./output/

    all_locus_summary = pd.read_csv(all_locus_path, sep='\t') # phenotype chrom start end locusID UID candidate_novel

    candidate_locus_summary = all_locus_summary.query("candidate_novel == True").copy()

    gwas_catalog = pd.read_csv(gwas_catalog_path, sep='\t', low_memory=False) # P-VALUE CHR_ID CHR_POS DISEASE/TRAIT
    # GWAS catalog also contain variants that did not reach the threshold 5e-8
    GWAS_catalog_replication_threshold = 5e-8
    GWAS_catalog_sig:pd.DataFrame = gwas_catalog.query('`P-VALUE` < @GWAS_catalog_replication_threshold').copy()
    def set_chrom(tx):
        try:
            return float(tx)
        except:
            return np.nan
    GWAS_catalog_sig['CHR_ID'] = GWAS_catalog_sig['CHR_ID'].replace({'X':23,'Y':24}).map(set_chrom)#.dropna(subset=['CHR_ID', 'CHR_POS'])
    GWAS_catalog_sig = GWAS_catalog_sig.dropna(subset=['CHR_ID', 'CHR_POS']).copy()
    GWAS_catalog_sig['CHR_ID'] = GWAS_catalog_sig['CHR_ID'].astype(int)
    GWAS_catalog_sig['CHR_POS'] = GWAS_catalog_sig['CHR_POS'].astype(int)

    ############################################################
    # 2. replication in gwas catalog
    ############################################################

    candidate_traits = candidate_locus_summary['phenotype'].unique()

    for trait in candidate_traits:
        # locus for this phenotype
        this_locus = candidate_locus_summary.query('phenotype == @trait')
    
        # for each locus
        for i, locus in this_locus.iterrows():
            start = locus.start - replication_window
            end = locus.end + replication_window
            chrom = locus.chrom

            # all reported gwas finding in this locus interval
            all_catalog_in_window = GWAS_catalog_sig.query('CHR_ID == @chrom and CHR_POS >= @start and CHR_POS <= @end')
            # stats the number of study within each label 
            study_ALL = all_catalog_in_window.groupby('DISEASE/TRAIT')['STUDY'].nunique().sort_values(ascending=False)
            candidate_locus_summary.loc[i, 'cata_ALL'] = format_trait_count(study_ALL)

    ############################################################
    # 3. save
    ############################################################

    candidate_locus_summary.to_csv(os.path.join(output_path, 'step7.candidate_locus_cata_ALL.csv'), sep='\t', index=False)

if __name__ == "__main__":
    main(sys.argv)
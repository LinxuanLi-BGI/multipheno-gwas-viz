#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np

Ancestry = {
    'EAS': [
        'East Asian',
        'Japanese',
        'Korean',
        'Chinese',
        'Taiwanese',
        'Singaporean'
    ],
    'EUR': [
        'European',
        'British',
        'Finnish',
        'Icelandic',
    ],
    'OTH': [
        'Pakistani',
        'South Asian',
        'Hispanic',
        'African',
        'Middle Eastern',
        'Latin American',
        'South East Asian',
        'Afro-Caribbean',
        'Latino'
    ]
}

def check_ancestry(A, ans='EAS'):
    target = Ancestry[ans]
    for i in target:
        if i in str(A):
            return True
    return False

# GWAS_catalog_label[0]
def format_trait_count(series: pd.Series) -> str:
    if series.empty:
        return ""
    parts = series.index + "[" + series.astype(str) + "]"
    return ", ".join(parts.tolist())

def main(args):
    ############################################################
    # 1. Load input
    ############################################################
    all_locus_path = args[1] # ./output/step4.all_locus.csv

    trait_2_catalog_path = args[2] # ./input/trait_2_catalog.csv

    gwas_catalog_path = args[3] # ./external/gwas_catalog_v1.0-associations_e116_r2026-07-19_full.zip

    replication_window = int(args[4]) # 500000

    output_path = args[5] # ./output/

    all_locus_summary = pd.read_csv(all_locus_path, sep='\t') # phenotype chrom start end locusID UID

    # phenotype catalog_labels ...
    trait_2_catalog:pd.Series = pd.read_csv(trait_2_catalog_path, sep='\t', index_col='phenotype')['catalog_labels'] 

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

    # mark ancenstry in gwas catalog
    GWAS_catalog_sig['is_EAS'] = GWAS_catalog_sig['INITIAL SAMPLE SIZE'].map(lambda x:check_ancestry(x,'EAS'))
    GWAS_catalog_sig['is_EUR'] = GWAS_catalog_sig['INITIAL SAMPLE SIZE'].map(lambda x:check_ancestry(x,'EUR'))
    GWAS_catalog_sig['is_OTH'] = GWAS_catalog_sig['INITIAL SAMPLE SIZE'].map(lambda x:check_ancestry(x,'OTH'))

    ############################################################
    # 2. replication in gwas catalog
    ############################################################

    # focus on the traits provided in trait_2_catalog table
    all_traits = trait_2_catalog.index.values

    for trait in all_traits:
        # DISEASE/TRAIT label for this phenotype
        this = trait_2_catalog[trait].split('; ')
        # locus for this phenotype
        this_locus = all_locus_summary.query('phenotype == @trait')
        # reported gwas finding in this phenotype
        this_catalog:pd.DataFrame = GWAS_catalog_sig.query('`DISEASE/TRAIT` in @this')
        # for each locus
        for i, locus in this_locus.iterrows():
            start = locus.start - replication_window
            end = locus.end + replication_window
            chrom = locus.chrom

            # select variants in locus interval
            catalog_in_window = this_catalog.query('CHR_ID == @chrom and CHR_POS >= @start and CHR_POS <= @end')
            # stats the number of unique study, in 3 ancenstry
            all_locus_summary.loc[i, 'rep_ALL'] = len(catalog_in_window['PUBMEDID'].unique())
            all_locus_summary.loc[i, 'rep_EAS'] = len(catalog_in_window.query('is_EAS')['PUBMEDID'].unique())
            all_locus_summary.loc[i, 'rep_EUR'] = len(catalog_in_window.query('is_EUR')['PUBMEDID'].unique())
            all_locus_summary.loc[i, 'rep_study'] = ', '.join(catalog_in_window['PUBMEDID'].unique().astype(str))

    ############################################################
    # 3. candidate_novel and save
    ############################################################

    # True: novel; False: known
    all_locus_summary['candidate_novel'] = (all_locus_summary['rep_ALL'] == 0)

    all_locus_summary.to_csv(os.path.join(output_path, 'step5.all_locus_cata.csv'), sep='\t', index=False)

if __name__ == "__main__":
    main(sys.argv)
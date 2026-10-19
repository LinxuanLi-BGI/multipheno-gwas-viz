#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import csv
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


def load_trait_mapping(path, catalog_traits, delimiter='; '):
    """Validate the TSV structure and check labels against the full Catalog."""
    def warn(message):
        print(f'Warning: trait_2_catalog {path!r}: {message}', file=sys.stderr)

    try:
        with open(path, encoding='utf-8-sig', newline='') as handle:
            rows = [row for row in csv.reader(handle, delimiter='\t', strict=True) if row]
        if (not delimiter or len(rows) < 2
                or not {'phenotype', 'catalog_labels'}.issubset(rows[0])
                or len(set(rows[0])) != len(rows[0])):
            raise ValueError
        header = rows[0]
        mapping = {}
        for row in rows[1:]:
            if len(row) != len(header):
                raise ValueError
            record = dict(zip(header, row))
            phenotype = record['phenotype']
            if not phenotype.strip() or phenotype in mapping:
                raise ValueError
            mapping[phenotype] = record['catalog_labels'].split(delimiter)
    except (csv.Error, UnicodeError, ValueError):
        warn("invalid mapping format. Use a UTF-8 tab-separated file with 'phenotype' "
             "and 'catalog_labels' columns, one row per phenotype, and a nonempty "
             'label delimiter. Step 5 cannot continue.')
        return None

    for phenotype, labels in mapping.items():
        missing = [label for label in labels if label not in catalog_traits]
        if missing:
            warn(f'phenotype={phenotype!r}: labels not found in GWAS Catalog '
                 f'DISEASE/TRAIT: {missing!r}. Check label names and the delimiter '
                 f'(currently {delimiter!r}; configurable via the optional final argument). '
                 'Unmatched labels will not contribute to replication counts.')
    return mapping


def main(args):
    ############################################################
    # 1. Load input
    ############################################################
    all_locus_path = args[1] # ./output/step4.all_locus.csv

    trait_2_catalog_path = args[2] # ./input/trait_2_catalog.csv

    gwas_catalog_path = args[3] # ./external/gwas_catalog_v1.0-associations_e116_r2026-07-19_full.zip

    replication_window = int(args[4]) # 500000

    output_path = args[5] # ./output/

    # Separator within catalog_labels; mapping-file columns always use tabs.
    catalog_label_delimiter = args[6] if len(args) > 6 else '; '

    all_locus_summary = pd.read_csv(all_locus_path, sep='\t') # phenotype chrom start end locusID UID

    gwas_catalog = pd.read_csv(gwas_catalog_path, sep='\t', low_memory=False) # P-VALUE CHR_ID CHR_POS DISEASE/TRAIT

    trait_2_catalog = load_trait_mapping(
        trait_2_catalog_path,
        set(gwas_catalog['DISEASE/TRAIT'].dropna()),
        catalog_label_delimiter,
    )
    if trait_2_catalog is None:
        return 1

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
    all_traits = trait_2_catalog.keys()

    unmapped = ~all_locus_summary['phenotype'].isin(all_traits)
    for trait in all_locus_summary.loc[unmapped, 'phenotype'].unique():
        print(
            f'Warning: phenotype={trait!r} is missing from trait_2_catalog '
            f'{trait_2_catalog_path!r}. Replication counts and candidate_novel '
            'will be NA. Add a mapping and rerun Step 5 to classify these loci.',
            file=sys.stderr,
        )

    # Keep unassessed loci missing, including when no input phenotype is mapped.
    for column in ('rep_ALL', 'rep_EAS', 'rep_EUR'):
        all_locus_summary[column] = np.nan
    all_locus_summary['rep_study'] = ''

    # Evidence from the most significant matching Catalog association.
    all_locus_summary['rep_trait'] = ''
    all_locus_summary['rep_variant'] = ''

    for trait in all_traits:
        # DISEASE/TRAIT label for this phenotype
        this = trait_2_catalog[trait]
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

            if not catalog_in_window.empty:
                # idxmin retains the first Catalog row when P-values are tied.
                top_match = catalog_in_window.loc[catalog_in_window['P-VALUE'].idxmin()]
                all_locus_summary.loc[i, 'rep_trait'] = top_match['DISEASE/TRAIT']
                all_locus_summary.loc[i, 'rep_variant'] = top_match['SNPS']

    ############################################################
    # 3. candidate_novel and save
    ############################################################

    # True: candidate novel; False: known; NA: not assessed.
    all_locus_summary['candidate_novel'] = (
        all_locus_summary['rep_ALL'].eq(0).astype('boolean')
        .mask(all_locus_summary['rep_ALL'].isna())
    )

    # Write an explicit NA for unknown classifications; preserve other output formatting.
    output = all_locus_summary.assign(
        candidate_novel=all_locus_summary['candidate_novel'].astype('string').fillna('NA')
    )
    output.to_csv(os.path.join(output_path, 'step5.all_locus_cata.csv'), sep='\t', index=False)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

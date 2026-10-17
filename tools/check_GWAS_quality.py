#!/usr/bin/env python3
"""Report variant-count and genomic-inflation warnings for one GWAS file.

Usage: python ./tools/check_GWAS_quality.py path/to/summary.tsv.gz
Uses Python only; does not modify the input or filter/output GWAS variants.
"""

import argparse
import csv
import gzip
import sys
from decimal import Decimal, InvalidOperation
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtri_exp


P_ALIASES = ('P', 'PVAL', 'PVALUE', 'P_VALUE', 'P-VALUE')
LOG_ALIASES = ('NEG_LOG10_P', 'LOG10P', 'LOG10_P', 'MLOG10P', 'NEGLOG10P')
NULL_MEDIAN_CHI2 = float(ndtri_exp(np.log(0.25)) ** 2)


def log_pvalues(raw):
    """Convert P strings to ln(P), preserving zero and scientific underflow."""
    numeric = pd.to_numeric(raw, errors='coerce').to_numpy(dtype=float)
    result = np.full(len(raw), np.nan)
    tiny = np.finfo(float).tiny
    normal = np.isfinite(numeric) & (numeric >= tiny) & (numeric <= 1)
    result[normal] = np.log(numeric[normal])
    # Decimal is only needed for zeros/subnormals: e.g. 1e-400 must not
    # silently become P=0, and -1e-400 must not become a valid negative zero.
    small = np.isfinite(numeric) & (np.abs(numeric) < tiny)
    for index in np.flatnonzero(small):
        try:
            value = Decimal(str(raw.iloc[index]).strip())
            if value == 0:
                result[index] = -np.inf
            elif 0 < value <= 1:
                result[index] = float(value.ln())
        except (InvalidOperation, ValueError, OverflowError):
            pass
    return result


def find_column(columns, explicit, aliases):
    if explicit is not None:
        if explicit not in columns:
            raise ValueError(f'Column {explicit!r} is not present in the input.')
        return explicit
    for alias in aliases:
        matches = [column for column in columns if column.upper() == alias]
        if len(matches) > 1:
            raise ValueError(f'Ambiguous columns for {alias}: {matches}. Specify a column explicitly.')
        if matches:
            return matches[0]
    return None


def check_quality(path, p_column=None, log10p_column=None, test_column=None,
                  keep_test='ADD', all_tests=False, chunk_size=200_000):
    """Read only needed columns; use exact in-place selection on numeric ln(P).

    lambda = median(chi2.isf(P, 1)) / chi2.ppf(0.5, 1). For even n, transform
    both central ln(P) values to chi-square BEFORE taking their average.
    """
    path = Path(path)
    compressed = path.suffix.lower() in ('.gz', '.gzip')
    opener = gzip.open if compressed else open
    with opener(path, 'rt', encoding='utf-8-sig') as handle:
        header = next((line.strip() for line in handle if line.strip()), '')
    if not header:
        raise ValueError('The summary-statistics file is empty.')
    separator = ',' if ',' in header else r'\s+'
    columns = next(csv.reader([header])) if separator == ',' else header.split()
    if len(columns) != len(set(columns)):
        raise ValueError('Duplicate column names in the input header.')
    p_column = find_column(columns, p_column, P_ALIASES)
    log10p_column = find_column(columns, log10p_column, LOG_ALIASES)
    test_column = find_column(columns, test_column, ('TEST',))
    if p_column is None and log10p_column is None:
        raise ValueError('No P or -log10(P) column found. Use --p-column or --log10p-column.')
    usecols = list(dict.fromkeys(column for column in
                                (p_column, log10p_column, test_column) if column))
    parts = []
    n_rows = n_variants = n_valid = 0
    with pd.read_csv(path, sep=separator, usecols=usecols, dtype=str,
                     keep_default_na=False, encoding='utf-8-sig',
                     compression='gzip' if compressed else None,
                     chunksize=chunk_size) as chunks:
        for chunk in chunks:
            n_rows += len(chunk)
            if test_column is not None and not all_tests:
                chunk = chunk.loc[chunk[test_column].str.strip().eq(keep_test)]
            n_variants += len(chunk)
            ln_p = (log_pvalues(chunk[p_column]) if p_column is not None
                    else np.full(len(chunk), np.nan))
            if log10p_column is not None:
                log10p = pd.to_numeric(chunk[log10p_column], errors='coerce').to_numpy(dtype=float)
                valid_log = np.isfinite(log10p) & (log10p >= 0)
                with np.errstate(over='ignore'):
                    ln_p[valid_log] = -log10p[valid_log] * np.log(10)
            valid = np.isfinite(ln_p) | np.isneginf(ln_p)
            n_valid += int(valid.sum())
            if valid.any():
                parts.append(ln_p[valid])
    inflation = None
    median_chi2 = None
    if parts:
        values = np.concatenate(parts)
        parts.clear()
        lo, hi = (n_valid - 1) // 2, n_valid // 2
        # Partition is exact selection, without sorting every value.
        values.partition((lo, hi))
        central_chi2 = ndtri_exp(values[[lo, hi]] - np.log(2)) ** 2
        median_chi2 = float(np.mean(central_chi2))
        inflation = median_chi2 / NULL_MEDIAN_CHI2
    return dict(path=str(path), n_rows=n_rows, n_variants=n_variants,
                n_valid=n_valid, n_invalid=n_variants - n_valid,
                lambda_gc=inflation, median_chi2=median_chi2,
                p_column=p_column, log10p_column=log10p_column,
                test_column=test_column,
                test_selection=keep_test if test_column and not all_tests else 'all rows')


def quality_warnings(n_variants, n_valid, inflation):
    messages = []
    if n_variants < 1_000_000 or n_variants > 10_000_000:
        messages.append('N variants is below 1 million or above 10 million; '
                        'check summary-statistics completeness and variant filtering.')
    if inflation is None:
        messages.append('Cannot calculate lambda: no valid P values.')
    elif inflation > 1.2:
        messages.append('Genomic inflation factor lambda is greater than 1.2; '
                        'check for potential statistical inflation.')
    if n_valid < n_variants:
        messages.append(f'Excluded {n_variants - n_valid:,} rows with missing or invalid '
                        'P values from lambda calculation.')
    return messages


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('summary_file', type=Path, help='Whitespace/TSV or CSV input; .gz/.gzip supported.')
    parser.add_argument('--p-column', help='P-value column; otherwise detect common names such as P/PVAL/PVALUE.')
    parser.add_argument('--log10p-column', help='-log10(P) column; otherwise detect NEG_LOG10_P/LOG10P and aliases.')
    parser.add_argument('--test-column', help='Test-type column; otherwise use TEST if present.')
    parser.add_argument('--keep-test', default='ADD', help='Test type to keep (default: ADD).')
    parser.add_argument('--all-tests', action='store_true', help='Use all rows even when TEST is present.')
    args = parser.parse_args(argv)
    try:
        result = check_quality(args.summary_file, args.p_column, args.log10p_column,
                               args.test_column, args.keep_test, args.all_tests)
    except (OSError, EOFError, ValueError, pd.errors.ParserError) as error:
        parser.exit(2, f'Error: {error}\n')
    print(f"File: {result['path']}")
    print(f"P column: {result['p_column'] or 'none'}; -log10(P) column: {result['log10p_column'] or 'none'}")
    selection = (f"{result['test_column']}={result['test_selection']}"
                 if result['test_selection'] != 'all rows' else 'all rows')
    print(f'Test selection: {selection}')
    print(f"Input rows: {result['n_rows']:,}")
    print(f"N variants: {result['n_variants']:,}")
    print(f"Valid P values: {result['n_valid']:,}")
    print(f"Missing/invalid P values excluded from lambda: {result['n_invalid']:,}")
    inflation = result['lambda_gc']
    print('lambda_GC: ' + ('unavailable' if inflation is None else f'{inflation:.16g}'))
    print('Rule: median(chi2.isf(P, df=1)) / 0.4549364231195727')
    for message in quality_warnings(result['n_variants'], result['n_valid'], inflation):
        print(f'[Warning] {message}')
    return 0


if __name__ == '__main__':
    sys.exit(main())

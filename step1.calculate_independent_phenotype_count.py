#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os
import sys
from scipy.stats import rankdata


def warn_phenotype_matrix(data, missing_threshold=0.20, input_path=None):
    """Report input concerns before calculation without changing data or QC.

    Rows should be unique participants and columns should be phenotypes.
    Complete, unique columns with no float-like values are possible IDs.
    This heuristic only warns; the original numeric selection and QC decide
    what enters the calculation. The input DataFrame is never modified.
    """
    def warn(message):
        print(f"[Warning] {message}", file=sys.stderr)

    # Read the original header because pandas renames duplicate/empty names.
    header = list(data.columns)
    if input_path is not None:
        header = pd.read_csv(
            input_path, sep="\t", header=None, nrows=1,
            dtype=str, keep_default_na=False,
        ).iloc[0].tolist()
    names = pd.Index([str(name).strip() for name in header])
    duplicates = names[names.duplicated()].unique().tolist()
    if duplicates:
        warn("Duplicate column names: " + ", ".join(repr(c) for c in duplicates)
             + ". Check the input header.")
    if any(not name for name in names):
        warn("Empty column names detected. Check the input header.")
    if data.shape[1] < 2:
        warn("Expected a tab-delimited matrix with at least two phenotype columns.")
    if not isinstance(data.index, pd.RangeIndex):
        warn("Unexpected row index; check tab delimiters and field counts.")
    if len(data) < 3:
        warn("Fewer than 3 participant rows; check matrix orientation and sample size.")
    if not np.isfinite(missing_threshold) or not 0 <= missing_threshold <= 1:
        warn("Missing-rate threshold is outside [0, 1]. Check the argument.")

    suspected_ids, nonnumeric, missing_text = [], [], []
    infinite, excluded, imputed, rare_binary = [], [], [], []
    numeric_count = 0
    retained_count = 0
    for column in data.columns:
        series = data[column]
        text = series.astype("string").str.strip()
        numeric = pd.to_numeric(text, errors="coerce")
        float_like = numeric.notna() & ~text.str.fullmatch(r"[+-]?\d+", na=False)
        if len(series) and series.notna().all() and series.is_unique and not float_like.any():
            suspected_ids.append(str(column))
        if not pd.api.types.is_numeric_dtype(series):
            nonnumeric.append(str(column))
            tokens = {"", "NA", "N/A", "NAN", "NULL", "NONE", "<NA>"}
            if (series.notna() & text.str.upper().isin(tokens)).any():
                missing_text.append(str(column))
            continue

        numeric_count += 1
        values = series.to_numpy(dtype=float)
        observed = values[~np.isnan(values)]
        n_observed = len(observed)
        missing_rate = float(np.isnan(values).mean()) if len(values) else np.nan
        with np.errstate(over="ignore", invalid="ignore", under="ignore"):
            sd = float(np.std(observed, ddof=1)) if n_observed >= 2 else np.nan
        reasons = []
        if np.isinf(values).any():
            infinite.append(str(column))
        if missing_rate > missing_threshold:
            reasons.append("high missing rate")
        if n_observed < 3:
            reasons.append("<3 observations")
        if n_observed and np.unique(observed).size == 1:
            reasons.append("constant")
        if n_observed >= 3 and not np.isfinite(sd):
            reasons.append("nonfinite SD")
        elif n_observed >= 3 and sd == 0 and "constant" not in reasons:
            reasons.append("zero SD")
        # Mirror the existing QC criteria for reporting only.
        retained = (
            missing_rate <= missing_threshold and n_observed >= 3
            and np.isfinite(sd) and sd > 0
        )
        if retained:
            retained_count += 1
            if np.isnan(values).any():
                imputed.append(str(column))
        elif reasons:
            excluded.append(f"{column} ({'; '.join(reasons)})")
        if np.array_equal(np.unique(observed), [0, 1]):
            cases = int(np.count_nonzero(observed == 1))
            if cases < 10:
                rare_binary.append(f"{column} ({cases})")

    if suspected_ids:
        warn("Possible ID columns (unique, non-float values): "
             + ", ".join(suspected_ids) + ". Remove IDs before analysis.")
    if nonnumeric:
        warn("Nonnumeric columns will be ignored: " + ", ".join(nonnumeric) + ".")
    if missing_text:
        warn("Unrecognized missing-value text: " + ", ".join(missing_text)
             + ". Normalize to NA or empty cells.")
    if infinite:
        warn("Infinite values detected: " + ", ".join(infinite) + ". Check input values.")
    if excluded:
        warn("QC will exclude: " + ", ".join(excluded) + ".")
    if imputed:
        warn("Missing values will be median-imputed: " + ", ".join(imputed) + ".")
    if rare_binary:
        warn("0/1 columns with <10 cases (1=case): "
             + ", ".join(rare_binary) + ". Review before analysis.")
    duplicate_rows = int(data.duplicated().sum())
    if duplicate_rows:
        warn(f"Repeated rows: {duplicate_rows}; confirm participant uniqueness.")
    if numeric_count < 2:
        warn("Fewer than 2 numeric phenotype columns are available.")
    elif retained_count < 2:
        warn("Fewer than 2 phenotypes will remain after QC.")


def calculate_effective_markers(
    data,
    marker_cols=None,
    id_cols=None,
    missing_threshold=0.20,
    cor_method="spearman",
    alpha=5e-8,
    pc_cutoffs=(0.95, 0.99, 0.995),
):
    """
    Calculate effective number of markers using Li-Ji and Galwey methods.

    Parameters
    ----------
    - data : pandas.DataFrame
    - marker_cols : list or None
    - id_cols : list or None
    - missing_threshold : float
    - cor_method : "spearman" or "pearson"
    - alpha : float, default 5e-8
        Genome-wide significance level before correction across phenotypes.
    - pc_cutoffs : iterable
    """
    cor_method = cor_method.lower()
    if cor_method not in ("pearson", "spearman"):
        raise ValueError("cor_method must be 'pearson' or 'spearman'.")
    df = data.copy()

    ############################################################
    # 1. Marker selection
    ############################################################

    if marker_cols is None:
        if id_cols is None:
            id_cols = []
        candidate_cols = [
            c for c in df.columns
            if c not in id_cols
        ]
        marker_cols = [
            c for c in candidate_cols
            if pd.api.types.is_numeric_dtype(df[c])
        ]
    if len(marker_cols) < 2:
        raise ValueError("At least two markers are required.")
    missing_cols = list(set(marker_cols) - set(df.columns))
    if len(missing_cols) > 0:
        raise ValueError(
            f"Marker columns not found: {missing_cols}"
        )
    X = df[marker_cols].to_numpy(dtype=float)
    n_subjects = X.shape[0]
    n_markers_raw = X.shape[1]
    if n_subjects < 3:
        raise ValueError("Sample size must be >=3.")
    
    ############################################################
    # 2. QC
    ############################################################

    missing_rate = np.mean(np.isnan(X), axis=0)
    non_missing_n = np.sum(~np.isnan(X), axis=0)
    marker_sd = np.nanstd(
        X,
        axis=0,
        ddof=1
    )
    keep = (
        (missing_rate <= missing_threshold)
        &
        (non_missing_n >= 3)
        &
        np.isfinite(marker_sd)
        &
        (marker_sd > 0)
    )
    marker_qc = pd.DataFrame({
        "marker": marker_cols,
        "missing_rate": missing_rate,
        "non_missing_n": non_missing_n,
        "standard_deviation": marker_sd,
        "retained": keep
    })
    X = X[:, keep]
    retained_markers = np.array(marker_cols)[keep].tolist()
    n_markers_qc = X.shape[1]
    if n_markers_qc < 2:
        raise ValueError(
            "Less than two markers remain after QC."
        )
    print(f"Number of subjects: {n_subjects}")
    print(f"Raw marker number: {n_markers_raw}")
    print(f"Markers retained after QC: {n_markers_qc}")
    print(f"Correlation method: {cor_method}")

    ############################################################
    # 3. Median imputation
    ############################################################

    for j in range(X.shape[1]):
        mask = np.isnan(X[:, j])
        if np.any(mask):
            median = np.nanmedian(X[:, j])
            X[mask, j] = median
    
    ############################################################
    # 4. Spearman ranking
    ############################################################

    if cor_method == "spearman":
        X = np.apply_along_axis(
            lambda x: rankdata(
                x,
                method="average"
            ),
            axis=0,
            arr=X
        )
    
    ############################################################
    # 5. Standardization
    ############################################################

    mean = X.mean(axis=0)
    sd = X.std(
        axis=0,
        ddof=1
    )
    X = (X - mean) / sd

    ############################################################
    # 6. Eigenvalues
    ############################################################

    n, p = X.shape
    if p <= n:
        correlation_matrix = X.T @ X / (n - 1)
        eigenvalues = np.linalg.eigvalsh(
            correlation_matrix
        )
    else:
        subject_matrix = X @ X.T / (n - 1)
        nonzero = np.linalg.eigvalsh(
            subject_matrix
        )
        eigenvalues = np.concatenate([
            nonzero,
            np.zeros(p - len(nonzero))
        ])

    ############################################################
    # Numerical correction
    ############################################################

    eigenvalues[
        (eigenvalues < 0)
        &
        (eigenvalues > -1e-8)
    ] = 0
    if np.any(eigenvalues < -1e-8):
        print(
            "Warning: negative eigenvalues detected."
        )
    eigenvalues[eigenvalues < 0] = 0
    eigenvalues = eigenvalues[::-1]
    eigenvalue_sum = eigenvalues.sum()

    ############################################################
    # 7. Li-Ji
    ############################################################

    me_li_ji = np.minimum(
        eigenvalues,
        1
    ).sum()

    ############################################################
    # 8. Galwey
    ############################################################

    me_galwey = (
        np.sum(np.sqrt(eigenvalues)) ** 2
    ) / eigenvalue_sum

    ############################################################
    # 9. PCA
    ############################################################

    variance_proportion = (
        eigenvalues
        /
        eigenvalue_sum
    )
    cumulative_variance = np.cumsum(
        variance_proportion
    )
    pc_numbers = []
    for cutoff in pc_cutoffs:
        idx = np.argmax(
            cumulative_variance >= cutoff
        )
        if cumulative_variance[idx] >= cutoff:
            pc_numbers.append(idx + 1)
        else:
            pc_numbers.append(len(eigenvalues))
    
    ############################################################
    # 10. Multiple testing
    ############################################################

    bonferroni_raw = alpha / n_markers_qc
    bonferroni_li_ji = alpha / me_li_ji
    bonferroni_galwey = alpha / me_galwey
    sidak_raw = 1 - (1 - alpha) ** (
        1 / n_markers_qc
    )
    sidak_li_ji = 1 - (1 - alpha) ** (
        1 / me_li_ji
    )
    sidak_galwey = 1 - (1 - alpha) ** (
        1 / me_galwey
    )

    ############################################################
    # 11. Summary
    ############################################################

    summary = pd.DataFrame({
        "method": [
            "Raw markers after QC",
            "Li-Ji effective markers",
            "Galwey effective markers",
        ],
        "effective_marker_number": [
            n_markers_qc,
            me_li_ji,
            me_galwey,
        ],
        "bonferroni_threshold": [
            bonferroni_raw,
            bonferroni_li_ji,
            bonferroni_galwey,
        ],
        "sidak_threshold": [
            sidak_raw,
            sidak_li_ji,
            sidak_galwey,
        ],
    })
    principal_components = pd.DataFrame({
        "cumulative_variance_cutoff": pc_cutoffs,
        "number_of_PCs": pc_numbers,
        "pc_based_threshold": [
            alpha / x
            for x in pc_numbers
        ],
    })
    eigenvalue_table = pd.DataFrame({
        "component": np.arange(
            1,
            len(eigenvalues) + 1
        ),
        "eigenvalue": eigenvalues,
        "variance_proportion": variance_proportion,
        "cumulative_variance": cumulative_variance,
    })
    return {
        "summary": summary,
        "principal_components": principal_components,
        "marker_qc": marker_qc,
        "eigenvalues": eigenvalue_table,
        "retained_markers": retained_markers,
        "n_subjects": n_subjects,
        "n_markers_raw": n_markers_raw,
        "n_markers_qc": n_markers_qc,
    }

def main(args):

    ############################################################
    # 1. Load input
    ############################################################

    phenotype_matrix_path = args[1] # './input/phenotype_matrix.csv'

    output_path = args[2] # './output/'

    missing_rate_threshold = float(args[3]) # 0.2

    correlation_method = args[4] # spearman

    phenotype_matrix = pd.read_csv(phenotype_matrix_path, sep='\t')
    warn_phenotype_matrix(
        phenotype_matrix, missing_threshold=missing_rate_threshold,
        input_path=phenotype_matrix_path,
    )

    ############################################################
    # 2. Calculation
    ############################################################

    res = calculate_effective_markers(
        data=phenotype_matrix,
        cor_method=correlation_method,
        missing_threshold=missing_rate_threshold,
        alpha=0.05,
    )
    
    print("Results:")
    print(res['summary'])

    M_eff = res['summary'].loc[1, 'effective_marker_number']
    adjusted_threshold = res['summary'].loc[1, 'bonferroni_threshold']
    print('M_eff: %.2f'%M_eff)
    print('Adjusted study-wide significant p-value: %.2e'%adjusted_threshold)

    ############################################################
    # 3. Visualization and output
    ############################################################

    res['summary'].to_csv(os.path.join(output_path, 'step1.summary.csv'), sep='\t', index=False)
    res['principal_components'].to_csv(os.path.join(output_path, 'step1.PC.csv'), sep='\t', index=False)
    res['marker_qc'].to_csv(os.path.join(output_path, 'step1.markers.csv'), sep='\t', index=False)
    res['eigenvalues'].to_csv(os.path.join(output_path, 'step1.eigenvalues.csv'), sep='\t', index=False)

if __name__ == "__main__":
    main(sys.argv)

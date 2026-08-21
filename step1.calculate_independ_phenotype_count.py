#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os
import sys
from scipy.stats import rankdata

def calculate_effective_markers(
    data,
    marker_cols=None,
    id_cols=None,
    missing_threshold=0.20,
    cor_method="spearman",
    alpha=0.05,
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
    - alpha : float
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
    print('M_eff: %.2f'%M_eff)
    print('Adjusted study-wide significant p-value: %.2e'%(5e-8/M_eff))

    ############################################################
    # 3. Visualization and output
    ############################################################

    res['summary'].to_csv(os.path.join(output_path, 'step1.summary.csv'), sep='\t', index=False)
    res['principal_components'].to_csv(os.path.join(output_path, 'step1.PC.csv'), sep='\t', index=False)
    res['marker_qc'].to_csv(os.path.join(output_path, 'step1.markers.csv'), sep='\t', index=False)
    res['eigenvalues'].to_csv(os.path.join(output_path, 'step1.eigenvalues.csv'), sep='\t', index=False)

if __name__ == "__main__":
    main(sys.argv)

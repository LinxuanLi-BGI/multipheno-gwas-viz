#!/usr/bin/env bash
# Linux/WSL: bash example/run_example_e2e.sh
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."

PYTHON=${PYTHON:-python3}
OUT=./tests/
GENES=external/KnownCanonicalGene.hg38.txt
CATALOG=external/gwas_catalog_v1.0-associations_e116_r2026-07-19_full.zip
mkdir -p ./tests/
export MPLBACKEND=Agg MPLCONFIGDIR="$OUT/.matplotlib" PYTHONDONTWRITEBYTECODE=1

"$PYTHON" step1.calculate_independent_phenotype_count.py example/example_for_Phenotype_matrix.csv "$OUT" 0.2 spearman
# Fixed threshold for the five-summary GWAS example.
# The synthetic Step1 phenotype matrix is a separate demonstration.
P_THRESHOLD=1e-8
"$PYTHON" step2.extract_significant_snp.py example/example_for_GWAS_manifest.csv "$OUT" example/summary_config.json "$P_THRESHOLD" "${WORKERS:-5}"
"$PYTHON" step3.mark_gene.py "$OUT/step2.all_sig.csv" "$GENES" "$OUT"
"$PYTHON" step4.snp2locus.py "$OUT/step3.all_sig_gene.csv" "$GENES" 500000 "$OUT"
"$PYTHON" step5.mark_known_GWAS_catalog.py "$OUT/step4.all_locus.csv" example/example_for_trait_2_catalog.csv "$CATALOG" 500000 "$OUT"
"$PYTHON" step6.draw_circos_plot.py "$OUT/step5.all_locus_cata.csv" example/color_config.json "$OUT" png candidate_novel
"$PYTHON" step7.mark_GWAS_catalog_background.py "$OUT/step5.all_locus_cata.csv" "$CATALOG" 500000 "$OUT"

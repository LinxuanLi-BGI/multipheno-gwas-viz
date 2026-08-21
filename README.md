# Multipheno GWAS Visualization

A Python pipeline for identifying significant GWAS variants, annotating genes and loci, comparing results with the GWAS Catalog, and generating a multi-phenotype circos plot.

## Requirements

- Python 3
- `pandas`
- `numpy`
- `scipy`
- `matplotlib`
- `awk` (used by step 2)

Install Python dependencies with:

```bash
pip install -r requirements.txt
```

## Pipeline

Run the steps in order:

1. Calculate effective phenotype marker counts.
2. Extract significant SNPs from GWAS summary statistics.
3. Annotate SNPs with known genes.
4. Group SNPs into loci and summarize cross-phenotype recurrence.
5. Mark loci reported in the GWAS Catalog.
6. Generate the circos plot.
7. Add GWAS Catalog background annotations.

## Project Structure

```text
config/     Configuration files
external/   Reference data
input/      Input data
output/     Generated results
step*.py    Pipeline scripts
```

## Usage

```bash
python step1.calculate_independ_phenotype_count.py <phenotype_matrix> <output_dir> <missing_rate> <correlation_method>
python step2.extract_significanct_snp.py <manifest> <output_dir> <summary_config> <pvalue_threshold> <workers>
python step3.mark_gene.py <input_file> <known_gene_file> <output_dir>
python step4.snp2locus.py <input_file> <known_gene_file> <locus_window> <output_dir>
python step5.mark_known_GWAS_catalog.py <locus_file> <trait_mapping> <gwas_catalog> <replication_window> <output_dir>
python step6.draw_circos_plot.py <locus_file> <color_config> <output_dir> <format> <candidate_novel>
python step7.mark_GWAS_catalog_background.py <locus_file> <gwas_catalog> <replication_window> <output_dir>
```

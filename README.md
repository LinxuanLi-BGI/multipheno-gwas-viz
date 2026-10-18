# Multipheno GWAS Visualization

A Python pipeline for identifying significant GWAS variants, annotating genes and loci, comparing results with the GWAS Catalog, and generating a multi-phenotype circos plot.

## Requirements

| Dependency | Version |
| --- | --- |
| Python | 3.15.0 |
| pandas | 3.0.6 |
| NumPy | 2.5.4 |
| SciPy | 1.18.1 |
| Matplotlib | 3.11.2 |

Create and activate a dedicated Conda environment from the repository root:

```bash
conda env create -f environment.yml
conda activate multipheno-gwas-viz
```

Alternatively, with Python 3.15.0 already installed, install the Python dependencies using pip:

```bash
python -m pip install -r requirements.txt
```

Step 2 also requires `awk`, `cat`, and `zcat` on `PATH` in a Unix-like shell environment (for example, linux or WSL on Windows). These system tools are separate from the Python dependencies in `environment.yml` and `requirements.txt`.

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
example/    Example inputs and JSON configuration files
external/   Reference data
input/      User-prepared input data (create as needed)
output/     Generated results (create before running)
step*.py    Pipeline scripts
tools/      Standalone utility scripts
```

## Usage

```bash
python step1.calculate_independent_phenotype_count.py <phenotype_matrix> <output_dir> <missing_rate> <correlation_method>
python step2.extract_significant_snp.py <manifest> <output_dir> <summary_config> <pvalue_threshold> <workers>
python step3.mark_gene.py <input_file> <known_gene_file> <output_dir>
python step4.snp2locus.py <input_file> <known_gene_file> <locus_window> <output_dir>
python step5.mark_known_GWAS_catalog.py <locus_file> <trait_mapping> <gwas_catalog> <replication_window> <output_dir> [catalog_label_delimiter]
python step6.draw_circos_plot.py <locus_file> <color_config> <output_dir> <format> [highlight_column]
python step7.mark_GWAS_catalog_background.py <locus_file> <gwas_catalog> <replication_window> <output_dir>
```

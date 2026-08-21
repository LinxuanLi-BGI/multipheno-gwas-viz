#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json
import pandas as pd
import numpy as np
import subprocess
import shlex
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

builtin_summary_config = {
    'PLINK2':{
        "chrom" : "#CHROM",
        "pos" : "POS",
        "A1" : "A1",
        "A2" : "OMITTED",
        "freq" : "A1_FREQ",
        "size" : "OBS_CT",
        "beta" : "BETA",
        "or" : "OR",
        "se" : 'SE',
        "pvalue" : "P",
        "log10p" : "NEG_LOG10_P",
        "test" : "TEST",
        "keep" : "ADD"
    },
    'PLINK2_logistics':{
        "chrom" : "#CHROM",
        "pos" : "POS",
        "A1" : "A1",
        "A2" : "OMITTED",
        "freq" : "A1_FREQ",
        "size" : "OBS_CT",
        "beta" : "",
        "or" : "OR",
        "se" : 'LOG(OR)_SE',
        "pvalue" : "P",
        "log10p" : "NEG_LOG10_P",
        "test" : "TEST",
        "keep" : "ADD"
    },
    'Regenie':{
        "chrom" : "CHROM",
        "pos" : "GENPOS",
        "A1" : "ALLELE1",
        "A2" : "ALLELE0",
        "freq" : "A1FREQ",
        "size" : "N",
        "beta" : "BETA",
        "or" : "",
        "se" : "SE",
        "pvalue" : "",
        "log10p" : "LOG10P",
        "test" : "TEST",
        "keep" : "ADD"
    }
}

def execute_command(command):
    """
    Start a parallel worker and execute the command

    Parameters
    ----------
    - command: string, command to execute
    """
    try:
        # execute the command
        result = subprocess.run(
            command,
            shell=True,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding='utf-8',
            timeout=None
        )
        # return the results
        return {
            'command': command,
            'pid': os.getpid(),
            'returncode': result.returncode,
            'stdout': result.stdout,
            'stderr': result.stderr,
            'status': 'SUCCESS' if result.returncode == 0 else 'FAILED'
        }
    except Exception as e:
        # catch the exception
        return {
            'command': command,
            'pid': os.getpid(),
            'returncode': -1,
            'stdout': '',
            'stderr': str(e),
            'status': 'EXCEPTION'
        }

def load_config(input_str: str) -> dict:
    # check if input matches builtin config keys
    if input_str in builtin_summary_config:
        return builtin_summary_config[input_str]

    path_obj = Path(input_str)
    suffix = path_obj.suffix.lower()

    # parse json file
    if suffix == ".json":
        with open(path_obj, "r", encoding="utf-8") as f:
            return json.load(f)

    # parse csv file, use column 0 as index, column 1 as value to build dict
    if suffix == ".csv":
        df = pd.read_csv(path_obj, index_col=0)
        # take the first data column as dict value
        first_col_name = df.columns[0]
        df[first_col_name] = df[first_col_name].fillna("").astype(str)
        return df[first_col_name].to_dict()

    raise ValueError(
        f"Unsupported input: {input_str}. "
        "Must be PLINK2/Regenie builtin key, .json path or .csv path."
    )

# format chrom
def reformat_chrom(x:str)->int:
    try:
        # int type
        if isinstance(x, (int, np.integer, float, np.floating)):
            return int(x)
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

def main(args):

    ############################################################
    # 1. Load input
    ############################################################

    GWAS_manifest_file = args[1] # ./input/GWAS_manifest.csv

    output_path = args[2] # ./output/

    summary_config_key = args[3] # PLINK

    significance_threshold = float(args[4]) # 1e-10

    PROCESS_POOL_SIZE = int(args[5]) # 24

    # GWAS_manifest col: phenotype_name, phenotype_category, summary_statistics_path
    GWAS_manifest = pd.read_csv(GWAS_manifest_file, sep='\t')

    # Path to store significant SNPs extracted from each phenotype
    try:
        os.makedirs(os.path.join(output_path, 'sig'), exist_ok=True)
    except Exception as E:
        print(E)

    # summary_config
    summary_config = load_config(summary_config_key)

    # Significance threshold
    log_t = -np.log10(significance_threshold)

    ############################################################
    # 2. Extract significant SNP in batch
    ############################################################

    # Use the header to find columns, allowing input files to have different
    # column orders while retaining their original header in the output.
    awk_script = (
        'BEGIN { FS="[ \\t]+"; OFS="\\t" } '
        'NR==1 { for (i=1; i<=NF; i++) col[$i]=i; '
        'p=col[pname]; l=col[lname]; t=col[tname]; print; next } '
        '((pname == "" && lname == "") || '
        '(pname != "" && p && $(p)+0 < threshold) || '
        '(lname != "" && l && $(l)+0 > log_threshold)) '
        '&& (tname == "" || keep_value == "" || (t && $(t) == keep_value)) '
        '{ print }'
    )

    def make_task(row):
        input_path = str(row['summary_statistics_path'])
        input_file = shlex.quote(input_path)
        output_file = shlex.quote(os.path.join(
            output_path, 'sig', f"{row['phenotype_name']}.sig"
        ))
        options = [
            f"-v threshold={shlex.quote(str(significance_threshold))}",
            f"-v log_threshold={shlex.quote(str(log_t))}",
            f"-v pname={shlex.quote(str(summary_config.get('pvalue') or ''))}",
            f"-v lname={shlex.quote(str(summary_config.get('log10p') or ''))}",
            f"-v tname={shlex.quote(str(summary_config['test']))}",
            f"-v keep_value={shlex.quote(str(summary_config['keep']))}",
        ]
        # Stream gzip-compressed inputs through zcat; awk reads regular files
        # through cat so both input formats use the same filtering logic.
        reader = 'zcat' if input_path.lower().endswith(('.gz', '.gzip')) else 'cat'
        return f"{reader} -- {input_file} | awk {' '.join(options)} {shlex.quote(awk_script)} > {output_file}"

    GWAS_manifest['task'] = GWAS_manifest.apply(make_task, axis=1)

    task_list = GWAS_manifest['task'].values

    print(f"Starting tasks; process pool size: {PROCESS_POOL_SIZE}")
    print(f"Total tasks: {len(task_list)}\n")

    # Create a process pool with the specified maximum number of workers
    with ProcessPoolExecutor(max_workers=PROCESS_POOL_SIZE) as executor:
        # Submit all tasks to the process pool and map futures to commands
        future_to_command = {executor.submit(execute_command, cmd): cmd for cmd in task_list}

        # Iterate over completed tasks in completion order
        for future in as_completed(future_to_command):
            command = future_to_command[future]
            try:
                # Get the task result
                result = future.result()
                # Print the execution result
                print(f"[{result['status']}] Process ID: {result['pid']} | Command: {command}")
                if result['stderr']:
                    print(f"  Error output: {result['stderr']}")
            except Exception as e:
                # Catch exceptions raised while retrieving the result
                print(f"[RESULT EXCEPTION] Command: {command} | Exception: {str(e)}")

    print("\nAll tasks completed!")

    ############################################################
    # 3. Load all sig, format and concat
    ############################################################

    sig_frames = []
    sig_dir = Path(output_path) / 'sig'

    for _, manifest_row in GWAS_manifest.iterrows():
        sig_file = sig_dir / f"{manifest_row['phenotype_name']}.sig"
        if not sig_file.exists():
            continue

        sig_df = pd.read_csv(sig_file, sep='\t')
        formatted = pd.DataFrame(index=sig_df.index)

        # Map source summary-statistics columns to the standard output names.
        for output_name, source_name in summary_config.items():
            if output_name in ('pvalue', 'log10p', 'beta', 'or', 'test', 'keep') or not source_name:
                continue
            if source_name in sig_df.columns:
                formatted[output_name] = sig_df[source_name]

        # Keep only beta; fill missing beta values from OR using exp(OR).
        beta_name = summary_config.get('beta') or ''
        or_name = summary_config.get('or') or ''
        if beta_name in sig_df.columns:
            beta = pd.to_numeric(sig_df[beta_name], errors='coerce')
        elif or_name in sig_df.columns:
            odds_ratio = pd.to_numeric(sig_df[or_name], errors='coerce')
            beta = np.log(odds_ratio)
        else:
            print('[Warning]: No beta in summary!')
            beta = pd.Series(np.nan, index=sig_df.index)
        formatted['beta'] = beta

        # Handle the pvalue or log10(pvalue) col
        p_name = summary_config.get('pvalue') or ''
        log_name = summary_config.get('log10p') or ''

        # Always output both forms. A zero p-value is represented as 308.
        if log_name in sig_df.columns:
            log10p = pd.to_numeric(sig_df[log_name], errors='coerce')
            pvalue:pd.Series = np.power(10,-log10p)
            pvalue = pvalue.mask(log10p > 308, 0)
        elif p_name in sig_df.columns:
            pvalue = pd.to_numeric(sig_df[p_name], errors='coerce')
            log10p:pd.Series = -np.log10(pvalue.where(pvalue > 0))
            log10p = log10p.mask(pvalue.eq(0), 308)
        else:
            print('[Warning]: No pvalue in summary!')
            pvalue = pd.Series(np.nan, index=sig_df.index)
            log10p = pd.Series(np.nan, index=sig_df.index)
        formatted['pvalue'] = pvalue
        formatted['log10p'] = log10p


        # format chrom code to int type
        if 'chrom' in formatted.columns:
            formatted['chrom'] = formatted['chrom'].map(reformat_chrom)
        else:
            pass

        formatted['phenotype'] = manifest_row['phenotype_name']
        formatted['category'] = manifest_row['phenotype_category']
        sig_frames.append(formatted)

    if sig_frames:
        all_sig = pd.concat(sig_frames, ignore_index=True)
        all_sig.to_csv(Path(output_path) / 'step2.all_sig.csv', sep='\t', index=False)
    else:
        print('No significant SNP extracted.')

if __name__ == "__main__":
    main(sys.argv)

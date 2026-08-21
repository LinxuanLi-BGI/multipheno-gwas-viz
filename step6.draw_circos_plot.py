#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib.colors import to_hex
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

bulitin_palettes = {
    'Accent': 8 ,
    'Dark2' : 8 , 
    'Set1'  : 9 , 
    'Set2'  : 8 , 
    'Set3'  : 12, 
    'tab10' : 10,
}

# Distance between chromosomes
Chr_sep = 10000000
# Arc range drawn around the circle
Dim = 1.5 * np.pi
# max gene mark on plot
N_gene_mark = 30

Chr_len = {
    1: 248956422,
    2: 242193529,
    3: 198295559,
    4: 190214555,
    5: 181538259,
    6: 170805979,
    7: 159345973,
    8: 145138636,
    9: 138394717,
    10: 133797422,
    11: 135086622,
    12: 133275309,
    13: 114364328,
    14: 107043718,
    15: 101991189,
    16: 90338345,
    17: 83257441,
    18: 80373285,
    19: 58617616,
    20: 64444167,
    21: 46709983,
    22: 50818468,
    23: 156040895,
    24: 57227415,
    25: 16569,
}

def load_config(input_str: str) -> dict:
    # check if input matches builtin config keys
    if input_str in bulitin_palettes:
        return {}

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
        "Must be Palettes builtin key, .json path or .csv path."
    )

def rotate_text(x_loc):
    if x_loc <= np.pi/Dim :
        return x_loc
    else:
        return x_loc - np.pi/Dim

def strip_arrow(s: str) -> str:
    """
    Remove leading '<-' and trailing '->' from string.
    """
    if type(s) is not str:
        return s
    if s.startswith("<-"):
        s = s[2:]
    if s.endswith("->"):
        s = s[:-2]
    return s

def main(args):
    ############################################################
    # 1. Load input
    ############################################################

    all_locus_path = args[1] # ./output/step5.all_locus_cata.csv
    palette_path = args[2] # Dark2 
    output_path = args[3] # ./output/
    output_format = args[4] # png
    try:
        highlight_marker = args[5] # candidate_novel
    except:
        highlight_marker = 'NA'

    all_locus_summary:pd.DataFrame = pd.read_csv(all_locus_path, sep='\t') # phenotype category chrom top UID

    if output_format not in ['png','svg','pdf']:
        print(f"Unsupported output format: {output_format}; using png instead.")
        output_format = 'png'
    # highlight_marker
    if highlight_marker == 'NA':
        pass
    elif highlight_marker not in all_locus_summary.columns:
        print(f"Highlight marker not found: {highlight_marker}")
        highlight_marker = 'NA'
    elif all_locus_summary[highlight_marker].nunique() != 2:
        print(f"Highlight marker N: {all_locus_summary[highlight_marker].nunique()}")
    else:
        # OK
        pass

    ############################################################
    # 2. phenotype-level summary
    ############################################################

    all_sig_phenotypes = all_locus_summary['phenotype'].unique()
    phenotype_summary = pd.DataFrame(index=all_sig_phenotypes)
    phenotype_summary['phenotype'] = phenotype_summary.index
    phenotype_summary['category'] = all_locus_summary.groupby('phenotype')['category'].unique().map(lambda x:x[0])
    phenotype_summary['N_SNP'] = all_locus_summary.groupby("phenotype")['UID'].count()
    phenotype_summary['N_locus'] = all_locus_summary.groupby("phenotype")['UID'].nunique()
    if highlight_marker != 'NA':
        phenotype_summary['N_hight'] = all_locus_summary.groupby("phenotype")[highlight_marker].sum()

    ############################################################
    # 3. Category-level summary
    ############################################################

    # order by mean locus count
    category_order = phenotype_summary.groupby('category')['N_locus'].mean().sort_values()
    category_summary = pd.DataFrame(
        {'order': range(len(category_order))},
        index=category_order.index,
    )
    category_summary['category'] = category_summary.index
    # height of each category
    category_summary['N_pheno'] = phenotype_summary.groupby('category')['phenotype'].count()
    category_summary['N_locus'] = phenotype_summary.groupby('category')['N_locus'].sum()
    # y of the bottom of each category
    category_summary['bottom'] = category_summary['N_pheno'].cumsum()-category_summary['N_pheno']

    # palette_path
    if palette_path in bulitin_palettes.keys():
        if len(category_summary) > bulitin_palettes[palette_path]:
            print(
                f"Palette {palette_path} supports at most "
                f"{bulitin_palettes[palette_path]} categories."
            )
            category_summary['color'] = '#888888'
        else:
            # load cmap from name
            cmap = plt.get_cmap(palette_path)
            # Convert RGB/RGBA values to readable hexadecimal colors. # revert order
            category_summary['color'] = [to_hex(color) for color in list(cmap.colors)[:len(category_summary)]][::-1]
    else:
        # considered as a path to json, and load
        category_color:dict = load_config(palette_path)
        category_summary['color'] = category_summary.index.map(
            lambda color: to_hex(category_color[color]) if color in category_color else '#888888'
        )

    ############################################################
    # 4. x axis, chrom length ; y axis, phenotype order
    ############################################################
    
    # length of chrom
    chr_len = pd.Series(Chr_len,name='len')
    # end x of chrom
    chr_len_pre = chr_len.cumsum() + np.linspace(1,25,25,dtype=int) * Chr_sep
    chr_len_pre[0] = 0
    chr_len_pre = chr_len_pre.sort_index()

    uf_global_pos = np.frompyfunc((lambda chrom,pos:pos+chr_len_pre[chrom-1]),2,1)

    # total_length of x
    max_chr = all_locus_summary['chrom'].max()
    chr_len_total = chr_len_pre[max_chr]

    # phenotype y axis
    phenotype_summary['c1'] = phenotype_summary['category'].map(category_summary['order'])
    phenotype_summary:pd.DataFrame = phenotype_summary.sort_values(['c1','N_locus'])
    ymax = len(phenotype_summary)
    phenotype_summary['y'] = range(ymax)

    ############################################################
    # 5. plot coordinate of locus
    ############################################################

    all_locus_summary['y'] = all_locus_summary['phenotype'].map(phenotype_summary.set_index('phenotype')['y'])
    # global pos
    all_locus_summary['global_pos'] = uf_global_pos(all_locus_summary['chrom'],all_locus_summary['top'])
    # x -> theta
    all_locus_summary['theta'] = all_locus_summary['global_pos']/chr_len_total*Dim
    # color
    all_locus_summary['color'] = all_locus_summary['category'].map(category_summary['color'])

    # pleiotropy gene to mark

    # remove 
    all_locus_summary['locus_top_gene'] = all_locus_summary['top_snp_nearest_gene'].map(strip_arrow)
    # gene-N repeat
    repete_gene = all_locus_summary['locus_top_gene'].value_counts().to_frame()
    # gene-pos
    repete_gene['theta'] = all_locus_summary.dropna(subset='locus_top_gene').groupby('locus_top_gene')['theta'].mean()
    # keep the top N_gene_mark
    repete_gene = repete_gene.iloc[:N_gene_mark].copy()
    # sort by theta
    repete_gene = repete_gene.sort_values('theta')

    # Circumferential histogram statistics
    sections = np.linspace(0,Dim,101)
    bins_count_stats = None
    for i,c in enumerate(all_locus_summary['category'].unique()):
        category_this = all_locus_summary[all_locus_summary['category'] == c]
        bins_count_this = pd.cut(category_this.theta,sections).value_counts()
        bins_count_this.name = c
        if bins_count_stats is None:
            bins_count_stats = pd.concat(
                [pd.Series(0,index=bins_count_this.index,name='Init'), bins_count_this]
                ,axis=1)
        else:
            bins_count_stats = pd.concat([bins_count_stats,bins_count_this],axis=1)
    # Sort in chromosome order
    bins_count_stats = bins_count_stats.sort_index()
    bins_max = bins_count_stats.sum(axis=1).max()

    ############################################################
    # 6. draw plot
    ############################################################

    fig = plt.figure(figsize=(30,30))
    fig.set_facecolor('white')

    ax=plt.subplot(111, polar=True)
    ax.set_rorigin(-ymax / 3)
    ax.axis('off')

    yend = ymax
        
    #bar：endpoint

    for i in range(1, max_chr+1): 
        x_loc = (chr_len_pre[i-1] + chr_len[i]/2)/chr_len_total
        chrom_label_rotate = rotate_text(x_loc)

        # box of chrom; layer5
        ax.bar(x=x_loc * Dim,
            height= yend * 0.05,
            width= chr_len[i]/chr_len_total*Dim,
            bottom = yend*1.025, ec='grey',color='white',lw=3,
            alpha=1,zorder=5) # 10 + 0.2
        # name of chrom; layer6
        chrom = str(i) if i <=22 else "X" if i==23 else "Y" if i==24 else "M" if i==25 else ""
        ax.text(x=(chr_len_pre[i-1]+chr_len[i]/2)/chr_len_total*Dim,
                y=yend*1.05,
                s=chrom,
                fontsize=18, rotation=chrom_label_rotate*Dim/(2*np.pi)*360-90,
                horizontalalignment='center',verticalalignment='center',zorder=6)
        # background of category; layer1

        for cate,values in category_summary.iterrows():
            ax.bar(x=x_loc*Dim,
                height = values.N_pheno,
                width = chr_len[i]/chr_len_total*Dim,
                bottom = values.bottom,
                color=(values['color']),
                alpha=0.2,
                zorder=1
                )
        # Histogram background
        ax.bar(x=x_loc * Dim,
                height = yend/4,
                width = chr_len[i]/chr_len_total*Dim,
                bottom = -yend /4,
                color=(0.92,0.92,0.92,1),
                zorder=1
                )
    # x grid; layer2
    for i in phenotype_summary['y']:
        ax.plot(np.linspace(0,Dim,360), i*np.ones(360), lw=1, color='w', zorder=2)

    bin_ymax = bins_max*1.2
    for i in range(1,100,1):
        if bin_ymax//i < 6:
            break
    yzoom = -yend/4/bin_ymax
    ticks = np.arange(0, bin_ymax, i).astype(int)
    # grid of histgram
    for d,i in enumerate(ticks):
        if i == 0:
            ax.plot(np.linspace(0,Dim,360),i*yzoom*np.ones(360), lw=3, color='w', zorder=4)
        else:
            ax.plot(np.linspace(0,Dim,360),i*yzoom*np.ones(360), lw=1, color='w', zorder=2)
        
        ax.text(
            x=Dim, y=i*yzoom,
            s=i, fontsize=16, zorder=2,
            horizontalalignment='left',verticalalignment='center'
        )
    # layer:4
    # Scale marker size with the radial range: larger markers for compact plots
    # and smaller markers for plots with a larger y range.
    marker_size = 100 - (np.clip(ymax, 30, 100) - 30) * 75 / 70
    marker_lw = np.sqrt(marker_size) / 10
    highlight_size = marker_size * 2
    highlight_lw = marker_lw * 2

    if highlight_marker == 'NA':
        ax.scatter(
            x=all_locus_summary.theta,
            y=(all_locus_summary.y+0.5),
            s=marker_size,
            color=all_locus_summary.color,
            marker='o',
            ec='k',
            lw=marker_lw,
            zorder=4,
            alpha=0.7
        )
    else:
        all_locus_false = all_locus_summary[all_locus_summary[highlight_marker] == False]
        all_locus_true  = all_locus_summary[all_locus_summary[highlight_marker] == True]
        # known
        ax.scatter(
            x=all_locus_false.theta,
            y=(all_locus_false.y+0.5),
            s=marker_size,
            color=all_locus_false.color,
            marker='o',
            ec='k',
            lw=marker_lw,
            zorder=4,
            alpha=0.7,
        )
        # novel
        for key,values in all_locus_true.iterrows():
            ax.scatter(
                x=values.theta,
                y=(values.y+0.5),
                s=highlight_size,
                color=values.color,
                marker=(4,0,values.theta/(2*np.pi)*360),
                ec='k',
                lw=highlight_lw,
                zorder=4,
                alpha=1
            )
        # legend
        ax.legend(
            handles=[
                Line2D([], [], marker='o', linestyle='None', markersize = 10,
                       markerfacecolor='grey', markeredgecolor='k', linewidth = 1,
                       label='False'),
                Line2D([], [], marker=(4, 0, 0), linestyle='None', markersize = 12 ,
                       markerfacecolor='grey', markeredgecolor='k', linewidth = 2,
                       label='True'),
            ],
            loc='upper left',
            frameon=False,
            title=highlight_marker,
            title_fontproperties={'size': 18, 'weight': 'bold'},
            fontsize=18,
        )
    # mark gene; use radial proportions so annotations remain aligned when
    # the number of phenotypes (and therefore yend) changes.

    word_wide = 0.0138 * yend
    last_theta = -1
    for key,values in repete_gene.iterrows():
        # theta project distance 0.02
        if values.theta - last_theta > 0.02:
            this_theta = values.theta
        else:
            this_theta = last_theta + 0.02
        # line layer:3
        ax.plot(
            [values.theta, values.theta, this_theta, this_theta],
            [0, 1.08 * yend, 1.11 * yend, 1.125 * yend],
            color='grey', zorder=3, lw=1,
        )
        # gene marker
        if this_theta < np.pi*0.5 or this_theta > np.pi*1.5:
            rotation = this_theta/(2*np.pi)*360
        else:
            rotation = this_theta/(2*np.pi)*360+180
        ax.text(
            x=this_theta, y=len(key)*word_wide + (1.13)*yend,
            s=key, font={'style':'italic'},
            fontsize=18, rotation=rotation,
            horizontalalignment='center',verticalalignment='center'
        )
        # update
        last_theta = this_theta

    # draw center hist

    for i,c in enumerate(all_locus_summary['category'].unique()):
        x_width = Dim/100
        x_pos = sections[:-1] + x_width/2
        y_bottom = yzoom * bins_count_stats.cumsum(axis=1).iloc[:,i]
        y_height = yzoom * bins_count_stats.loc[:,c]
        c_bar = category_summary.loc[c,'color']
        ax.bar(x=x_pos,width=x_width,bottom=y_bottom,height=y_height,color=c_bar,zorder=3)

    
    # Radial histogram
    ax_c = ax.inset_axes([0.61,0.060,0.35,0.333])
    ax_c.axis('off')
    
    # Radial histogram

    ax_c.barh(
        y=phenotype_summary['y'],
        height=1,
        width=phenotype_summary['N_locus'],
        left=0,
        alpha=0.5,
        color=phenotype_summary['category'].map(category_summary['color'])
    )
    if highlight_marker != 'NA':
        ax_c.barh(
            y=phenotype_summary['y'],
            height=1,
            width=phenotype_summary['N_hight'],
            left=0,
            alpha=1,
            color=phenotype_summary['category'].map(category_summary['color'])
        )
        ax_c.legend(
            handles=[
                Patch(facecolor='grey', edgecolor='none', alpha=0.5, label='False'),
                Patch(facecolor='grey', edgecolor='none', alpha=1, label='True'),
            ],
            loc='upper right',
            frameon=False,
            title=highlight_marker,
            title_fontproperties={'size': 18, 'weight': 'bold'},
            fontsize=18,
        )
    # Category names
    for cate,values in category_summary.iterrows():
        y_this = values.bottom+values.N_pheno/2
        ax_c.text(
            x=-0.5,y=y_this,s=cate,
            color=values.color, font={'weight':'bold'},
            fontsize=18,
            horizontalalignment='right',verticalalignment='center')
    # Axes
    pheno_max = phenotype_summary['N_locus'].max()
    bin_xmax = pheno_max*1.1
    for i in range(1,100,1):
        if bin_xmax//i < 6:
            break
    ticks = np.arange(0, bin_xmax, i).astype(int)
    ax_c.plot([0,0,bin_xmax],[-0.1,yend,yend],color='black',lw=2)
    for i in ticks:
        ax_c.plot([i,i],[yend, yend*1.01],color='black',lw=2)
        ax_c.text(x=i,y=yend*1.02,s=str(i), fontsize=18, horizontalalignment='center',verticalalignment='top')
    ax_c.text(x=bin_xmax/2,y=yend*1.05,s='Number of identified loci', fontsize=18,horizontalalignment='center',verticalalignment='top')
    ax_c.set_ylim(1.02*yend, -1)

    ############################################################
    # 7. save
    ############################################################

    fig.savefig(os.path.join(output_path, 'step6.circos_plot.'+output_format), dpi=200)

    category_summary.to_csv(os.path.join(output_path, 'step6.category_summary.csv'), sep='\t', index=False)

if __name__ == "__main__":
    main(sys.argv)

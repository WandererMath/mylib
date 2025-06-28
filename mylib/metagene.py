import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from .utils import *


def bedgraph2pd(bedgraph):
    col_names = ['chrom', 'start', 'end', 'value']
    df = pd.read_csv(
        bedgraph,
        sep='\t',
        skiprows=1,
        header=None,
        names=col_names
    )

    return df

def parse_gtf(gtf):
    # (chrom, direction) : (start, end)
    chrom_start_end={}
    for gene_id in gtf.all_genes():
        start, stop, direction, chrom= gtf.get_gene_positions_strand_chrom(gene_id)
        start+=1
        # 1-based []
        if not gtf.is_protein(gene_id):
            continue
        if (chrom, direction) not in chrom_start_end:
            chrom_start_end[(chrom, direction)] = []
        chrom_start_end[(chrom, direction)].append((start, stop))
    return chrom_start_end

def _metagene_1_chrom(chrom, df, start_end_list, direction):
    """
    Perform metagene analysis for one chromosome.
    Return:
        dict {offset: intensity}
    """
    result={}
    df.sort_values(by="start", inplace=True)
    starts_bed = df["start"].tolist()
    values_bed = df["value"].tolist()

    if direction=='-':
        pass

    counter=0
    ptr_gtf=0
    result_tmp={}

    for pos, v in zip(starts_bed, values_bed):
        next_gene_flag= False
        try:
            while pos > start_end_list[ptr_gtf][1]:
                ptr_gtf+=1
                next_gene_flag = True
            if next_gene_flag:
                # Save previous gene's results
                n= counter/gene_length *100
                for offset in result_tmp:
                    if counter<100:
                        break
                    if offset not in result:
                        result[offset] = 0
                    result[offset] += result_tmp[offset] / n
                # Clear cache
                result_tmp={}
                counter=0

            target_start, target_end = start_end_list[ptr_gtf]
            if pos >= target_start:
                counter+=v
            offset = pos - target_start
            gene_length = target_end - target_start + 1
            result_tmp[offset]=v
            

        except IndexError:
            break

    return result

def _plot_metagene(result, chrom, direction):
    '''
    Input:
        dict {offset: intensity}
    '''
    WITHIN=100
    OUTSIDE=50
    X=list(range(-OUTSIDE, WITHIN+1))
    Y=[]
    for i in X:
        if i not in result:
            Y.append(0)
        else:
            Y.append(result[i])
    plt.figure(figsize=(10, 4))
    plt.scatter(X, Y, marker='x')
    plt.plot(X, Y)

    ######################################
    # Formatting
    xticks = np.arange((min(X)//3)*3, max(X)+1, 3)
    _xticks = np.arange((min(X)//3+1)*3, max(X)+1, 6)
    plt.xticks(_xticks)

    # Draw vertical dashed lines at x-tick positions
    for xtick in xticks:
        plt.axvline(x=xtick, color='gray', linestyle='--', linewidth=0.7)
    plt.gca().set_xlim([min(X), max(X)])

    plt.title('5\' End')
    plt.savefig(f'{chrom}{direction}.pdf')
    plt.close()

    

def metagene_analysis(gtf, bedgraph):
    print(f"Performing metagene analysis for group {bedgraph}")
    if 'plus' in bedgraph:
        DIRECTION= '+'
    elif 'minus' in bedgraph:
        DIRECTION= '-'
    else:
        raise ValueError("Bedgraph file name must contain 'plus' or 'minus' to indicate strand direction.")
    
    print('Parsing bedgraph file...')
    df= bedgraph2pd(bedgraph)
    df_chrom_grouped = df.groupby('chrom')
    # for chrom, subdf in df_chrom_grouped:
    #     print(f"Chrom: {chrom}")
    #     print(subdf.head())
    
    print('Parsing GTF file...')
    chrom_start_end = parse_gtf(gtf)

    for chrom, subdf in df_chrom_grouped:
        print(f"Metagene Analysis for Chrom: {chrom} Direction: {DIRECTION}")
        result=_metagene_1_chrom(chrom, subdf, chrom_start_end[(chrom, DIRECTION)], DIRECTION)
        _plot_metagene(result, chrom, DIRECTION)

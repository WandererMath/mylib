import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from .utils import *
LARGEST_INT32=2147483647
LARGEST_INT32=0

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
    for feature in gtf.db.features_of_type('CDS'):
        gene_id=feature.id
        start, stop, direction, chrom= gtf.get_gene_positions_strand_chrom(gene_id)
        start+=1
        # 1-based []
        if (chrom, direction) not in chrom_start_end:
            chrom_start_end[(chrom, direction)] = []
        if direction=='+':
            chrom_start_end[(chrom, direction)].append((start, stop))
        elif direction=='-':
            chrom_start_end[(chrom, direction)].append((LARGEST_INT32-stop, LARGEST_INT32-start))
    for key in chrom_start_end:
        chrom_start_end[key].sort(key=lambda x: x[0])
    return chrom_start_end

def _metagene_1_chrom(df, start_end_list, direction, Five_or_Three='5'):
    """
    Perform metagene analysis for one chromosome.
    Return:
        dict {offset: intensity}
    """
    result={}
    df.sort_values(by="start", inplace=True)
    if Five_or_Three == '5':
        Five_Prime = True
    elif Five_or_Three == '3':
        Five_Prime = False
    else:
        raise ValueError("Five_or_Three must be either '5' or '3'")
    starts_bed=df['start'].tolist() 
    values_bed = df["value"].tolist()

    if direction=='-':
        starts_bed = [LARGEST_INT32 - s for s in starts_bed[::-1]]

    counter=0
    ptr_gtf=0
    result_tmp={}

    for pos, v in zip(starts_bed, values_bed):
        next_gene_flag= False
        try:
            while (pos > start_end_list[ptr_gtf][1] and Five_Prime) \
                or (pos > start_end_list[ptr_gtf+1][0] and not Five_Prime):
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
            if (pos >= target_start and Five_Prime)\
                or (pos <= target_end and not Five_Prime):
                counter+=v
            offset = pos - target_start if Five_Prime\
                else pos- target_end
            gene_length = target_end - target_start + 1
            result_tmp[offset]=v
            

        except IndexError:
            break

    return result

def _plot_metagene(result, chrom, direction, Five_or_Three='5'):
    '''
    Input:
        dict {offset: intensity}
    '''
    WITHIN=100
    OUTSIDE=50
    X=list(range(-OUTSIDE, WITHIN+1)) if Five_or_Three == '5' \
        else list(range(-WITHIN, OUTSIDE+1))
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

    plt.title(f'{Five_or_Three}\' End')
    plt.savefig(f'{Five_or_Three}Prime_{chrom}{direction}.pdf')
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
        for Five_or_Three in ['5', '3']:
            print(f"Processing {Five_or_Three} Prime...")
            result=_metagene_1_chrom(subdf, chrom_start_end[(chrom, DIRECTION)], DIRECTION, Five_or_Three)
            _plot_metagene(result, chrom, DIRECTION, Five_or_Three)


def metagene_analysis_all(gtf, folder_path):
    folder_path=os.path.abspath(folder_path)
    samples=get_paths_ends_with_something(folder_path, '.bedgraph')
    output_folder=os.path.join(os.path.dirname(folder_path), 'metagene2')
    os.makedirs(output_folder, exist_ok=True)
    for sample in samples:
        sample_name=os.path.basename(sample).split('.')[0]
        print(f"Processing sample: {sample_name}")
        sample_output_folder= os.path.join(output_folder, sample_name)
        os.makedirs(sample_output_folder, exist_ok=True)
        os.chdir(sample_output_folder)
        metagene_analysis(gtf, sample)

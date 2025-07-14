from .gtf import GTF
from .plotting_tools import *

import os
import pysam
from intervaltree import IntervalTree, Interval

def parse_gtf(gtf: GTF):
    '''
    Return:
        Dict (chrom, direction): list of (start, end)
        Dict (chrom, direction): {start: type (transcript or CDS)}
    '''
    # (chrom, direction) : (start, end)
    chrom_start_end={}
    chrom_start_type={}
    for feature in gtf.db.features_of_type(('CDS', 'transcript')):
        start, stop, direction, chrom= gtf.get_gene_positions_strand_chrom(feature.id)
        # start+=1
        # 1-based []
        if (chrom, direction) not in chrom_start_end:
            chrom_start_end[(chrom, direction)] = []
            chrom_start_type[(chrom, direction)] = {}
        chrom_start_end[(chrom, direction)].append((start, stop))
        chrom_start_type[(chrom, direction)][start]=feature.featuretype

    return chrom_start_end, chrom_start_type


def parse_sam(sam_path):
    '''
    Return:
        Dict (chrom, direction): [start (1-based), length]
    '''
    result={}
    with pysam.AlignmentFile(sam_path, "r") as samfile:
        for read in samfile.fetch():
            if not isinstance(read.reference_length, int):
                continue
            direction='-' if read.is_reverse else '+'
            chrom, start,  length, end = read.reference_name, read.reference_start, read.reference_length, read.reference_end

            start+=1
            if (chrom, direction) not in result:
                result[(chrom, direction)]=[]
            result[(chrom, direction)].append((start, end,  length))
    return result

def match2start(tree, x):
    matches = tree[x]
    if matches:
        for m in matches:
            return m.begin
    else:
        return None


def process_sam(sam_path, gtf):
    chrom_start_end, chrom_start_type=parse_gtf(gtf)
    interval_trees={k: IntervalTree(Interval(start, end) for start, end in chrom_start_end[k]) for k in chrom_start_end}
    print("GTF parsed and interval tree built.\n")

    
    sam_data=parse_sam(sam_path)
    print(f'SAM File {sam_path} parsed.')
    protein_coding=[]
    rna=[]
    others=[]
    for k in sam_data:
        for read_start, read_end,  length in sam_data[k]:
            matched_start=match2start(interval_trees[k], read_start)
            if matched_start is None:
                matched_start=match2start(interval_trees[k], read_end)
            if matched_start is None:
                others.append(length)
            elif chrom_start_type[k][matched_start]=='CDS':
                protein_coding.append(length)
            elif chrom_start_type[k][matched_start]=='transcript':
                rna.append(length)
            else:
                raise Exception(f'feature type {chrom_start_type[k][matched_start]} not supported')
    stacked_hist(['protein_coding', 'rRNA, tRNA, etc.', 'unmapped'][::-1], [protein_coding, rna, others][::-1], f"{os.path.basename(sam_path).split('.')[0]}.pdf")
    return protein_coding, rna, others

    
import os
import matplotlib.pyplot as plt
import pandas as pd

from .gtf import GTF
from .utils import *
from .utils import _auto_grouping
# Use Python 3.12. Pysam doesn't support 3.13 yet. 
# Lazy import 
# import pysam

# 
# Indexing rules:
#     SAM: 1-based  [ and length
#     However, pysam API gives 0-based [)
#     GTF: 1-based []
#     Bedgraph: 0-based [)

# 
# Bedgraph example, separated by tabs:
#  track type=bedGraph
#  *	6	7	701
#  *	7	8	20061
#  NC_003063.2	1763	1764	1
#  NC_003063.2	1772	1773	1
#  NC_003063.2	1810	1811	1


class Feature:
    '''
        Class to handle 1 featureCount file
    '''
    name: str
    counts: list    # int
    gene_ids: list # string
    gene_count: dict # gene_id: count
    def __init__(self, path, ribo_seq=False):
        df = pd.read_csv(path, sep="\t", comment="#")
        FIELD_COUNTS=df.columns.tolist()[-1]
        self.counts=df[FIELD_COUNTS].tolist()
        self.gene_ids=df['Geneid'].tolist()
        self.gene_count={}
        for g, c in zip( self.gene_ids, self.counts):
            self.gene_count[g]=c
        self.name=path.split('/')[-1].split('.')[0]
        if ribo_seq:
            self.name='Ribo-'+self.name


def sam2bedgraph(sam_path, ribo=True, offset=14, cutoff=None, output_folder=None):
    import pysam
    sam_path=os.path.abspath(sam_path)
    FIRST_LINE='track type=bedGraph\n'
    BASENAME=os.path.basename(sam_path).split('.')[0]
    DIR_NAME=os.path.join( os.path.dirname(sam_path),\
                          '..',\
                        f"coverage" if output_folder is None else output_folder)
    os.makedirs(DIR_NAME, exist_ok=True)
    prefix='Ribo-' if ribo else 'RNA-'
    OUTPUT_PLUS=os.path.join(DIR_NAME, prefix+BASENAME + '-plus.bedgraph')
    OUTPUT_MINUS=os.path.join(DIR_NAME, prefix+BASENAME + '-minus.bedgraph')

    DIR_NAME_NORM=os.path.join( os.path.dirname(sam_path),\
                          '..',\
                        f"coverage_norm" if output_folder is None else f"{output_folder}_norm_merged")
    os.makedirs(DIR_NAME_NORM, exist_ok=True)
    OUTPUT_PLUS_NORM=os.path.join(DIR_NAME_NORM, prefix+BASENAME + '-plus.bedgraph')
    OUTPUT_MINUS_NORM=os.path.join(DIR_NAME_NORM, prefix+BASENAME + '-minus.bedgraph')

    result={'+':{}, '-':{}} # chrom: dict
                #     # dict :=   start: value
    with pysam.AlignmentFile(sam_path, "r") as samfile:
        for read in samfile.fetch():
            if cutoff is not None:
                if not isinstance(read.reference_length, int):
                    continue
                if read.reference_length < cutoff:
                    continue
            chrom, start, end = read.reference_name, read.reference_start, read.reference_end
            
            strand = '-' if read.is_reverse else '+'   
            if read.is_unmapped:
                continue
            if chrom not in result[strand]:
                result[strand][chrom] = {} 
                
            if ribo:
                start+=1
                if strand=='+':
                    pos= end-offset
                else:
                    pos= start+offset 
            else:
                pos= (start + end)//2

            if pos not in result[strand][chrom]:
                result[strand][chrom][pos] = 1
            else:
                result[strand][chrom][pos] += 1

    total=0
    with open(OUTPUT_PLUS, 'w') as f:
        f.write(FIRST_LINE)
        for chrom in result['+']:
            for pos in result['+'][chrom]:
                total+=result['+'][chrom][pos]
                f.write(f"{chrom}\t{pos}\t{pos+1}\t{result['+'][chrom][pos]}\n")

    with open(OUTPUT_MINUS, 'w') as f:
        f.write(FIRST_LINE)
        for chrom in result['-']:
            for pos in result['-'][chrom]:
                total+=result['-'][chrom][pos]
                f.write(f"{chrom}\t{pos}\t{pos+1}\t{result['-'][chrom][pos]}\n")
    
    # Normalize 
    total/=(1E6)
    with open(OUTPUT_PLUS_NORM, 'w') as f:
        f.write(FIRST_LINE)
        for chrom in result['+']:
            for pos in result['+'][chrom]:
                result['+'][chrom][pos]/=total
                f.write(f"{chrom}\t{pos}\t{pos+1}\t{result['+'][chrom][pos]}\n")

    with open(OUTPUT_MINUS_NORM, 'w') as f:
        f.write(FIRST_LINE)
        for chrom in result['-']:
            for pos in result['-'][chrom]:
                result['-'][chrom][pos]/=total
                f.write(f"{chrom}\t{pos}\t{pos+1}\t{result['-'][chrom][pos]}\n")
    return result



def bedgraph_for_all_samples(path, ribo=True, gtf: GTF=None, cutoff=None, output_folder='coverage_norm_merged'):
    path=os.path.abspath(path)
    DIR_NAME=os.path.dirname(path)
    OUTPUT_DIR=os.path.join(DIR_NAME, output_folder)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    
    paths= get_paths_ends_with_something(path, '.txt')
    basenames=[os.path.basename(p).split('.')[0] for p in paths]
    grouped= _auto_grouping(basenames)
    for group in grouped:
        samples= grouped[group]
        num_samples=len(samples)
        print(f"Processing group: {group}")
        results=[]
        merged_result={}
        # Merge bedgraphs
        for sample in samples:
            sam_path=paths[basenames.index(sample)]
            # Normalized bedgraph data
            results.append(sam2bedgraph(sam_path, ribo=ribo, cutoff=cutoff, output_folder=output_folder.replace('_merged', '') ))

        for result in results:
            for strand in result:
                if strand not in merged_result:
                    merged_result[strand] = {}
                for chrom in result[strand]:
                    if chrom not in merged_result[strand]:
                        merged_result[strand][chrom] = {}
                    for pos in result[strand][chrom]:
                        if pos not in merged_result[strand][chrom]:
                            merged_result[strand][chrom][pos] = result[strand][chrom][pos]
                        else:
                            merged_result[strand][chrom][pos] += result[strand][chrom][pos]
        # Write merged results
        for strand in merged_result:
            output_file = os.path.join(OUTPUT_DIR, f"{'Ribo-' if ribo else 'RNA-'}{group}{strand}.bedgraph")
            with open(output_file, 'w') as f:
                f.write('track type=bedGraph\n')
                for chrom in merged_result[strand]:
                    for pos in merged_result[strand][chrom]:
                        f.write(f"{chrom}\t{pos}\t{pos+1}\t{merged_result[strand][chrom][pos]/num_samples}\n")


def sam2distribution(path):
    import pysam
    path=os.path.abspath(path)

    basename=os.path.basename(path).split('.')[0]
    output_dir=os.path.join( os.path.dirname(path),'..', "length_distribution")
    os.makedirs(output_dir, exist_ok=True)
    output_file=os.path.join(output_dir, basename + '.pdf')

    lengths=[]
    with pysam.AlignmentFile(path, "r") as samfile:
        for read in samfile.fetch():
            if isinstance(read.reference_length, int):
                lengths.append(read.reference_length)
    plt.hist(lengths, bins=max(lengths)-min(lengths)+1)
    plt.title(f"Length Distribution of {basename}")
    plt.savefig(output_file)
    plt.close()
    print(f"Length distribution saved to {basename}")

def sam2distribution_all_files(path):
    for p in get_paths_ends_with_something(path, '.txt'):
        sam2distribution(p)

def filter_1_sam(path, threshold=23, high_threshold=None):
    import pysam

    # Input/output files
    path = os.path.abspath(path)
    if high_threshold is None:
        OUTPUT_DIR= os.path.join(os.path.dirname(path),'..', f'bowtie{threshold}')
    else:
        OUTPUT_DIR= os.path.join(os.path.dirname(path),'..', f'bowtie_{threshold}_{high_threshold}')
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    output_sam = os.path.join(OUTPUT_DIR, os.path.basename(path))

    # Open SAM input and output
    with pysam.AlignmentFile(path, "r") as infile, \
        pysam.AlignmentFile(output_sam, "w", template=infile) as outfile:

        if high_threshold is None:
            for read in infile:
                if not read.is_unmapped and read.query_length >= threshold:
                    outfile.write(read)
        else:
            for read in infile:
                if not read.is_unmapped and read.query_length >= threshold\
                    and read.query_length<=high_threshold:
                    outfile.write(read)         


def filter_sams(path, t, ht=None):
    for p in get_paths_ends_with_something(path, '.txt'):
        filter_1_sam(p, t, ht)

if __name__ == "__main__":
    print("This is the mylib/bedgraph.py module.")
    SAM_FOLDER='/fs/ess/PAS2967/S21/RNA-seq/bowtie23'
    FNA='/fs/ess/PAS2967/S21/reference/GCF_000092025.1_ASM9202v1_genomic.fna'
    GTF_FILE='/fs/ess/PAS2967/S21/reference/genomic.gtf'
    gtf= GTF(GTF_FILE, FNA)
    bedgraph_for_all_samples(SAM_FOLDER, t=23)
    #bedgraph_for_all_samples(SAM_FOLDER, ribo=True, gtf=gtf, cutoff=23)

    #gtf=GTF('/users/PAS0291/dengyw144/Fredrick008/ref/genomic.gtf', '/users/PAS0291/dengyw144/Fredrick008/ref/GCF_000750555.1_ASM75055v1_genomic.fna')
    #sam2bedgraph(SAM)
    #bedgraph_for_all_samples('/fs/ess/PAS2967/S21/Ribo-seq/bowtie_renamed/', gtf=gtf)
    #bedgraph_for_all_samples('/users/PAS0291/dengyw144/Fredrick008/data08/bowtie', gtf=gtf)

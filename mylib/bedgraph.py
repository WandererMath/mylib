import os

from gtf import GTF
from utils import *
# Use Python 3.12. Pysam doesn't support 3.13 yet. 
# Lazy import 
# import pysam

'''
Indexing rules:
    SAM: 1-based  [ and length
    However, pysam API gives 0-based [)
    GTF: 1-based []
    Bedgraph: 0-based [)

'''
# Bedgraph example, separated by tabs:
#  track type=bedGraph
#  *	6	7	701
#  *	7	8	20061
#  NC_003063.2	1763	1764	1
#  NC_003063.2	1772	1773	1
#  NC_003063.2	1810	1811	1


def sam2bedgraph(sam_path, ribo=True, offset=15):
    import pysam
    sam_path=os.path.abspath(sam_path)
    FIRST_LINE='track type=bedGraph\n'
    BASENAME=os.path.basename(sam_path).split('.')[0]
    DIR_NAME=os.path.join( os.path.dirname(sam_path),\
                          '..',\
                        "coverage")
    os.makedirs(DIR_NAME, exist_ok=True)
    prefix='Ribo-' if ribo else 'RNA-'
    OUTPUT_PLUS=os.path.join(DIR_NAME, prefix+BASENAME + '-plus.bedgraph')
    OUTPUT_MINUS=os.path.join(DIR_NAME, prefix+BASENAME + '-minus.bedgraph')

    DIR_NAME_NORM=os.path.join( os.path.dirname(sam_path),\
                          '..',\
                        "coverage_norm")
    os.makedirs(DIR_NAME_NORM, exist_ok=True)
    OUTPUT_PLUS_NORM=os.path.join(DIR_NAME_NORM, prefix+BASENAME + '-plus.bedgraph')
    OUTPUT_MINUS_NORM=os.path.join(DIR_NAME_NORM, prefix+BASENAME + '-minus.bedgraph')

    result={'+':{}, '-':{}} # chrom: dict
                #     # dict :=   start: value
    with pysam.AlignmentFile(sam_path, "r") as samfile:
        for read in samfile.fetch():
            chrom, start, end = read.reference_name, read.reference_start, read.reference_end
            
            strand = '-' if read.is_reverse else '+'   
            if read.is_unmapped:
                continue
            if chrom not in result[strand]:
                result[strand][chrom] = {} 
                
            if ribo:
                if strand=='+':
                    pos= end-offset -1
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


def _auto_grouping(samples):
    result={}
    for s in samples:
        if s[0] not in result:
            result[s[0]]=[s]
        else:
            result[s[0]].append(s)
    return result

def bedgraph_for_all_samples(path, ribo=True):
    path=os.path.abspath(path)
    DIR_NAME=os.path.dirname(path)
    OUTPUT_DIR=os.path.join(DIR_NAME, 'coverage_norm_merged')
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    paths= get_paths_ends_with_something(path, '.txt')
    basenames=[os.path.basename(p).split('.')[0] for p in paths]
    grouped= _auto_grouping(basenames)
    for group in grouped:
        samples= grouped[group]
        print(f"Processing group: {group}")
        results=[]
        merged_result={}
        # Merge bedgraphs
        for sample in samples:
            sam_path=paths[basenames.index(sample)]
            # Normalized bedgraph data
            results.append(sam2bedgraph(sam_path, ribo=ribo))
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
            output_file = os.path.join(OUTPUT_DIR, f"{'ribo-' if ribo else 'RNA-'}{group}{strand}.bedgraph")
            with open(output_file, 'w') as f:
                f.write('track type=bedGraph\n')
                for chrom in merged_result[strand]:
                    for pos in merged_result[strand][chrom]:
                        f.write(f"{chrom}\t{pos}\t{pos+1}\t{merged_result[strand][chrom][pos]}\n")


if __name__ == "__main__":
    print("This is the mylib/bedgraph.py module.")
    SAM='/fs/ess/PAS2967/S21/Ribo-seq/bowtie_renamed/A1.txt'
    FNA='/fs/ess/PAS2967/S21/reference/GCF_000092025.1_ASM9202v1_genomic.fna'
    GTF_FILE='/fs/ess/PAS2967/S21/reference/genomic.gtf'
    gtf= GTF(GTF_FILE, FNA)
    #sam2bedgraph(SAM)
    bedgraph_for_all_samples('/fs/ess/PAS2967/S21/Ribo-seq/bowtie_renamed/')

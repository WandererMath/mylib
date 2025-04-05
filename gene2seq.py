from Bio import SeqIO
from Bio.Seq import Seq
import gffutils

import os

def get_start_and_direction(gene_id, FILE_GTF):
    gtf_path_splitted=FILE_GTF.split('.')
    gtf_path_splitted[-1]='db'
    db_path=".".join(gtf_path_splitted)
    if os.path.exists(db_path):
        db = gffutils.FeatureDB(db_path, keep_order=True)
    else:
        db = gffutils.create_db(FILE_GTF, dbfn=db_path, force=True, keep_order=True, merge_strategy='merge')

    # Retrieve the gene entry by its ID
    try:
        gene = db[gene_id]
        start_position= gene.start
        return start_position,db[gene_id].strand
        #print(f"The start position of gene {gene_id} is {start_position}")
    except KeyError:
        print(f"Gene ID {gene_id} not found in the GTF file.")

def get_seq(start, end, strand, FILE_FNA):
    with open(FILE_FNA, 'r') as fna_file:
        for record in SeqIO.parse(fna_file, 'fasta'):
            # Get the sequence from start to end position
            sequence = record.seq[start - 1:end]  # Adjust to 0-based indexing
            #print(f"Extracted sequence from position {start} to {end}:")
            #print(sequence)
            if strand=="-":
                s=Seq(sequence)
                s=s.reverse_complement()
                sequence=str(s)
            return sequence
        
def get_seq_from_gene_id(gene_id, offset1, offset2, fna_path, gtf_path):
    # offset2> offest1
    start, direction=get_start_and_direction(gene_id, gtf_path)
    return get_seq(start-offset2, start-offset1, direction,fna_path)

if __name__=='__main__':
    FILE_FNA='GCF_000750555.1_ASM75055v1_genomic.fna'
    FILE_GTF='genomic.gtf'
    seq=get_seq_from_gene_id("BW25113_RS00035", 1, 50, FILE_FNA, FILE_GTF)
    print(seq)
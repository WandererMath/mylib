from Bio import SeqIO
from Bio.Seq import Seq
import gffutils

import os

class GTF:
    def __init__(self, gtf_file):
        gtf_path_splitted=gtf_file.split('.')
        gtf_path_splitted[-1]='db'
        db_path=".".join(gtf_path_splitted)
        if os.path.exists(db_path):
            self.db = gffutils.FeatureDB(db_path, keep_order=True)
        else:
            self.db = gffutils.create_db(gtf_file, dbfn=db_path, force=True, keep_order=True, merge_strategy='merge')
    def all_genes(self):
        """
        Return:
            All gene IDs
        """
        gene_symbols=[]
        for feature in self.db.features_of_type('gene'):
            gene_symbol = feature.attributes.get('gene_id', [None])[0]  # Using [None] as fallback in case the attribute is missing
            if gene_symbol:
                gene_symbols.append(gene_symbol)
        return gene_symbols
    def id2name(self, gene_id):
        gene_feature = self.db[gene_id]  # Lookup gene feature by ID
        # Common keys for gene symbol in GTF files: 'gene_name', 'gene_symbol'
        return gene_feature.attributes.get('gene', [None])[0] 
    
    @staticmethod
    def get_seq(start, end, strand, FILE_FNA):
        with open(FILE_FNA, 'r') as fna_file:
            for record in SeqIO.parse(fna_file, 'fasta'):
                seq=record.seq[start - 1:end]
                if strand=='+':
                    return seq # Adjust to 0-based indexing
                elif strand=="-":
                    #print('minus')
                    s=Seq(seq)
                    s=s.reverse_complement()
                    return str(s)
                else:
                    raise Exception("strand must be \"+\" or \"-\" ")

    def get_start_and_direction(self, gene_id):

        # Retrieve the gene entry by its ID
        try:
            gene = self.db[gene_id]
            start_position= gene.start
            return start_position, self.db[gene_id].strand
            #print(f"The start position of gene {gene_id} is {start_position}")
        except KeyError:
            print(f"Gene ID {gene_id} not found in the GTF file.")


    def get_seq_from_gene_id(self, gene_id, offset1, offset2, fna_path):
        '''
        offset2> offest1
        '''
        start, direction=self.get_start_and_direction(gene_id)
        if direction=='+':
            return self.get_seq(start-offset2, start-offset1, direction,fna_path)
        elif direction=='-':
            return self.get_seq(start+offset1, start+offset2, direction,fna_path)



######## Legacy

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
        return start_position, db[gene_id].strand
        #print(f"The start position of gene {gene_id} is {start_position}")
    except KeyError:
        print(f"Gene ID {gene_id} not found in the GTF file.")
    



### !!!! Problematic for minus strand
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
    '''
      offset2> offest1
    '''
    start, direction=get_start_and_direction(gene_id, gtf_path)
    return get_seq(start-offset2, start-offset1, direction,fna_path)

if __name__=='__main__':
    FILE_FNA='GCF_000007825.1_ASM782v1_genomic.fna'
    FILE_GTF='genomic.gtf'
    #seq=get_seq_from_gene_id("BW25113_RS00035", 1, 50, FILE_FNA, FILE_GTF)
    #print(seq)
    gtf=GTF(FILE_GTF)
    #all_genes=gtf.all_genes()
    #print(gtf.id2name(all_genes[0]))
    seq=gtf.get_seq_from_gene_id("BC_RS21510", 1, 50, FILE_FNA)
    print(seq)
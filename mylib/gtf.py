import os
from dataclasses import dataclass

from Bio import SeqIO
from Bio.Seq import Seq
import gffutils
import gffutils.attributes
import gffutils.feature


#gffutils.attributes.Attributes
#gffutils.FeatureDB.execute

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
    
    def get_info_from_id(self, gene_id):
        '''
            Return:
                A dict
        '''
        gene_feature = self.db[gene_id] 
        gene_feature: gffutils.feature.Feature
        #breakpoint()
        return gene_feature.attributes.__dict__['_d']

    def save_IDs_info(self, IDs, output, fields=['gene_id', 'gene', 'gene_biotype', 'product', 'go_function']):
        def get_value(id, key):
            info=self.get_info_from_id(id)
            try:
                # value is a list
                value=info[key]
                if len(value)==0:
                    return ''
                else:
                    return ';'.join(v for v in value)
            except KeyError:
                return ''
        with open(output, 'w') as f:
            f.write(','.join([k for k in fields])+'\n')
            for id in IDs:
                line=','.join([get_value(id, key) for key in fields])
                f.write(line+'\n')

    @staticmethod
    def get_seq(start, end, strand, FILE_FNA):
        with open(FILE_FNA, 'r') as fna_file:
            for record in SeqIO.parse(fna_file, 'fasta'):
                seq=record.seq[start - 1:end]
                s=Seq(seq)
                if strand=='+':
                    return s.transcribe() 
                elif strand=="-":
                    #print('minus')
                    s=s.reverse_complement_rna()
                    return str(s)
                else:
                    raise Exception("strand must be \"+\" or \"-\" ")

    def get_start_and_direction(self, gene_id):

        # Retrieve the gene entry by its ID
        try:
            gene = self.db[gene_id]
            start_position= gene.start
            strand=self.db[gene_id].strand
            if strand=="+":
                return start_position, strand
            else:
                return gene.end, strand
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

    def name2id(self, name):
        for id in self.all_genes():
            if self.id2name(id)==name:
                return id


'''
    ASD
L   
O   
N   
G   
'''


@dataclass
class Aligner:
    ASD: str
    pairing_rule: list # list of sets
    def _nt_pairs(self, a, b)->bool:
        if any(set([a, b])==rule for rule in self.pairing_rule):
            return True
        return False
    def _pair1(self, seq, n):
        # align at n-th of the long 
        # return staring index (of the long) of the longest pairing region and the pairing length
        tmp={}
        j=n
        for i in range(n, len(self.ASD)+n):
            if self._nt_pairs(seq[i],self.ASD[i-n]):
                if j in tmp:
                    tmp[j]+=1
                else:
                    tmp[j]=1
            else:
                j=n+1
        max_key = max(tmp, key=tmp.get)
        max_value = tmp[max_key]
        return max_key, max_value
    
    def pair(self, seq):
        tmp=[]
        for i in range(len(seq)-len(self.ASD)+1):
            r=self._pair1(seq, i)
            tmp.append(r)
        return max(tmp, key=lambda x: x[1])
    def output_pairing(self, seq):
        i, l=self.pair(seq)
        return f"{seq[:i]}[{seq[i:i+l]}]{seq[i+l:]}", i, l
            
asd_ba_sub=Aligner('UCACCUCCUUUCU', [{'A', 'U'},{'C', 'G'}, {'G', "U"}])

if __name__=='__main__':
    FILE_FNA='GCF_000007825.1_ASM782v1_genomic.fna'
    FILE_GTF='genomic.gtf'

    print(asd_ba_sub.output_pairing("CCCCUUUUCCC"))
    
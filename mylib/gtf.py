import os
from dataclasses import dataclass

from Bio import SeqIO
from Bio.Seq import Seq
import gffutils
import gffutils.attributes
import gffutils.feature


START_CODONS=['AUG', 'GUG', 'UUG'] # AUG is the only one that codes for Methionine

#gffutils.attributes.Attributes
#gffutils.FeatureDB.execute

class GTF:
    db: gffutils.FeatureDB
    fna: Seq
    def __init__(self, gtf_file, fna_file=None):
        gtf_path_splitted=gtf_file.split('.')
        gtf_path_splitted[-1]='db'
        db_path=".".join(gtf_path_splitted)
        if os.path.exists(db_path):
            self.db = gffutils.FeatureDB(db_path, keep_order=True)
        else:
            self.db = gffutils.create_db(gtf_file, dbfn=db_path, force=True, keep_order=True, merge_strategy='merge')
    
        if fna_file is not None:
            with open(fna_file, 'r') as fna_file:
                for record in SeqIO.parse(fna_file, 'fasta'):
                    self.fna=Seq(record.seq.transcribe())  # Transcribe the sequence to RNA
                    break

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
    
    def all_start_positions(self):
        '''
            Return:
                {'+':[int, ...], '-':[]}
        '''
        results = {'+': [], '-': []}
        for gene_id in self.all_genes():
            start, strand = self.get_start_and_direction(gene_id)
            if strand == '+':
                results['+'].append(start)
            elif strand == '-':
                results['-'].append(start)
        return results

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
    
    def save_IDs_info(self, IDs, output, fields=['gene_id', 'gene', 'gene_biotype', 'product', 'go_function', 'go_process']):        
        to_write=[]
        for feature in self.db.features_of_type('CDS'):
            gene_id = feature.attributes.get('gene_id', [''])[0]  # Using [None] as fallback in case the attribute is missing
            if gene_id in IDs:
                values=[]
                for f in fields:
                    values_list_serialized=';'.join(feature.attributes.get(f, ['']))
                    values.append(values_list_serialized)
                to_write.append(','.join([a for a in values]))
        with open(output, 'w') as f:
            f.write(','.join([k for k in fields])+'\n')
            for row in to_write:
                f.write(row+'\n')

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
    
    def get_low_high_direction(self, gene_id):
        try:
            gene = self.db[gene_id]
            return gene.start, gene.end, gene.strand
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

    def get_gene_whole_seq(self, id, fna, extend=0):
        # For gtf file, start end positions are inclusive
        l, h, s= self.get_low_high_direction(id)
        if s=='+':
            h+=extend
        else:
            l-=extend
        with open(fna, 'r') as fna_file:
            for record in SeqIO.parse(fna_file, 'fasta'):
                seq=record.seq[l:h+1]
                seq=Seq(seq)
        if s=="+":
            return seq.transcribe()
        else:
            return seq.reverse_complement_rna()

    def name2id(self, name):
        for id in self.all_genes():
            if self.id2name(id)==name:
                return id
    
    def get_all_candidate_start_positions_plus(self):
        return [i for i in range(len(self.fna) - 2)
           if self.fna[i:i+3] in START_CODONS]
    
    def get_all_candidate_start_positions_minus(self):
        START_CODONS_REVERSE_COMPLEMENT = [Seq(codon).reverse_complement_rna() for codon in START_CODONS]
        return [i+3 for i in range(len(self.fna) - 2)
           if self.fna[i:i+3] in START_CODONS_REVERSE_COMPLEMENT]

    def get_tir_by_position(self, position, strand, l=-18, h=9):
        """
        Get the TIR sequence around a given position.
        Args:
            position (int): The position in the sequence.
            l (int): Length to extend to the left.
            h (int): Length to extend to the right.
        Returns:
            str: The TIR sequence.
            None: if not enough sequence available.
        """
        if strand == '+':
            start = max(0, position + l)
            end = min(len(self.fna), position + h)
            if end-start!=h-l:
                return None
            return str(self.fna[start:end])
        elif strand == '-':
            start = max(0, position - h)
            end = min(len(self.fna), position - l)
            if end-start!=h-l:
                return None
            seq = self.fna[start:end]
            if any(x in seq for x in ['Y', 'N', 'R']):
                return None
            return str(seq.reverse_complement_rna())
    def prepare_data(self):
        """
        Return:
            [seqs], [labels]
            label: 1 Gene, 0 Not Gene
        """
        positions_plus=self.get_all_candidate_start_positions_plus()
        positions_minus=self.get_all_candidate_start_positions_minus()
        tmp=self.all_start_positions()
        gene_positions_plus=tmp['+']
        gene_positions_minus=tmp['-']
        seqs=[self.get_tir_by_position(p, '+') for p in positions_plus if self.get_tir_by_position(p, '+') is not None]\
            +[self.get_tir_by_position(p, '-') for p in positions_minus if self.get_tir_by_position(p, '+') is not None]
        labels=[1 if p in gene_positions_plus else 0 for p in positions_plus if self.get_tir_by_position(p, '+') is not None]\
            +[1 if p in gene_positions_minus else 0 for p in positions_minus if self.get_tir_by_position(p, '+') is not None]
        assert len(labels)==len(seqs), f"Labels and sequences length mismatch: {len(labels)} != {len(seqs)}"
        return seqs, labels

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
                j=i+1
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
            

if __name__=='__main__':
    FILE_FNA='GCF_000007825.1_ASM782v1_genomic.fna'
    FILE_GTF='genomic.gtf'
    gtf=GTF(FILE_GTF, FILE_FNA)
    #print(gtf.fna[:100])
    #print(gtf.all_start_positions()['-'][-10:])
    #for p in gtf.get_all_candidate_start_positions_plus()[:100]:
    #    print(gtf.get_tir_by_position(p, '+'))
    seqs, labels=gtf.prepare_data()
    breakpoint()
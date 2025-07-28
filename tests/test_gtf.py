import unittest
import random

from mylib.gtf import GTF
from mylib.bedgraph import filter_sams, bedgraph_for_all_samples
from mylib.metagene import metagene_analysis_all

class TestGTF(unittest.TestCase):
    def setUp(self):
        self.gtf = GTF('/fs/ess/PAS2967/mylib/data/genomic.gtf', "/fs/ess/PAS2967/mylib/data/GCF_000092025.1_ASM9202v1_genomic.fna")
    def test_get_seq(self):
        genes=self.gtf.all_genes()
        gene= random.choice(genes)
        seq= self.gtf.get_seq_by_gene_and_offset(gene, -3, 3)
        print(seq)

    def test_get_gene_info(self):
        self.gtf.save_IDs_info(self.gtf.all_genes(), 'all_genes.csv')
    
    def test_get_info_from_id(self):
        gene = random.choice(self.gtf.all_genes())
        # gene="ATU_RS02540"
        info = self.gtf.get_info_from_id(gene)
        print(info)
        
    def test_id2protein(self):
        print(self.gtf.id2protein("ATU_RS02540", 'gene'))
        
    def test_prep_data(self):
        seqs, labels=self.gtf.prepare_data()
        breakpoint()
    
    def test_prep_data_encoded(self):
        seqs, labels=self.gtf.prep_data_encoded()
        breakpoint()
    
    def test_transcript(self):
        for transcript in self.gtf.db.features_of_type('transcript'):
            info=self.gtf.get_info_from_id(transcript.id)
            print(info['transcript_biotype'][0])

    def test_ml(self):
        self.gtf.prepare_data()
    
    def test_sam(self):
        SAM_FOLDER='/fs/ess/PAS2967/S21/Ribo-seq/bowtie_gene'
        filter_sams(SAM_FOLDER, 23, 32)
    
    def test_bedgraph(self):
        bedgraph_for_all_samples('/fs/ess/PAS2967/S21/Ribo-seq/bowtie_23_32', ribo=True, output_folder='coverge_norm_merged_23_32')
    def test_metagene(self):
        metagene_analysis_all(self.gtf, \
            '/fs/ess/PAS2967/S21/Ribo-seq/coverge_23_32',\
                'metagene_23_32')
import unittest
import random

from mylib.gtf import GTF

class TestSeq(unittest.TestCase):
    def setUp(self):
        self.gtf = GTF('/fs/ess/PAS2967/mylib/data/genomic.gtf', "/fs/ess/PAS2967/mylib/data/GCF_000092025.1_ASM9202v1_genomic.fna")
    def test_get_seq(self):
        genes=self.gtf.all_genes()
        gene= random.choice(genes)
        seq= self.gtf.get_seq_by_gene_and_offset(gene, -3, 3)
        print(seq)

    def test_get_gene_info(self):
        self.gtf.save_IDs_info(self.gtf.all_genes(), 'all_genes.csv')

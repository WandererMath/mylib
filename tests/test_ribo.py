import unittest
import random

from mylib.gtf import GTF
from mylib.ribo_length_dist import *
from mylib.plotting_tools import *
from mylib.metagene import *
from mylib.create_rrna_ref import *
class TestGTF(unittest.TestCase):
    def setUp(self):
        self.project_base_path=f"{os.environ.get('HOMEP')}/S21"
        self.bedgraph_base_path=os.path.join(self.project_base_path, 'Ribo-seq/coverage')
        self.sam_path='data/W2.txt'
        self.gtf = GTF('/fs/ess/PAS2967/mylib/data/genomic.gtf', "/fs/ess/PAS2967/mylib/data/GCF_000092025.1_ASM9202v1_genomic.fna")
    def test_parse_gtf(self):
        parse_gtf(self.gtf)
    def test_process_sam(self):
        # process_sam(self.sam_path, self.gtf)
        plot_stacked_hist_all_separate_files(os.path.join(self.project_base_path, 'Ribo-seq/bowtie_rrna'),\
            os.path.join(self.project_base_path, 'Ribo-seq/bowtie_gene'),\
            self.gtf)
    def test_metagene(self):
        # metagene_analysis(GTF('/users/PAS0291/dengyw144/Fredrick008/ref/genomic.gtf','/users/PAS0291/dengyw144/Fredrick008/ref/GCF_000750555.1_ASM75055v1_genomic.fna'), '/users/PAS0291/dengyw144/Fredrick008/data08/coverage/Ribo-A01-bowtie1-minus.bedgraph')
        metagene_analysis_all(self.gtf, self.bedgraph_base_path)
    def test_rrna_ref(self):
        create_rrna_ref(self.gtf)
        
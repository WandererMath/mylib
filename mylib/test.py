from .gtf import GTF
from .bedgraph import bedgraph_for_all_samples

import numpy as np

FNA='/fs/ess/PAS2967/S21/reference/GCF_000092025.1_ASM9202v1_genomic.fna'
GTF_FILE='/fs/ess/PAS2967/S21/reference/genomic.gtf'
gtf= GTF(GTF_FILE, FNA)
bedgraph_for_all_samples('/fs/ess/PAS2967/S21/Ribo-seq/bowtie', lt=True, cutoff=18, output_folder='coverage_norm_merged_lt_18')
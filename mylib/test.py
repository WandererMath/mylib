from .gtf import GTF
from .bedgraph import bedgraph_for_all_samples
from .sd import run_fs
import numpy as np
import RNA, ray, multiprocessing

FNA='/fs/ess/PAS2967//dengyw144/S21/reference/GCF_000092025.1_ASM9202v1_genomic.fna'
GTF_FILE='/fs/ess/PAS2967/dengyw144/S21/reference/genomic.gtf'
gtf= GTF(GTF_FILE, FNA)
# gene_16S='ATU_RS00260'
# gene=gtf.all_genes()[2]
# print(gtf.is_protein(gene))
# # seq=gtf.get_gene_seq(gene)
# start, _, strand, chrom=gtf.get_gene_positions_strand_chrom(gene)
# seq=gtf.get_tir_by_position(start,chrom,  strand, -150,150)
# structure, mfe = RNA.fold(seq)  # predict MFE structure for one sequence

# print(seq)
# print(structure)  # dot-bracket notation
# print("MFE (kcal/mol):", mfe)
N_CPU=multiprocessing.cpu_count()
print(N_CPU)
ray.init(num_cpus=N_CPU)
seqs, energies, labels=gtf.prep_structure_data()
print(seqs.shape, energies.shape, labels.shape)
ray.shutdown()
from mylib.bedgraph import bedgraph_for_all_samples, filter_sams
from mylib.metagene import bedgraph2pd, metagene_analysis
from mylib import GTF
if __name__ == "__main__":
    print("This is the mylib/bedgraph.py module.")
    SAM_FOLDER='/fs/ess/PAS2967/S21/Ribo-seq/bowtie'
    FNA='/fs/ess/PAS2967/S21/reference/GCF_000092025.1_ASM9202v1_genomic.fna'
    GTF_FILE='/fs/ess/PAS2967/S21/reference/genomic.gtf'
    GTF_ECOLI='/users/PAS0291/dengyw144/Fredrick008/ref/genomic.gtf'
    gtf= GTF(GTF_FILE, FNA)
    
    TEST_BEDGRAPH='/fs/ess/PAS2967/S21/Ribo-seq/coverage/Ribo-A1-plus.bedgraph'
    metagene_analysis(gtf, TEST_BEDGRAPH)
    gtf.gene_chrom_summary()
    # filter_sams(SAM_FOLDER, 18)
    # bedgraph_for_all_samples(SAM_FOLDER, ribo=True, gtf=gtf, cutoff=18)

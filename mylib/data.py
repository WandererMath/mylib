from mylib.gtf import Aligner, GTF
A=-7
B=25
C=7
# A irrelavant
TEMP=B-C
# s'= TEMP-s
# B-s' = index (staring from 0) of core ASD-SD left to right 
# n_NT extend to right = seq.find(']')- above -1


WEAK_RNA_RULE=[{'A', 'U'},{'C', 'G'}, {'G', "U"}]
#ASD_FJO='ACACCUCCUUUCU'
ASD_BA_SUB='UCACCUCCUUUCU'[::-1]
ASD_FJO='UCUUUCCUCCACA' # from published figure
asd_fjo=Aligner(ASD_FJO, WEAK_RNA_RULE)
asd_ba_sub=Aligner(ASD_BA_SUB, WEAK_RNA_RULE)


# Bacillus
FNA_FILE='../GCF_000009045.1/GCF_000009045.1_ASM904v1_genomic.fna'
GTF_FILE='../GCF_000009045.1/genomic.gtf'


#asd_ba_sub_strict=Aligner('UCACCUCCUUUCU', [{'A', 'U'},{'C', 'G'}])
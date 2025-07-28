from .gtf import GTF
import numpy as np
FNA='/fs/ess/PAS2967/S21/reference/GCF_000092025.1_ASM9202v1_genomic.fna'
GTF_FILE='/fs/ess/PAS2967/S21/reference/genomic.gtf'
gtf= GTF(GTF_FILE, FNA)
seqs, labels=gtf.prep_data_in_CDS_encoded()
labels=np.array(labels)
print(len(labels))
breakpoint()
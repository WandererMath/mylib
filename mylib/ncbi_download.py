import os
from subprocess import run




def download(acc):
    os.makedirs("download", exist_ok=True)
    os.makedirs("reference", exist_ok=True)
    cmd=f"datasets download genome accession {acc} --include genome,gtf,gff3 --filename download/{acc}.zip &&\
        unzip download/{acc}.zip -d download"
    run(cmd, shell=True, check=True)
    run(f"mv download/ncbi_dataset/data/{acc} reference/{acc}", shell=True, check=False)
    run(f"rm -r download", shell=True, check=False)

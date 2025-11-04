import os
from subprocess import run

def download(species, acc):
    os.makedirs("download", exist_ok=True)
    os.makedirs("reference", exist_ok=True)
    cmd=f"datasets download genome accession {acc} --include genome,gtf,gff3 --filename download/{acc}.zip &&\
        unzip download/{acc}.zip -d download"
    run(cmd, shell=True, check=True)
    run(f"mv download/ncbi_dataset/data/{acc} reference/{species}", shell=True, check=False)
    run(f"rm -r download", shell=True, check=False)

def batch_download(name_to_accession):
	for species, acc in name_to_accession.items():
		if os.path.isdir(f"reference/{species}"):
			continue
		print(f"Downloading {species} ({acc})")
		download(species, acc)

import os
import pandas as pd
import matplotlib.pyplot as plt

class Feature:
    '''
        Class to handle 1 featureCount file
    '''
    name: str
    counts: list    # int
    gene_ids: list # string
    def __init__(self, path, ribo_seq=False):
        df = pd.read_csv(path, sep="\t", comment="#")
        FIELD_COUNTS=df.columns.tolist()[-1]
        self.counts=df[FIELD_COUNTS].tolist()
        self.gene_ids=df['Geneid'].tolist()
        self.name=os.path.basename(path).split('.')[0]
        if ribo_seq:
            self.name='Ribo-'+self.name
    def scatter_plot(self, other, base_path):
        os.system(f'mkdir -p {base_path}')
        plt.scatter(self.counts, other.counts)
        plt.xscale('log')
        plt.yscale('log')
        plt.savefig(os.path.join(base_path, f"{self.name}-{other.name}.png"), dpi=300)
        plt.cla()
from qc_module import kmer_pca
import pandas as pd

pca_call = kmer_pca()
pca_df = pd.read_csv(snakemake.input[0], index_col=0)
lib_stats = pd.read_csv(snakemake.input[1], index_col=0)
meta_fields = [m for m in pca_df.columns if ~m.startswith('PC')]

pca_call.callModel(Y='PC1',X='Operator ID', data=pca_df)
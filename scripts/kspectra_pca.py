from qc_module import kmer_pca

if __name__ == '__main__':
    qc = kmer_pca()
    qc.generatePCA(kmer_spectra=qc.combine_kmerSpectra(snkmk_in=snakemake.input[0:-1]),
                   snkmk_out=snakemake.output, metadata=snakemake.input[-1])
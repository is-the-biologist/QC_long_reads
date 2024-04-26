import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.decomposition import PCA
import sys

class kmer_pca:
    """
    Module for performing necessary analysis for PCA generation from Jellyfish k-mer tables.
    """
    def __init__(self):
        pass

    def combine_kmerSpectra(self, snkmk_in):
        """
        Combine jellyfish kmer spectra tables
        :param snkmk_in:
        :return:
        """
        kmer_spectra = pd.concat([pd.read_csv(T, sep=' ', index_col=0, header=None) for T in snkmk_in], axis=1)
        kmer_spectra.columns = [name.split('/')[-1] for name in snkmk_in]
        kmer_spectra = kmer_spectra.fillna(0)

        return kmer_spectra

    def plotPCA(self, pc_comp, snkmk_out, meta=None):
            #generate the plots for the PCA of first two PCs
            with sns.axes_style('whitegrid'):
                if meta == None:
                    sns.scatterplot(data=pc_comp, x='PC1', y='PC2', edgecolor=None, alpha=0.5, color='black')
                else:
                    sns.scatterplot(data=pc_comp, x='PC1', y='PC2', edgecolor=None, alpha=0.5, hue=meta)
                plt.xticks(fontsize=15)
                plt.yticks(fontsize=15)
                plt.xlabel('PC1', fontsize=20)
                plt.ylabel('PC2', fontsize=20)
                plt.tight_layout()
                fname = snkmk_out[0].replace("None", str(meta))
                plt.savefig(fname, dpi=300)
                plt.close()

    def generatePCA(self, kmer_spectra, snkmk_out, metadata=None):
        """

        Performs PCA on the combined 21-mer content of each library from the Jellyfish output and returns a PCA of the first
        two components and the variance explained of each component (up to 100).

        :param kmer_spectra: input kmer spectra table from combined jellyfish input
        :param snkmk_out: snakemake.output
        :param metadata: Currently unused until I can get my hands on metadata files that I can use to plot against PCs
        :return:
        """

        #remove all k-mers where SD is 0 -- ie identical values:
        s = np.std(kmer_spectra.values, axis=1)
        kmer_spectra = kmer_spectra[s != 0]

        #Z-normalize the kmer values
        m = np.mean(kmer_spectra.values, axis=1)
        s = np.std(kmer_spectra.values, axis=1)

        z_kmer_spectra = ((kmer_spectra.values.transpose() - m) / s).transpose()

        #Run PCA to see if k-mer representations are very different
        pca_model = PCA(n_components=min([100, z_kmer_spectra.shape[1]])).fit(z_kmer_spectra.T)
        components = pca_model.transform(z_kmer_spectra.T)

        pc_comp = pd.DataFrame(data=components, columns=[f"PC{n+1}" for n in range(components.shape[0])],
                               index=[sname.split(".")[0] for sname in kmer_spectra.columns])

        if metadata != None: #if specified metadata exists
            meta_df = pd.read_csv(metadata, index_col=0)
            #list comprehension to include only metadata where fields are non-identical across samples but not unique for each sample
            meta_fields = [col for col in meta_df.columns if (len(set(meta_df[col].values)) > 1) and (len(set(meta_df[col].values)) < len(meta_df.index) )]
            meta_df = meta_df[meta_fields]
            meta_df.columns = [m.replace(" ", "_").replace("/", "_").replace("#", "Number") for m in meta_df.columns] #clean annoying characters
            pc_comp = pd.concat([pc_comp, meta_df], axis=1)
            #generate generic PCA
            self.plotPCA(pc_comp, snkmk_out)
            #generate colored PCA by meta fields:
            for meta in meta_df.columns:
                self.plotPCA(pc_comp, snkmk_out, meta)

        else: #ignore metadata
            #generate the plots for the PCA of first two PCs
            self.plotPCA(pc_comp, snkmk_out)

        #output variance explained by each component
        varExp = pd.DataFrame(data={'PC':[n+1 for n in range(components.shape[0])],
                                    'VarExplained':[f"{ev*100:0.2f}" for ev in pca_model.explained_variance_ratio_]})
        varExp.to_csv(snkmk_out[1], index=None)
        pc_comp.to_csv(snkmk_out[2])


class library_statistics:
    """
    Module for generating library statistic tables and plots directly from the BAM files.
    """
    def __init__(self):
        self.bin_num=100

    def read_in_bam_stream(self):
        """
        We are going to read in the stream from the BAM file directly and use an iterator function defined here which
        links back to the invoker script to generate a numpy array of the read lengths. This should have relatively small
        memory footprint even at >10M reads.
        """
        for line in sys.stdin:
            read_length = len(line.strip("\n"))
            yield read_length
        sys.stdin.close()

    def generateStatsTable(self, RL_array):
        """
        Compute simple summary statistics from the read length numpys to output as a table for ease of reference
        :param RL_array:
        :return:
        """
        T = np.sum(RL_array)
        mean_RL = T / len(RL_array)
        sd = np.std(RL_array)
        RL_percentiles = np.percentile(RL_array, [1 ,10, 25, 50, 75, 90, 99])
        stat_row = np.concatenate(([T, mean_RL, sd], RL_percentiles))
        return stat_row

    def histogram(self, RL_array, snk_out):
        """
        Generates histograms of read length for each library
        :param RL_array:
        :param snk_out:
        :return:
        """
        with sns.axes_style('whitegrid'):
            sns.histplot(x=RL_array, bins=self.bin_num, color='black')
            plt.xlabel('Read lengths', fontsize=20)
            plt.ylabel('Counts', fontsize=20)
            plt.xticks(fontsize=15, rotation=45)
            plt.yticks(fontsize=15)
            plt.tight_layout()
            plt.savefig(f"readlengths/{snk_out.split('/')[-1].split('npy')[0]}histogram.png", dpi=300)
            plt.close()

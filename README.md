# QC_long_reads

## Purpose:
This is a SnakeMake pipeline that will take in sample names in a config file and output a series of tables and plots that I think will be useful QC summary statistics for HiFi reads. 

## Dependencies:
All dependencies should be included in the `environment.yml` file, with the exception of `Jellyfish` (https://github.com/gmarcais/Jellyfish). To activate the environment you need conda or mamba installed and then simply do: `mamba env create --file=environment.yml`

## Inputs:
The Snakefile requires unaligned BAM files as input (in our case from HiFi, but conceivably from any sequencing platform). These can be specified within the `config.yaml` file. Simply modify the name and filepath of the entries within the `config.yaml` file with the BAM files in question before running. There is no-deduplication or adapter removal step in this workflow, but this can easily be modified by adding in those steps to before the calculation of k-mer spectra and read length statistics. 

## Outputs:
    plots/kmer_PCA.png
    plots/{sample}.histogram.png
    qc_tables/pca_var_explained.csv
    qc_tables/library_stats.csv
    

`kmer_PCA.png` provides a plot of the first two principal components of the k-mer spectra as computed by Jellyfish. The PCA is generated from the table of individual 21-mer counts normalized across samples. 21-mer was chosen as it is a common sized k-mer used for spectra statistics but could be changed to smaller or larger k-mer size. There is additionally the `pca_var_explained.csv` which provides the % variance explained by each PC in the PCA decomposition. Once metadata is available this output will be updated to mark samples by technical designations, which will allow us to interpret plots more easily. For example, if we find there is a clustering of samples in k-mer content by sequencing run this can indicate sequencing biases that will need to be accounted for downstream.

library_stats.csv is a summary table of read length and library size statistics of the sequencing libraries. The fields it contains are: 
    
    Total bp: size of library in bp
    Mean: average read length size
    SD: standard deviation of read lengths
    1-ile, 10-ile,25-ile,50-ile,75-ile,90-ile,99-ile: 1st to 99th percentile of read length sizes in the library.
    
A companion to this output are the {sample}.histogram.png plots that show the histogram of read length sizes of the libraries. Each library included in the config file will generate its own histogram.

## Intermediary files:
A number of intermediary files are generated to speed up re-running the pipeline when new samples are added. These are the `.jf` files and the `.npy` files. These contain the k-mer spectra and the read length distributions respectively.

# Running the pipeline:

To run it should be as simple as:

`mamba env create --file=environment.yml`

`mamba activate LongReadQC`

`Snakemake --cores 1`

# QC_long_reads

## Purpose:
This is a SnakeMake pipeline that will take in sample names in a config file and output a series of tables and plots that I think will be useful QC summary statistics for HiFi reads. 
Primarily, this pipeline focuses on analyses of k-mer spectra of the libraries. k-mers are representative of the sequence composition of libraries and there is evidence that different sequencer platforms have biases in k-mer content. This may prove useful to explore as we continue to assemble and analyze genomes. Particularly, regions of low complexity that can be biased by k-mer content more easily than single copy regions will be more strongly affected. Tandem repeats may also be biased by k-mer drop-out rates. Ultimately understanding at a broad scale the k-mer spectra can serve to QC libraries.

## Dependencies:
All dependencies should be included in the `environment.yml` file, with the exception of `Jellyfish` (https://github.com/gmarcais/Jellyfish). To activate the environment you need conda or mamba installed and then simply do: `mamba env create --file=environment.yml`

## Inputs:

### BAM files:
The Snakefile requires unaligned BAM files as input (in our case from HiFi, but conceivably from any sequencing platform). In the toy example I include HiFi reads from a human genome and from D melanogaster genome as well as Illumina reads from D melanogaster to illustrate stark k-mer differences. There is no-deduplication or adapter removal step in this workflow, but this can easily be modified by adding in those steps to before the calculation of k-mer spectra and read length statistics. 

### metadata.csv:
The metadata csv file contains the various meta-fields that may be important in assessing batch or technical effects in the data. In the toy example included in this repo I have added two additional columns one for "organism" and "chemistry". This won't be included in real metadata but is useful for the example workflow.

### config.yaml
.bam inputs can be specified in the `config.yaml` file. Simply modify the name and filepath of the entries within the `config.yaml` file with the BAM file paths in the "sample" section. Importantly, place the BAM files in "[]" such that each Shipping ID has all corresponding BAM files within the same ID. 
You must also specify a metadata.csv file to use for your samples, which MUST have identical Shipping ID column entries to the IDs within the 'sample' data structure. Finally, Jellyfish takes in an argument of has table size which can be modified by 'kmer_mem_usage'. `jellyfish -s 100M` sets 100 million entries in the hash table. The recommended size of the table is: _(G + k ∗ n)/0.8_, where _G_ is genome size _k_ is k-mer size and _n_ is the number of reads in the library. Increasing hash table size will use more memory and increase run time and may not be necessary as I exclude k-mers that are singletons in the library. 


## Outputs:
    plots/kmer_PCA.png
    plots/{sample}.histogram.png
    qc_tables/pca_var_explained.csv
    qc_tables/library_stats.csv
    qc_tables/meta_r2_effects.csv
    plots/meta_varExp.png


`kmer_PCA.png` provides a plot of the first two principal components of the k-mer spectra as computed by Jellyfish. The PCA is generated from the table of individual 21-mer counts normalized across samples. 21-mer was chosen as it is a common sized k-mer used for spectra statistics but could be changed to smaller or larger k-mer size. There is additionally the `pca_var_explained.csv` which provides the % variance explained by each PC in the PCA decomposition. Once metadata is available this output will be updated to mark samples by technical designations, which will allow us to interpret plots more easily. For example, if we find there is a clustering of samples in k-mer content by sequencing run this can indicate sequencing biases that will need to be accounted for downstream.

library_stats.csv is a summary table of read length and library size statistics of the sequencing libraries. The fields it contains are: 
    
    Total bp: size of library in bp
    Mean: average read length size
    SD: standard deviation of read lengths
    1-ile, 10-ile,25-ile,50-ile,75-ile,90-ile,99-ile: 1st to 99th percentile of read length sizes in the library.
    
A companion to this output are the {sample}.histogram.png plots that show the histogram of read length sizes of the libraries. Each library included in the config file will generate its own histogram.

The `meta_r2_effects.csv` and `meta_varExp.png` as the result of linear mixed-effect models run with PCs of Jellyfish spectra as dependent variables and metadata fields as dependent variables.
The idea behind this is quite simple, the PCs are a representation of the variation in k-mer spectra in the libraries and they should be independent of technical effects. Therefore, we use a linear
mixed-effect model where we model the PCs as the dependent variable and technical covariates (metadata) as random-effects and estimate the proportion of the variance in the PCs that can be explained by
the random effects. We additionally add in the fixed effect, total library size in bp (Z_total_bp), and estimate the effect of that as well. Additional, fixed-effects could be included as determined through
future analyses.

## Intermediary files:
A number of intermediary files are generated to speed up re-running the pipeline when new samples are added. These are the `.jf` files and the `.npy` files. These contain the k-mer spectra and the read length distributions respectively.

# Running the pipeline:

To run it should be as simple as:

`mamba env create --file=environment.yml`

`mamba activate LongReadQC.v1.1`

`snakemake --cores 1`

## DAG of rules:
![dag](https://github.com/is-the-biologist/QC_long_reads/assets/20618833/1112b665-9b02-494d-8938-9d30821cb2a6)

## Toy example:
Running this Snakefile as is will run a toy example wherein we generate k-mer spectra and library statistics for human HiFi, D melanogaster HiFi, and D melanogaster Illumina BAM files. This example illustrates how regression on the PCs of k-mer spectra can generate results that allow us to find covariates driving variation in our data.

![meta_varExp](https://github.com/is-the-biologist/QC_long_reads/assets/20618833/5894d147-341e-46fc-94e2-5b30fd5831e5)

From the result of the regression on the PCs we can see that most of the variation on PC1 is explained by "organism", followed by "chemistry" and then "Z_total_bp". You can see this pretty clearly in the PCA as well.

![kmer_organism_PCA](https://github.com/is-the-biologist/QC_long_reads/assets/20618833/47f94425-b173-4d1d-999a-a01061878209)

PC2 is largely explained by "chemistry" which may indicate the differences between Illumina and HiFi sequencing or adapter content in the Illumina library.

![kmer_chemistry_PCA](https://github.com/is-the-biologist/QC_long_reads/assets/20618833/404f0bb0-f4a5-47a9-bd5b-3882a2e9c238)

PC3 (which is not visualized in the scatterplot). Is largely driven by "Z_total_bp" which is just the Z-score normalized total bps of the library. This makes sense as greater sequencing depth will give greater abundance of common k-mers and greater variety of k-mer sequences.

Further down the PCs the interpretations are less clear as the % variance explained by these PCs is minimal.

It is important to keep in mind that this QC analysis is not meant to provide a perfect model of all the aspects of the HiFi data and completely diagnose problems. Rather, it is meant to serve as a starting point for exploratory analysis of batch effects of HiFi (or really any genomic data). 

## Parameters and modifications to consider:
There are a few modifications that we may need to do to efficiently run on large samples. If storage of read length numpy files and I/O operations are too costly then we should modify the scripts to not save read lengths to disk or store to memory and use file streaming to generate a simplified histogram. Secondly, 




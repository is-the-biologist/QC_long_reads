# QC_long_reads

##Purpose:
This is a SnakeMake pipeline that will take in sample names in a config file and output a series of tables and plots that I think will be useful QC summary statistics for HiFi reads. 

## Dependencies:
All dependencies should be included in the `environment.yml` file, with the exception of `Jellyfish` (https://github.com/gmarcais/Jellyfish).

## Inputs:
The Snakefile requires unaligned BAM files as input (in our case from HiFi, but conceivably from any sequencing platform). These can be specified within the `config.yaml` file. Simply modify the name and filepath of the entries within the `config.yaml` file with the BAM files in question before running. 

##Outputs:

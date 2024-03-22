configfile: "config.yaml" #modify config to add more samples in analysis
rule all:
    input:
        "plots/kmer_PCA.png",
        "qc_tables/pca_var_explained.csv",
        'qc_tables/library_stats.csv'

rule generate_stats:
    input:
        "bam_files/{sample}.bam"
    output:
        "readlengths_numpy/{sample}.npy"
    shell:
        "samtools view {input} -h | cut -f10 | python scripts/generateReadLengths.py {output}"

rule compute_library_stats:
    input:
        expand("readlengths_numpy/{sample}.npy", sample=config["samples"])
    output:
        "qc_tables/library_stats.csv",
    script:
        "scripts/libStats.py"

rule generate_jellyfish_spectra:
    input:
        "bam_files/{sample}.bam"
    output:
        "jellyfish/{sample}.jf"
    shell:
        "samtools fastq {input} | jellyfish count /dev/fd/0 -m21 -s 100M -L 2 -o {output}"

rule dump_jelly:
    input:
        "jellyfish/{sample}.jf"
    output:
        "jellyfish/{sample}.jf.txt"
    shell:
        "jellyfish dump -c {input} > {output}"

rule plotSpectraPCA:
    input:
        expand("jellyfish/{sample}.jf.txt", sample=config["samples"])
    output:
        "plots/kmer_PCA.png",
        "qc_tables/pca_var_explained.csv"
    script:
        "scripts/kspectra_pca.py"




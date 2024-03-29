configfile: "config.yaml"  #modify config to add more samples in analysis


def get_bam(wildcards):
    return config["samples"][wildcards.sample]

def concat_bam_stream(input, out):

    shell_commands =  "cat " + " ".join([f'<(samtools view {bam} -h | cut -f10)'for bam in input]) + f" | python scripts/generateReadLengths.py {out}"
    return shell_commands

def concat_jellyfish_stream(input, out):

    shell_commands =  "cat " + " ".join([f'<(samtools fastq {bam})'for bam in input]) + f" | jellyfish count /dev/fd/0 -m21 -s 100M -L 2 -C -o {out}"
    return shell_commands

rule all:
    input:
        "plots/kmer_None_PCA.png",
        "qc_tables/pca_var_explained.csv",
        "qc_tables/library_stats.csv",
        "qc_tables/meta_effects.csv"


rule generate_stats:
    input:
        get_bam,
    output:
        "readlengths_numpy/{sample}.npy",
    run:
        shell(concat_bam_stream(input=input, out=output))


rule compute_library_stats:
    input:
        expand("readlengths_numpy/{sample}.npy", sample=config["samples"]),
    output:
        "qc_tables/library_stats.csv",
    script:
        "scripts/libStats.py"


rule generate_jellyfish_spectra:
    input:
        get_bam,
    output:
        "jellyfish/{sample}.jf",
    run:
        shell(concat_jellyfish_stream(input, output))


rule dump_jelly:
    input:
        "jellyfish/{sample}.jf",
    output:
        "jellyfish/{sample}.jf.txt",
    shell:
        "jellyfish dump -c {input} > {output}"


rule plotSpectraPCA:
    input:
        expand("jellyfish/{sample}.jf.txt", sample=config["samples"]),
        config['metadata']
    output:
        "plots/kmer_None_PCA.png", #also will generate a variable # of metadata colored PCAs
        "qc_tables/pca_var_explained.csv",
        "qc_tables/pca_components.csv"
    script:
        "scripts/kspectra_pca.py"

rule spectraVarModels:
    input:
        "qc_tables/pca_components.csv",
        "qc_tables/library_stats.csv"
    output:
        "qc_tables/meta_effects.csv"
    script:
        'scripts/varExpModels.py'
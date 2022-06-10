configfile:"config.yaml"
cwd = os.getcwd()
accession=config["accession"]
rule all:
    input:
      "sleuth_object.so",
      "gene_table_results.txt"

rule fetch_FASTQ_from_SRA:
    output:
        temp("reads/{accession}/{accession}_1.fastq"),
        temp("reads/{accession}/{accession}_2.fastq")
    params:
        args = "--split-files --progress --details",
        accession = "{accession}"
    log:
        "reads/{accession}.log"
    conda:
        "env.yml"
    shell:
        'mkdir -p reads/{params.accession} && '
        'fasterq-dump {params.args} {params.accession} -O reads/{params.accession}'

rule downloadKallistoIndex:
    output:
      "homo_sapiens/Homo_sapiens.GRCh38.96.gtf",
      "homo_sapiens/Homo_sapiens.GRCh38.cdna.all.fa",
      "homo_sapiens/transcriptome.idx",
      "homo_sapiens/transcripts_to_genes.txt"
    conda:
      "env.yml"
    shell:
      'wget https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/homo_sapiens.tar.gz && \
      tar –xvzf homo_sapiens.tar.gz'

rule kallistoQuant:
    input:
        "homo_sapiens/Homo_sapiens.GRCh38.96.gtf",
        "homo_sapiens/Homo_sapiens.GRCh38.cdna.all.fa",
        "homo_sapiens/transcripts_to_genes.txt",
        index="homo_sapiens/transcriptome.idx",
        r1="reads/{accession}/{accession}_1.fastq",
        r2="reads/{accession}/{accession}_2.fastq",
    output:
        "quant/{accession}/abundance.h5",
        "quant/{accession}/abundance.tsv",
        "quant/{accession}/run_info.json"
    params:
        bootstrap="100",
        folder = "quant/{accession}"
    conda:
        "env.yml"
    shell:
        "mkdir -p {params.folder} && \
        kallisto quant -b {params.bootstrap} -i {input.index} -o {params.folder} {input.r1} {input.r2}"

rule mergeQuant:
    input:
        expand("quant/{accession}/abundance.h5",accession=config["accession"])
    output:
        out1="sleuth_object.so",
        out2="gene_table_results.txt"
    params:
      wd=cwd
    conda:
      "r.yml"
    shell:
      "Rscript scripts/sleuthR.R {params.wd}"

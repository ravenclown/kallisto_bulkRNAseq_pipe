configfile:"config.yaml"

rule all:
    input:
      "quant/abundance.h5",
      "quant/abundance.tsv",
      "quant/run_info.json"

rule fetch_FASTQ_from_SRA:
    output:
        temp("reads/{accession}_1.fastq"),
        temp("reads/{accession}_2.fastq")
    params:
        args = "--split-files --progress --details",
        accession = "{accession}"
    log:
        "reads/{accession}.log"
    conda:
        "env.yml"
    shell:
        'mkdir -p {params.accession}_reads && '
        'fasterq-dump {params.args} {params.accession} -O reads/'

rule downloadKallistoIndex:
    output:
      "homo_sapiens/Homo_sapiens.GRCh38.96.gtf",
      "homo_sapiens/Homo_sapiens.GRCh38.cdna.all.fa",
      "homo_sapiens/transcriptome.idx",
      "homo_sapiens/transcripts_to_genes.txt"
    conda:
      "env.yml"
    shell:
      "curl https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/homo_sapiens.tar.gz && \
      gunzip homo_sapiens.tar.gz"

rule kallistoQuant:
    input:
      "homo_sapiens/Homo_sapiens.GRCh38.96.gtf",
      "homo_sapiens/Homo_sapiens.GRCh38.cdna.all.fa",
      "homo_sapiens/transcriptome.idx",
      "homo_sapiens/transcripts_to_genes.txt",
      expand("reads/{accession}_1.fastq", accession=config["accession"]),
      expand("reads/{accession}_2.fastq", accession=config["accession"])
    output:
        "abundance.h5",
        "abundance.tsv",
        "run_info.json"
    conda:
      "env.yml"
    shell:
        "mkdir -p quant && kallisto quant -i homo_sapiens/transcriptome.idx -o quant/ reads/*.fastq"

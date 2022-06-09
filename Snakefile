configfile:"config.yaml"
cwd = os.getcwd()

rule all:
    input:
      "sleuth_object.so"
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
      "curl https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/homo_sapiens.tar.gz && \
      gunzip homo_sapiens.tar.gz"

rule kallistoQuant:
    input:
      "homo_sapiens/Homo_sapiens.GRCh38.96.gtf",
      "homo_sapiens/Homo_sapiens.GRCh38.cdna.all.fa",
      "homo_sapiens/transcriptome.idx",
      "homo_sapiens/transcripts_to_genes.txt",
      r1="reads/{accession}/{accession}_1.fastq",
      r2="reads/{accession}/{accession}_2.fastq"
    output:
        "abundance.h5",
        "abundance.tsv",
        "run_info.json"
    params:
      bootstrap="100"
    conda:
      "env.yml"
    shell:
        "mkdir -p quant/{wildcards.accession} && kallisto quant -b {params.bootstrap} -i homo_sapiens/transcriptome.idx -o quant/{wildcards.accession} {input.r1} {input.r2}

rule mergeQuant:
    input:
        expand("quant/{accession}/abundance.h5",accession=config["accession"])
    output:
        "sleuth_object.so"
        "gene_table_results.txt"
    params:
      wd=cwd
    conda:
      "r.yml"
    shell:
      "Rscript scripts/sleuthR.R {params.wd}"

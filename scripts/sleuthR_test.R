#!/usr/bin/env Rscript
#library calls
devtools::install_github('pachterlab/sleuth')
library("sleuth")
library("biomaRt")
library("dplyr")
#arguments
args = commandArgs(trailingOnly=TRUE)
base_dir <- "~/Downloads/kallisto_bulkRNAseq_pipe-main/"
sample_type <- "sample_conditions.txt"
setwd(base_dir)
#get sample id and directories
sample_id <- dir(file.path(base_dir,"quant"))
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir,"quant", id))
condition <- read.table(sample_type, header = TRUE, stringsAsFactors=FALSE) 
#get gene data from ensemble biomart
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
#Create a data frame for samples
s2c <- data.frame(sample=sample_id,condition=sample_type,path=kal_dirs,row.names = NULL)
#Run sleuth
so <- sleuth_prep(s2c,~condition,target_mapping = t2g)
so <- sleuth_fit(so)
so <- sleuth_wt(so,  'conditionscramble')
#Create table
gene_table <- sleuth_gene_table(so, test_type = "wt")
#Write table
write.table(gene_table, paste("gene_table_results.txt"), sep="\t")
#Save table
save(so, file=paste("sleuth_object.so"))

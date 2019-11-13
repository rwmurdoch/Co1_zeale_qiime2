## This is a simple script for combining outputs from qiime2 into a single OTU table which includes taxonomy and sequences
## The inputs are generated from qiime tools exports of sequences, feature table (via biom convert), and taxonomy

### This version of the script is based on a pipeline that does not do any feature taxonomic classification ### <-----

library(dplyr)
library(tidyverse)
library(Biostrings)

#read in the feature table (this is already sorted high to low by combined abundance)
feature.table <- read.csv("exports/otu_table.txt", skip = 1, comment.char = "", sep = "\t")
colnames(feature.table)[1] <- "OTU_ID"

## this section ignored in this version ##

#read in the taxonomy 
#taxonomy <- read.csv("exports/taxonomy.tsv",sep="\t", stringsAsFactors = F)
#colnames(taxonomy) [1] <- "OTU_ID"

#parse out the taxonomy using tidyverse functions as described here:
#https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121
#taxonomy <- taxonomy %>% as.tibble() %>% separate(Taxon, sep= ";", c("Kingdom","Phylum","Class","Order","Family","Genus","species"))

##

#read in the sequence fasta
q2.fasta <- Biostrings::readDNAStringSet("exports/dna-sequences.fasta")
OTU_ID <- names(q2.fasta)
sequence <- paste(q2.fasta)
q2.seq.tab <- data.frame(OTU_ID,sequence)

#join it all together

OTU_table <- dplyr::left_join(feature.table,q2.seq.tab,by="OTU_ID")

#write it
write.csv(OTU_table,"exports/combined_OTU_table.csv")

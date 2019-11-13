#This is a simple script for importing and size filtering a qiime2 ASV representative sequences list generated using 
#"qiime tools export" on the dada2 seqs file
#1. import with BioStrings
#2. make a simple table of "name" and "seq" and "length"
#3. filter the table by length
#4. remove length data from table
#5. write fasta using custom script
#This is then re-imported under the same name using "qiime tools import --type FeatureData[Sequence] (not by this script,
#rather by the master qiime2 script)

library(Biostrings)

#custom write fasta function
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

#import multifasta
seq.table <- Biostrings::readDNAStringSet("unfiltered.ASV.reps.fa/dna-sequences.fasta")

#extract relevant info
name <- names(seq.table)
seq <- paste(seq.table)
len <- width(seq.table)
seq.table <- data.frame(name, seq, len)

#filter by size
seq.table.filter <- seq.table[seq.table$len >= 147,]
seq.table.filter <- seq.table.filter[seq.table.filter$len <= 167,]

#write out the new fasta file
new.seq.table <- seq.table.filter[c(1,2)]
writeFasta(new.seq.table, "new.seqs.fa")

#prepare a qiime2 metadata file which contains id's to include in feature table
colnames(seq.table.filter)[1] <-"featureid"
write.csv(seq.table.filter[1], "include_names.csv", row.names = F)

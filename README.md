# Co1_zeale_qiime2

A simplified qiime2 pipeline for processing Zeale-primer Co1 amplicon data.

## Setup

this script was designed to work with qiime.2018.11

* place all files into the working directory
* create a folder called 'reads'
* create a folder called 'exports'
* this pipeline does NOT employ metadata

### optional extraction for nested fastq files

* create a folder called 'nested_reads'

An initial optional script will extract fastq files from a nested file structure and place in the 'reads' folder; place all subfolders into the 'reads_truncate' directory and execute the 'extract.nested.files.sh' script

## Pipeline overview

This pipeline executes a classic vsearch clustering approach

steps with (X) produce a log or visualization file

1. truncate all reads to 140bp
2. import truncated reads (X)
3. remove primers (X)
4. join pairs (X)
5. filter out low quality sequences (X)
6. detect and remove chimeras (check this step and adjust sensitivity if needed) (X)
7. cluster reads at 98.5%
8. apply a size filter, allowing only reads 147 - 167bp, via the Rscript 'size.filter.seqs.q2.R'.  Note that this script is run automatically from the shell script, but must be in the same folder (described in setup) (X)
9. remove OTUs which occur 10 or fewer times throughout the study (X)
10. export OTU sequence list and OTU tables <- these are the primary results, intended to be sent to taxonomic classification via BOLD.

For particular settings, please see the master qiime2 script itself ('qiime2willcox_classic.sh')

## Identification

The BOLD_results.Rmd r-markdown file is provided. The script itself has its own documentation. It will take sequence table and send each to the BOLD website for identification.  The result is an unwieldy xlsx file, which requires manual curation.

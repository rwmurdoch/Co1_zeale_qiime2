#!/bin/bash

#########################################
## This is an adaptation of the initial cprhd script
## with some reductions in the steps used, it uses a classic simliarity-based
## clustering pipeline via vsearch
## custom R scripts are used to apply length filter
## and produce a final OTU table
## no classification or secondary analyses are performed
## Nov 13, 2019
#########################################


## read truncation ##
# The nature of the overlap requires that the reads be
# truncated a bit (they extend through the entire amplicon)
# this is handled prior to importing into qiime2
# via a cutadapt loop through the reads folder

mkdir reads_truncate

for f in reads/*
	do
	name=${f##*/}
	cutadapt -l 140 -j 8 -o reads_truncate/"$name" "$f"
done

## conda activate qiime2-2018.11 ##

## data import ##

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path reads_truncate \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path cprhdreads.qza

qiime demux summarize \
--i-data cprhdreads.qza \
--o-visualization cprhd-demux.qzv \
--verbose

qiime cutadapt trim-paired \
--i-demultiplexed-sequences cprhdreads.qza \
--p-cores 8 \
--p-front-f AGATATTGGAACWTTATATTTTATTTTTGG \
--p-front-r WACTAATCAATTWCCAAATCCTCC \
--o-trimmed-sequences cprhd-data-cutadapt \
--verbose > cutadapt_log.txt

qiime demux summarize \
--i-data cprhd-data-cutadapt.qza \
--o-visualization cutadapt.qzv \
--verbose

#1#
#join pairs with some strictness and QC filtering#
#this is used only for the VSearch  and deblur pipelines

qiime vsearch join-pairs \
--i-demultiplexed-seqs cprhd-data-cutadapt.qza \
--p-minovlen 20 \
--o-joined-sequences cprhd-data-joined \
--verbose

qiime demux summarize \
--i-data cprhd-data-joined.qza \
--o-visualization joined.qzv \
--verbose

#2#
#this is a strict quality filtering#

qiime quality-filter q-score-joined \
--i-demux cprhd-data-joined.qza \
--p-min-quality 25 \
--o-filtered-sequences cprhd-seqs-strict-qtrim \
--o-filter-stats cprhd-seqs-strict-qtrim-stats \
--verbose

qiime demux summarize \
--i-data cprhd-seqs-strict-qtrim.qza \
--o-visualization qtrim.qzv \
--verbose

#3#
#this is required, creates a feature table#
#very time consuming step#

qiime vsearch dereplicate-sequences \
  --i-sequences cprhd-seqs-strict-qtrim.qza \
  --o-dereplicated-table cprhd-denovo-derep-table \
  --o-dereplicated-sequences cprhd-denovo-derep-seqs

#4#
# chimera removal #

qiime vsearch uchime-denovo \
  --i-sequences cprhd-denovo-derep-seqs.qza \
  --i-table cprhd-denovo-derep-table.qza \
  --o-chimeras cprhd-denovo-chimera-seqs \
  --o-nonchimeras cprhd-denovo-nonchimera-seqs \
  --o-stats cprhd-denovo-chimera-stats

# 6 #
# after separating chimeras, i have to filter the input data, removing chimeras#

qiime feature-table filter-features \
  --i-table cprhd-denovo-derep-table.qza \
  --m-metadata-file cprhd-denovo-chimera-seqs.qza \
  --p-exclude-ids \
  --o-filtered-table cprhd-denovo-filtered-table.qza

qiime feature-table filter-seqs \
  --i-data cprhd-denovo-derep-seqs.qza \
  --m-metadata-file cprhd-denovo-chimera-seqs.qza \
  --p-exclude-ids \
  --o-filtered-data cprhd-denovo-filtered-seqs.qza

qiime feature-table summarize \
  --i-table cprhd-denovo-filtered-table.qza \
  --o-visualization cprhd-denovo-filtered-table.qzv

#3#
#this is the clustering but does not actually create a taxonomy table#
#essentially, clusters are created primarily using the ref database as seed

qiime vsearch cluster-features-de-novo \
--i-sequences cprhd-denovo-filtered-seqs.qza \
--i-table cprhd-denovo-filtered-table.qza \
--p-perc-identity 0.985 \
--p-threads 8 \
--o-clustered-table cprhd-denovo-OTU-table \
--o-clustered-sequences cprhd-denovo-OTU-seqs

###########################################
# Size filtering of reads/ASVs #
#########################################

# this is a custom step which is not built into qiime2
# basic strategy is to export the vsearch_seqs, apply a size filter in R
# and re-import under the same name

qiime tools export \
  --input-path cprhd-denovo-OTU-seqs.qza \
  --output-path unfiltered.ASV.reps.fa

#this R script handles size filtering, removing sequences
#from the "unfiltered.ASV.reps.fa" file and creating a "new.seqs.fa"
#which is then moved in and replaces *-vsearch.seqs
#upper and lower size limits must be altered in the script directly

Rscript size.filter.seqs.q2.R

#re-import the size filtered seq set
qiime tools import \
  --type FeatureData[Sequence] \
  --input-path new.seqs.fa \
  --output-path cprhd-denovo_seqs

#filter the feature table by retained feature list
qiime feature-table filter-features \
  --i-table cprhd-denovo-filtered-table.qza \
  --m-metadata-file include_names.csv \
  --p-no-exclude-ids \
  --o-filtered-table cprhd-denovo-filtered-table

qiime feature-table summarize \
--i-table cprhd-denovo-filtered-table.qza \
--o-visualization cprhd-vsearch-seq-stats \
--verbose

#remove low frequency features
qiime feature-table filter-features \
  --i-table cprhd-denovo-filtered-table.qza \
  --p-min-frequency 11 \
  --o-filtered-table cprhd-denovo-filtered-table-filtered
  
#remove low frequency sequences
qiime feature-table filter-seqs \
	--i-data cprhd-denovo_seqs.qza \
	--i-table cprhd-denovo-filtered-table-filtered.qza \
	--o cprhd-denovo_seqs.qza

qiime feature-table summarize \
	--i-table cprhd-denovo-filtered-table-filtered.qza \
	--o-visualization cprhd-vsearch-seq-stats \
	--verbose

#######################################
# create OTU table and export seq and feature otu_table
# unique to this denovo pipeline
######################################

qiime tools export --input-path cprhd-denovo-filtered-table-filtered.qza --output-path exports
qiime tools export --input-path cprhd-denovo_seqs.qza --output-path exports

####################################
# sample visualization
###############################

## There is no metadata for this project so most of the remaining
## steps are not used

#qiime feature-table summarize \
#  --i-table cprhd-denovo-filtered-table-filtered.qza \
#  --o-visualization cprhd-denovo-filtered-table.qzv \
#  --m-sample-metadata-file cprhd_metadata.csv

qiime feature-table tabulate-seqs \
  --i-data cprhd-denovo-OTU-seqs.qza \
  --o-visualization cprhd-denovo_seqs.qzv

#######################################
# preparation for alpha div metrics
#######################################

#qiime alignment mafft \
#  --i-sequences cprhd-denovo-OTU-seqs.qza \
#  --o-alignment cprhd-denovo_seqs-aligned.qza \
#  --p-n-threads 8

#qiime alignment mask \
#--i-alignment cprhd-denovo-OTU-seqs-aligned.qza \
#--o-masked-alignment cprhd-denovo_seqs-aligned-masked.qza

#qiime phylogeny fasttree \
# --i-alignment cprhd-denovo-OTU-seqs-aligned-masked.qza \
# --o-tree cprhd-unrooted-tree-vsearch.qza

#qiime phylogeny midpoint-root \
#--i-tree cprhd-unrooted-tree-vsearch.qza \
#--o-rooted-tree cprhd-rooted-tree-vsearch.qza

################################################
# generate core diversity metrics #
# pay close attention to the sampling depth command #
###################################################

#qiime diversity core-metrics-phylogenetic \
#  --i-phylogeny cprhd-rooted-tree-vsearch.qza \
#  --i-table cprhd-denovo-filtered-table-filtered.qza \
#  --p-sampling-depth 30000 \
#  --m-metadata-file cprhd_metadata.csv \
#  --output-dir cprhd-core-metrics-results-vsearch

## visualization of diversity metrics ##

#qiime diversity alpha-rarefaction \
#  --i-table cprhd-denovo-filtered-table-filtered.qza \
#  --i-phylogeny cprhd-rooted-tree-vsearch.qza \
#  --p-min-depth 1 \
#  --p-max-depth 10000 \
#  --m-metadata-file cprhd_metadata.csv \
#  --o-visualization cprhd-denovo-alpha-rarefaction-vsearch.qzv

###############################################
# exporting data
# you can easily export seqs and taxonomy directly out of
# their corresponding .qza files using the qiime tools export feature.
# to get the feature table, you export the table in biom format and then convert
# into a simple text file
#################################

#qiime tools export --input-path cprhd-denovo-filtered-table-filtered.qza --output-path exports
#qiime tools export --input-path cprhd-denovo-OTU-seqs.qza --output-path exports
#qiime tools export --input-path cprhd-vsearch-taxonomy.qza --output-path exports

biom convert -i exports/feature-table.biom \
-o exports/otu_table.txt --to-tsv

## combine everything into a single table
## this uses the R script create.OTU.table.R
## script is available here: https://github.com/rwmurdoch/project.scripts/blob/master/create.OTU.table.R
## download/copy the script into your project directory; it will read and write in the "export" directory

#sudo Rscript create.OTU.table.R

## combine various outputs into a results package

#tar -cf results.tar.gz \
#exports \
#cprhd-vsearch-taxa-bar-plots.qzv \
#cprhd_metadata.csv \
#cprhd-vsearch-seq-stats.qzv \
#cprhd-demux.qzv

#cprhd-core-metrics-results-vsearch \
#cprhd-denovo-alpha-rarefaction-vsearch.qzv \

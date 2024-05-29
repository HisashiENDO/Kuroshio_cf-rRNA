#!/bin/sh
source /etc/profile.d/modules.sh

## This script aims to make ASV profile for KS-22-15 curated sequencing data

# 20240224 Update the PR2 fasta and tax files from ver.4.14.0 to ver.5.0.0, which greatly improve some taxonomic annotation.

## Prepare high-quality and length unified fastq file by the pipeline in the following directory

## The general procedures as follows
#1. Importing raw sequence data
#2. Trim primer sequences 
#3. DADA2 denoising, defining ASVs, and making table
#4. Remove sigleton
#5. Rarefy ASVs
 # ->Then make visual sumamries for the qiime2 view
#6. Make trained classifer for the taxonomic analysis
#7. Classify ASVs into taxonomy
#8. Make bar plot of summarizing taxonomy levels
#9. Ordination based on Unifrac distance etc.
#10. Output ASV and taxonomy tables as tsv format



# Activate QIIME2
# source activate qiime2-2020.2

# Load module
module load qiime2/2023.9
module load Python/3.11.7

# Set working directory

#this file can be gerenated with a script (0 MANIFEST), but manuall creation is recommended.
# tutorial: https://docs.qiime2.org/2018.11/tutorials/importing/
MANIFEST=manifest.csv

#This file is needed for a barplot and other plots and has to be cre ated manually
# tutorial: https://docs.qiime2.org/2018.11/tutorials/metadata/
METADATA_FILE=./sample_metadata.tsv


############ Import sequences as qiime2 "qza" format ############

## https://docs.qiime2.org/2020.2/tutorials/importing/#sequences-without-quality-information-i-e-fasta
echo "Importing sequence data using $MANIFEST"
qiime tools import \
 --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path $MANIFEST \
 --output-path seqs.qza \
 --input-format PairedEndFastqManifestPhred33

# Trim primer sequences from the raw reads
echo "Trim primer sequences from the raw reads"
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences seqs.qza \
  --o-trimmed-sequences seqs.trimmed.qza \
  --p-front-f CYGCGGTAATTCCAGCTC \
  --p-front-r AYGGTATCTRATCRTCTTYG \
  --p-error-rate 0.1 \
  --p-cores 12 



############ Denoizng by DADA2 ##############

# https://docs.qiime2.org/2023.2/plugins/available/dada2/denoise-paired/
echo "DADA2 denoising, defining ASVs, and making table"
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs seqs.trimmed.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 280\
  --p-trunc-len-r 250\
  --p-n-threads 12 \
  --o-table seqs.table.qza \
  --o-representative-sequences seqs.repsequence.qza \
  --o-denoising-stats seqs.dada2_denoising-stats.qza

# Convert to qzv format
echo "Converting DADA2 results to qzv format"
## Make visualization qza format
qiime metadata tabulate \
  --m-input-file seqs.dada2_denoising-stats.qza \
  --o-visualization seqs.dada2_denoising-stats.qzv
## Make ASV table
qiime metadata tabulate \
  --m-input-file seqs.table.qza \
  --o-visualization seqs.table.qzv


## Remove rare ASVs
echo "Filtering rare ASVs"
# Here I remove ASVs appear less than 10 times across samples, as well as those appear only 1 sample.
# These ASVs are not useful for the downstream analysis.
qiime feature-table filter-features \
  --i-table seqs.table.qza \
  --p-min-frequency 10 \
  --p-min-samples 2 \
  --o-filtered-table seqs.table.remrare.qza
## Make ASV table
qiime metadata tabulate \
  --m-input-file seqs.table.remrare.qza \
  --o-visualization seqs.table.remrare.qzv


##### Summarize feature-table #####
echo "Make visual sumamries of the data"
qiime feature-table summarize \
  --i-table seqs.table.remrare.qza \
  --o-visualization seqs.table.rarefy.sum.qzv \
  --m-sample-metadata-file $METADATA_FILE

# Make representative ASV sequence table with statistics
qiime feature-table tabulate-seqs \
  --i-data seqs.repsequence.qza \
  --o-visualization seqs.repsequence.qzv


###
# I do not conduct rafarication of reads across samples
# Unifying the read counts will be done by R vegan function "rarefy" 
###



############## Taxonomic assignment ################

echo "Train reference database used for the taxonomic annotation"


### Generete trained reference sequences for Comeau's V4 primer ###

# Full length silva dataset: https://docs.qiime2.org/2023.2/data-resources/
# Load taxonomic data
#wget "https://data.qiime2.org/2020.8/common/silva-138-99-seqs.qza"
#wget "https://data.qiime2.org/2020.8/common/silva-138-99-tax.qza"
##! If you use aptmp, you can link to "/aptmp/endo/Data/silva138/qiime2/"
# not included in qiime2 environent of own PC

# Define taxonomy data
#Silva138_seqs=/aptmp/endo/Data/silva138/qiime2/silva-138-99-seqs.qza
#Silva138_tax=/aptmp/endo/Data/silva138/qiime2/silva-138-99-tax.qza

# This time we use PR2 database to annotate sequences instead of using silva database.

# Prepare PR2 database for qimme2: https://rpubs.com/aschrecengost/794628
#get fasta and taxonomy info
#wget https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_mothur.fasta.gz
#wget https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_mothur.tax.gz
# Then gunzip these files


# Skip following part as it was already done and the curated clasifier was copied in Data directory


PR2_seqs=/aptmp/endo/Data/PR2/pr2_version_5.0.0_SSU_mothur.fasta
PR2_tax=/aptmp/endo/Data/PR2/pr2_version_5.0.0_SSU_mothur.tax

# Import reference sequences as qiime artifacts
qiime tools import \
  --input-path $PR2_seqs \
  --output-path pr2_5.0.0_fasta.qza \
  --type 'FeatureData[Sequence]'

# Import reference taxonomy as qiime artifacts
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path $PR2_tax \
  --output-path pr2_version_5.0.0_tax.qza

# Extract target seqience region to the training
# https://docs.qiime2.org/2023.2/tutorials/feature-classifier/
# https://docs.qiime2.org/2023.2/plugins/available/feature-classifier/extract-reads/
qiime feature-classifier extract-reads \
  --i-sequences pr2_5.0.0_fasta.qza \
  --p-f-primer CYGCGGTAATTCCAGCTC \
  --p-r-primer AYGGTATCTRATCRTCTTYG \
  --p-trunc-len 0 \
  --p-min-length 100 \
  --p-max-length 0 \
  --o-reads pr2_5.0.0_fasta.E572F_E1009R.curated.qza \
  --p-n-jobs 12

# https://rpubs.com/aschrecengost/794628
# この処理は時に必要な配列も削除してしまい、場合によってannotation結果を悪化させることもあるので注意


# Train references targeting the target sequences that was sequenced
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads pr2_5.0.0_fasta.E572F_E1009R.curated.qza \
  --i-reference-taxonomy pr2_version_5.0.0_tax.qza \
  --o-classifier pr2_5.0.0_fasta.E572F_E1009R.curated.tax.classifier.qza


PR2_classifier=/aptmp/endo/Data/PR2/pr2_5.0.0_fasta.E572F_E1009R.curated.tax.classifier.qza

### Assign taxonomy for each ASVs ###
echo "Classify taxonomy by q2-feature-classifier plugin"

# Classify taxonomy
qiime feature-classifier classify-sklearn \
  --i-classifier $PR2_classifier \
  --i-reads seqs.repsequence.qza \
  --o-classification seqs.repsequence.taxonomy.qza

# Output visualization
qiime metadata tabulate \
  --m-input-file seqs.repsequence.taxonomy.qza \
  --o-visualization seqs.repsequence.taxonomy.qzv



############## Filtered out uninterested taxa ################

echo "Here I remove metazoan, unassigned, and bacterial ASVs from the table"

qiime taxa filter-table \
  --i-table seqs.table.remrare.qza \
  --i-taxonomy seqs.repsequence.taxonomy.qza \
  --p-mode contains \
  --p-exclude 'Eukaryota;Obazoa;Opisthokonta;Metazoa, Unassigned, Bacteria;' \
  --o-filtered-table seqs.table.remrare.protist.qza



############## Generate some preliminary plots ################

# Bar plot
# Metazoa included
qiime taxa barplot \
  --i-table seqs.table.remrare.qza \
  --i-taxonomy seqs.repsequence.taxonomy.qza \
  --m-metadata-file $METADATA_FILE \
  --o-visualization seqs.taxa-bar-plots.all.qzv

# Metazoa excluded
qiime taxa barplot \
  --i-table seqs.table.remrare.protist.qza \
  --i-taxonomy seqs.repsequence.taxonomy.qza \
  --m-metadata-file $METADATA_FILE \
  --o-visualization seqs.taxa-bar-plots.protist.qzv



############## Generate class level (PR2 level5) profiles ################

## Make level 5 (class level)
# Metazoa included
qiime taxa collapse \
  --i-table seqs.table.remrare.qza \
  --i-taxonomy seqs.repsequence.taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table seqs.table.remrare.table.lv5.qza
# Output visualization
qiime metadata tabulate \
  --m-input-file seqs.table.remrare.table.lv5.qza \
  --o-visualization seqs.table.remrare.table.lv5.qzv

# Metazoa excluded
qiime taxa collapse \
  --i-table seqs.table.remrare.protist.qza \
  --i-taxonomy seqs.repsequence.taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table seqs.table.remrare.protist.table.lv5.qza
# Output visualization
qiime metadata tabulate \
  --m-input-file seqs.table.remrare.protist.table.lv5.qza \
  --o-visualization seqs.table.remrare.protist.table.lv5.qzv


## Make level 6 (order level)
# Metazoa included
qiime taxa collapse \
  --i-table seqs.table.remrare.qza \
  --i-taxonomy seqs.repsequence.taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table seqs.table.remrare.table.lv6.qza
# Output visualization
qiime metadata tabulate \
  --m-input-file seqs.table.remrare.table.lv6.qza \
  --o-visualization seqs.table.remrare.table.lv6.qzv

# Metazoa excluded
qiime taxa collapse \
  --i-table seqs.table.remrare.protist.qza \
  --i-taxonomy seqs.repsequence.taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table seqs.table.remrare.protist.table.lv6.qza
# Output visualization
qiime metadata tabulate \
  --m-input-file seqs.table.remrare.protist.table.lv6.qza \
  --o-visualization seqs.table.remrare.protist.table.lv6.qzv


############# output results ############

echo "Generate output files for abundance table, taxonomy, and repsequences"

mkdir general_outputs

## Make ASV table for curated table
qiime metadata tabulate \
  --m-input-file seqs.table.remrare.protist.qza \
  --o-visualization seqs.table.remrare.protist.tabulate.qzv
  # Compare with "seqs.table.remrare.qzv" to calcurate protist/(protist+metazoa)

# Output data table as a tsv format
qiime tools extract \
  --input-path seqs.table.remrare.protist.tabulate.qzv \
  --output-path general_outputs

# Output taxonomy table as a tsv format
qiime tools extract \
  --input-path seqs.repsequence.taxonomy.qzv \
  --output-path general_outputs

# Output representative sequences as a tsv format
qiime tools extract \
  --input-path seqs.repsequence.qzv \
  --output-path general_outputs



# Execute with
# qsub -q cdb -l select=1:ncpus=12:mem=64gb run_qiime2_add.sh
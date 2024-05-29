#!/bin/bash

# This scipt aims to do fastp QC for all the file listed in the text file

# !If you made list.txt in excel, you have to convert "Line feed code" from CR+LF to \n
##cf. nkf -Lu --overwrite seqhead_list.txt

cd ~/Library/CloudStorage/Dropbox/2.Analysis_DB/Endo_cfRNA/00_QC

# Make diretoies that the results will be stored
mkdir ./trimmed_seqs
mkdir ./trimmed_QC


# Run fastp
for file in `cat seqhead_list.txt`;
do
  fastp \
   --detect_adapter_for_pe \
   --qualified_quality_phred 30 \
   -i ../18S_raw_data/raw_data/HN00188688/${file}_1.fastq.gz \
   -I ../18S_raw_data/raw_data/HN00188688/${file}_2.fastq.gz \
   -o ./trimmed_seqs/${file}_1.trim.fastq.gz \
   -O ./trimmed_seqs/${file}_2.trim.fastq.gz \
   -h ./trimmed_QC/${file}_report.html \
   -j ./trimmed_QC/${file}_report.json \
   -w 6;
done

<< COMMENTOUT
COMMENTOUT

# Execute with sh _fastp_all.sh 
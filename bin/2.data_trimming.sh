#!/bin/bash

# This script is for data cleaning adapters using trimmomatic.
# Run this script from directory ~/bin/  and  the sequences are in ~/data/DE
# run trought conda


# make out directory for DE clean data

mkdir -p /projects/irosas/data/tree_3_clean


# Clean the DE sequences with trimmomatic (adaptor removal, trimming of low quality bases and reads)

for i in `ls /projects/irosas/data/raw_data/tree_3/ | grep ".fq.gz" | sed "s/_1.fastq.gz//"| sed "s/_2.fastq.gz//" | uniq` ; do
echo ${i}
                           trimmomatic PE -threads 8 -phred33 \
                          -trimlog ./projects/irosas/data/tree_3_clean/${i}_trimlog.txt \
                          /projects/irosas/data/raw_data/tree_3/${i}_1.fastq.gz /projects/irosas/data/raw_data/tree_3/${i}_2.fastq.gz \
                          /projects/irosas/data/tree_3_clean/${i}_1P.fastq.gz  /projects/irosas/data/tree_3_clean/${i}_1U.fastq.gz \
                          /projects/irosas/data/tree_3_clean/${i}_2P.fastq.gz  /projects/irosas/data/tree_3_clean/${i}_2U.fastq.gz \
                          ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
                          LEADING:15 TRAILING:15 MINLEN:75 SLIDINGWINDOW:4:25

done

PCDT3-10b_1.fastq.gz
PCDT3-10b_2.fastq.gz

/home/irosas


/projects/irosas/data/raw_data

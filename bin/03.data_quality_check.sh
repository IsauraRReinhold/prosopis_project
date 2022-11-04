#!/bin/bash
# In order to evaluate que quality of the trimmed data run this script
# script to evaluate the quality of trimm_samples using fastqQC and multiQC
# Run this script from directory ~/bin/  clean data is saved in ~/data/
# Prerequisites: fastQC and multiQC

for i in tree_3_clean tree_4_clean tree_5_clean; do

#make out directories

mkdir -p /projects/irosas/out/${i}_quality

#evaluate raw quality samples with fasqc and save the out in dir quality

for k in /projects/irosas/data/${i}_clean/*.fastq.gz; do

fastqc $k -o /projects/irosas/out/${i}_quality/

  done
done

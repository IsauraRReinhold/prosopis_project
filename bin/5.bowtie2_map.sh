#!/bin/bash
#script to map reads vs bowtie2_index
#run this script from bin/
#Prerequisites: bowtie2



#run bowtie2 alignmet using the index and the clean data

for i in data_2021; do
#create out dir
mkdir -p /projects/irosas/out/${i}_mapping

for k in `ls /projects/irosas/data/${i}_clean |  grep -Ev "(_1U.fastq.gz|_2U.fastq.gz|.txt)" | sed "s/_1P.fastq.gz//"| sed "s/_2P.fastq.gz//" | uniq ` ; do

echo "Processing sample ${k}"

bowtie2 --local --no-unal -p 16 -x /projects/irosas/out/prosopis_index/prosopis_index \
        -1 /projects/irosas/data/${i}_clean/${k}_1P.fastq.gz  -2 /projects/irosas/data/${i}_clean/${k}_2P.fastq.gz  \
        -S /projects/irosas/out/${i}_mapping/${k}_P.sam 2> ${k}_bowtie.log
done
 done

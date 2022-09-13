#!/bin/bash
#script to create count matrix using sorted.bam files
#run this script from bin/
#Prerequisites: subread

#run bowtie2 alignmet using the index and the clean data

#create out dir

mkdir -p /media/cris/Isaura/prosopis_project/out/count_matrix/

#use featureCounts to create the count matrix

featureCounts -T 4 -t "exon" -g "transcript_id" -p -s 2 -O -a /media/cris/Isaura/prosopis_project/data/annotation/PC_final_gene_all_function_stringtie.gtf -o /media/cris/Isaura/prosopis_project/out/count_matrix/prosopis_count_matrix_P2.txt /media/cris/Isaura/prosopis_project/out/sorted/*.sorted.bam
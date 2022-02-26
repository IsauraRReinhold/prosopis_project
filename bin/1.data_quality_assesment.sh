#!bin/bash
#In order to evaluate que quality of my raw data I run this script
#script to evaluate the quality of raw_samples using fastqc
# Run this script from directory ~/bin/  and  the raw data is in ~/data/

## folders with raw data are saved in ../data/raw_data

for i in tree_3 tree_4 tree_5;
do

#make out folders

  mkdir -p ../data/raw_data/${i}_quality

  for k in ../data/raw_data/${i}/*.gz; do
  echo ${k}

  fastqc ${k} -o ../data/raw_data/${i}_quality/ #run fastqc in every sample saved in prosipis samples and save the ouput

  done
done

#made multiqc analysis in tree_3_quality tree_4_quality tree_5_quality directories

for i in tree_3 tree_4 tree_5; do

cd ../data/raw_data/${i}_quality/

multiqc .  #run multiqc analysys inside folfers

done

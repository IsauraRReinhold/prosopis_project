

trimmomatic PE -threads 4 -phred33 \
-trimlog /media/cris/Isaura/prosopis_project/data/data_2021_clean/S21DT3_10_2_trimlog.txt \
/media/cris/Isaura/prosopis_project/data/raw_data/data_2021/movidos/S21DT3_10_2_1.fq.gz /media/cris/Isaura/prosopis_project/data/raw_data/data_2021/movidos/S21DT3_10_2_1.fq.gz \
/media/cris/Isaura/prosopis_project/data/data_2021_clean/S21DT3_10_2_1P.fastq.gz /media/cris/Isaura/prosopis_project/data/data_2021_clean/S21DT3_10_2_1U.fastq.gz \
/media/cris/Isaura/prosopis_project/data/data_2021_clean/S21DT3_10_2_2P.fastq.gz  /media/cris/Isaura/prosopis_project/data/data_2021_clean/S21DT3_10_2_2U.fastq.gz \
ILLUMINACLIP:/home/cris/Documentos/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
LEADING:15 TRAILING:15 MINLEN:75 SLIDINGWINDOW:4:25

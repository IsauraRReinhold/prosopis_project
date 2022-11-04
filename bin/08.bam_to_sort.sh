
#!/bin/bash
#script to map reads vs bowtie2_index
#run this script from bin/
#Prerequisites: bowtie2


#run bowtie2 alignmet using the index and the clean data

for i in tree_3 tree_4 tree_5; do
#create out dir

for k in `ls /projects/irosas/out/${i}_mapping |  grep "_P.bam" | sed "s/_P.bam//" | uniq` ; do

echo "Processing sample ${k}"

samtools sort /projects/irosas/out/${i}_mapping/${k}_P.bam -o /projects/irosas/out/${i}_mapping/${k}_P.sorted.bam

done
done

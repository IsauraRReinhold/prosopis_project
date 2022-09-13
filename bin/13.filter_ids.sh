


for i in tree_3 tree_4 tree_5; do
  for k in `ls /media/cris/Isaura/prosopis_project/out/${i}_mapping/ | grep ".fasta" | sed "s/.fasta//"| uniq` ; do
  echo ${k}
perl -pe '$i=$1if/^>(\S+)/;map$i{$_}++,split;$i{$i}or$_=""' /home/cris/Documentos/Prosopis_project/out/specific_ids.txt /media/cris/Isaura/prosopis_project/out/${i}_mapping/${k}.fasta > /media/cris/Isaura/prosopis_project/out/${i}_mapping/${k}_filter.fasta

done
done


#!/bin/bash

while read url; do
    wget $url
done < /home/cris/Documentos/Prosopis_project/metadata/Transcriptome_2021_data_download_link.txt

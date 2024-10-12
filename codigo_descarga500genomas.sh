#!/bin/bash                                                                    
#:Title: Genomic Vigilance Staphylococcus pseudintermedius                                                                           
#:Date: 07-08-2024                                                              
#:Author: "Gerald Moreno" <gmoreno993@gmail.com>                                
#:Version: 1.0                                                                  
#:Description : Genomic Vigilance Staphylococcus pseudintermedius                             
#:Options: None      

#Flujograma:
#1.Descargar genomas de NCBI (assemblies)
#Descargar genomas de Staphylococcus aureus Patric
cd /media/gerald/DC604FBB604F9B62/sa_congreso/patric

#Descargar herramienta datasets de NCBI de https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/
#Colocar el ejectuable en la carpeta donde esta tu archivo .txt de las accesiones
##extraer la fila de codigos de accesion manualmente y 

#Descargar cada muestra
# Archivo de entrada con los códigos de acceso (lista de las cepas con las que deseas trabajar
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets #datasets
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat #dataformat



input_file="patric_assembly_list.txt"
# Leer cada línea del archivo y descargar los archivos
while IFS= read -r accession; do
    echo "Descargando $accession..."
    ./datasets download genome accession $accession --filename "${accession}.zip"
done < "$input_file"

echo "Descargas completas."

#No olvidar descargar herramienta datasets de NCBI de https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/


cd /media/gerald/DC604FBB604F9B62/sa_congreso/patric/genomes_sa_patric_southamerica
for archivo in *.zip; do
  nombre_carpeta="${archivo%.zip}"
  mkdir -p "$nombre_carpeta"
  unzip "$archivo" -d "$nombre_carpeta"
done

#Sacar todos los .fna que esten en carpetas internas
find /media/gerald/DC604FBB604F9B62/sa_congreso/patric -type f -name "*.fna" -exec mv {} /media/gerald/DC604FBB604F9B62/sa_congreso/patric \;
mkdir genomes_sa_patric_southamerica
mv *.fna genomes_sa_patric_southamerica #guardar todos los .fna en una carpeta

#Se usará: Staphopia: https://github.com/staphopia/staphopia-sccmec/issues/7


conda create -n staphopia-sccmec -c conda-forge -c bioconda staphopia-sccmec

conda activate staphopia-sccmec 

#una cepa
staphopia-sccmec --assembly GCA_000442355.1_SA16_genome_assembly_genomic.fna > output.csv
mkdir output_staphophia
#para varias cepas
staphopia-sccmec --assembly . > output_final.tsv

#run control strain
cd /media/gerald/DC604FBB604F9B62/sa_congreso/patric/genomes_sa_patric_southamerica/controls
staphopia-sccmec --assembly . > output2.tsv

#amrfinderplus
conda create --name amrfinderplus
conda activate amrfinderplus
conda install bioconda::ncbi-amrfinderplus

cd /media/gerald/DC604FBB604F9B62/sa_congreso/patric/genomes_sa_patric_southamerica
python -m resfinder -o armfinderplus -s "Staphylococcus aureus" -l 0.6 -t 0.8 --acquired --point -ifa *.fna

#Prueba de amrfinderplus
#para una cepa
cd /media/gerald/DC604FBB604F9B62/sa_congreso/patric/genomes_sa_patric_southamerica/armfinderplus


amrfinder -n ./*.fna --threads 4 --plus > out.amrfinder


#Para todas las cepas

amrfinder -n ./*.fna --threads 4 --plus > out.amrfinder

for assembly in *.fna
do
    base=$(basename $assembly .fna)
    # note that we use the --name option to add the "base" as the first
    # column of the output
    amrfinder -n $assembly --threads 4 --plus --name=$base > $base.amrfinder
done

head -1 $(ls *.amrfinder | head -1) > combined.tsv
grep -h -v 'Protein identifier' *.amrfinder >> combined.tsv

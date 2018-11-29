#!/usr/bin/env bash

export PATH=$PATH:/Users/derek_newberger/sratoolkit.2.9.2-mac64/bin

for SRA_number in $(cut -f 7 data/metadata/SraRunTable.txt | tail -n +2)
do
   fastq-dump -v $SRA_number -O data/raw_data
done

#fastqc data/raw_data/*.fastq --outdir=output/fastqc


#for file in data/raw_data/*.fastq
#do
#TrimmomaticSE -threads 2 -phred33 $file data/trimmed/$(basename -s .fastq $file).trim.fastq LEADING:5 TRAILING:5 SLIDINGWINDOW:8:25 MINLEN:150
#done

#for file in data/trimmed/*.trim.fastq
#do
#bioawk -c fastx '{print ">"$name"\n"$seq}' $file > data/billy_fasta/$(basename -s .trim.fastq $file).fasta
#done

#for file in data/billy_fasta/ERR1942297.fasta
#do
#blastn -db /blast-db/nt -num_threads 2 -outfmt '10 sscinames std' -out output/billytwin_csv/$(basename -s .fasta $file).csv -max_target_seqs 1 -negative_gilist /blast-db/2018-09-19_environmental_sequence.gi -query $file
#done
#!/usr/bin/env bash

export PATH=$PATH:/Users/derek_newberger/sratoolkit.2.9.2-mac64/bin

for SRA_number in $(cut -f 7 data/metadata/SraRunTable.txt | tail -n +2)
do
   fastq-dump -v "$SRA_number" -O data/raw_data
done
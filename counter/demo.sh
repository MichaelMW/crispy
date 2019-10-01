#!/bin/bash

### 1. how to determine flank sequence
## eg. in the excel sheet for oligo design, "FMR1 oligo-final 11162017.xlsx"
## the final sequence follows the formula:
## ="CTTGGAGAAAAGCCTTGTTT"&A2&"GTTTAgagacg"&F2&"cgtctcACACC"&C2&"GTTTTAGAGCTAGAAATAGCAAGTT"
## where A2 and C2 are fields for sgRNA1 and sgRNA2. 
## this corresponds to sgRNA1 and sgRNA2 field in the "oligo.tsv" provided to Crispy
## therefore, the flanking sequences are "GTTT" and "GTTT" for sgRNA1, which in the formula flank A2. 
## the flanking sequences are "CACC" and "GTTT" for sgRNA2, which in the formula flank C2.

### 2. call reads: fastq -> reads
## here depending on if the reads are single reads or paired reads, the command looks different. 
## for FMR1, the reads are paired. Use this:

inFile="demo1.fq.gz"
zcat < $inFile| awk 'NR%4==2' | ./counter.py -f1 "GTTT" -f2 "GTTT" -i oligos.tsv 

inFile="demo2.fq.gz"
zcat < $inFile| awk 'NR%4==2' | ./counter.py -f1 "CACC" -f2 "GTTT" -i oligos.tsv -S 7 -R 1



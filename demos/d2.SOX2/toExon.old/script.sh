#!/bin/sh

## get sgRNA in exon region

intersectBed -a oligos.tsv -b exon.mm9.bed -u > oligos.exon.tsv
./filter.py > reads.exon.tsv

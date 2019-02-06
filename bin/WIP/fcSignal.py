#!/usr/bin/env python
# encoding: utf-8


## this script input oligo tsv file from call.gRNA (edgeR) and output 
## fold change (FC) signal, using averaged enrichment in target bins. 
## Averaging FC in target bins gives high resolution than pvalue, 
## as has been practiced by various previous studies.

from sys import argv
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)


## input 
#pvalTsv=$1
#inRegion=$2
#inSgrna=$3

pvalTsv = "results.demo/cis_loose.sgRNA.tsv"
INSGRNA = "demos/d1.Yarui/oligos.tsv"
SGRNAKW = "test"
fcLimit = -2

## index gID2fc: sgRNA id -> fold change
gID2fc = {}
inFile = pvalTsv
with open(inFile) as f:
	header = l.readline()
	for l in f.readlines():
		ls = l.strip().split()
		if ls[-1] == SGRNAKW:
			gID = ls[0]
			fc = float(ls[1])
			gID2fc[gID] = fc
			
## get genomic location of INSGRNA
gID2locus = {}
chrms = set()
inFile = INSGRNA
with open(inFile) as f:
	for l in f.readlines():
		header = l.readline()
		ls = l.strip().split()
		chrm, start, end, gID = ls[:4]
		chrms.add(chrm)
		gID2locus[gID] = [chrm, start, end]

## 

# cutoff, use FC>0?

## get FC
tail -n+2 $pvalTsv | awk -v kw=$SGRNAKW -v OFS="\t" '{if($NF==kw && $2>0){print $1, $2}}' > tmp.1

## map to genome loci
awk '{print $4"\t"$1"_"$2"_"$3}' $INSGRNA > tmp.2
./bin/rp tmp.1 tmp.2 | tr "_" "\t" > tmp.3

## average in inregion
intersectBed -a $INREGION -b tmp.3 -wa -wb | ./bin/collapseBed -c 7 -o "mean" -q  > tmp.bedgraph

## filter:
## 1.with pval cutoff

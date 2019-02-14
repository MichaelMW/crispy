#!/bin/sh
set -ue

## this script is a post-RRA region filter
## it asks the question: does each RRA region have enough positive sgRNAs in it?
## "enough" is defined by either a raw number (eg. c=2) or by a percentage (eg. 0.1). 

## this script input 
## 1. positive sgRNA track
## 2. all sgRNA track
## 3. input region track (of RRA output)
## output 
## region track, each has 

## these regions are filtered and then feed to rra.R. 

### manual input
## input files
#INREGION="demos/d2.SOX2/regions.bed"
#INSGRNA="demos/d2.SOX2/oligos.tsv"
#OUTDIR="results.test"
#PREFIX="Sox2_all_G"
## input parameter
#percCut=0.5 ## percentage cutoff. 

### automated
## input files
INREGION=$1
INSGRNA=$2
OUTDIR=$3
PREFIX=$4
## input parameter
percCut=$5 ## percentage cutoff. 

## derived input files
sgrnaSignal="$OUTDIR/$PREFIX.sgRNA.bedgraph"
sgrnaSignalAll="$OUTDIR/$PREFIX.sgRNA_all.bedgraph"
regionSignal="$OUTDIR/$PREFIX.peak.bedgraph"

## derived output file
mkdir -p "$OUTDIR.percFilter"
regionStats="$OUTDIR.percFilter/$PREFIX.regionStats.$percCut.bedgraph"
goodRegion="$OUTDIR.percFilter/$PREFIX.goodRegion.$percCut.bed"
goodRegionSignal="$OUTDIR.percFilter/$PREFIX.goodPeak.$percCut.bed"

## put sgRNA in region
intersectBed -a $INREGION -b $sgrnaSignal -wa -wb | ./bin/collapseBed -c 7 -o len | awk '{print $1"_"$2"_"$3"\t"$4}' | sort > percFilter.tmp.nPos
intersectBed -a $INREGION -b $sgrnaSignalAll -wa -wb | ./bin/collapseBed -c 7 -o len | awk '{print $1"_"$2"_"$3"\t"$4}' | sort > percFilter.tmp.nAll

## get stats for each region
## output locus, nPos, nAll, nPos/nAll > $regionStats
join percFilter.tmp.nPos percFilter.tmp.nAll | tr "_" "\t" | tr " " "\t" | awk '{print $0"\t"$4/$5}' > $regionStats

## regions that passed percentage cutoff.
awk -v percCut=$percCut -v OFS="\t" '{if($5>=percCut){print $1,$2,$3}}' $regionStats > $goodRegion

## intersect the RRA signal with good regions. 
intersectBed -a $goodRegion -b $regionSignal -wa -wb | cut -f 4- > $goodRegionSignal

## clean up
rm percFilter.tmp.*


#!/bin/bash

## input files
# target region
inRegion=$1

# map oligoID to oligo region
inOligo=$2

# sgRNA statistics output from edgeR
inSgRna=$3

# pvalue cutoff used to define the positive sgRNAs. 
pvalCut=$4

## test input
#inRegion="../mmr_ATAC_EXON_regions_unique.bed"
#inOligo="../MMR_unique_oligos_03272019_1.tsv"
#inSgRna="MMR_CRISPRi_XR058_59(unsorted)_XR066_67(6tg_80)_chr2.sgRNA.tsv" 
#pvalCut=0.039894

## parameters
# keyword for sgRNA group
kw="test"

## process temp file
awk '{print($4"\t"$1"_"$2"_"$3)}' $inOligo > "id2loci.tmp"

## get all sgRNA
tail -n+2 $inSgRna | awk -v kw=$kw '{if($NF==kw){print}}' | awk -v pvalCut=$pvalCut -v OFS="\t" '{print $1,$2,-log($5),$6}' | ./bin/rp - id2loci.tmp | tr "_" "\t" > allSgRNA.bedgraph

## get positive sgRNA is defined by edgeR pvalue.
tail -n+2 $inSgRna | awk -v kw=$kw '{if($NF==kw){print}}' | awk -v pvalCut=$pvalCut -v OFS="\t" '{if($5<pvalCut && $2>0){print $1,$2,-log($5),$6}}' | ./bin/rp - id2loci.tmp | tr "_" "\t" > posSgRNA.bedgraph

## put sgRNAs to region
intersectBed -a $inRegion -b allSgRNA.bedgraph -wa -wb | ./bin/collapseBed -e 1,2,3 -c 7,8,9,9 -o mean,mean,set,len -q | ./bin/mySortBed > allSgRNA.region.tsv
intersectBed -a $inRegion -b posSgRNA.bedgraph -wa -wb | ./bin/collapseBed -e 1,2,3 -c 7,8,9,9 -o mean,mean,set,len -q | ./bin/mySortBed > posSgRNA.region.tsv

### merge all and pos to one file
header="chrm\tstart\tend\t#all\tavgFC_all\tavgLogP_all\t#pos\tavgFC_pos\tavgLogP_pos\tperc"
echo -e $header
./bin/myJoin allSgRNA.region.tsv posSgRNA.region.tsv -e 1,2,3 -c 7,4,5 -f "NA" | tail -n+2 | grep -v -w "NA" | awk '{printf "%s\t%s\t%s\t%s\t%.3f\t%.3f\t%s\t%.3f\t%.3f\t%.3f\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$7/$4}' | sort -k10,10nr


## remove temp files
rm id2loci.tmp allSgRNA.region.tsv posSgRNA.region.tsv allSgRNA.bedgraph posSgRNA.bedgraph


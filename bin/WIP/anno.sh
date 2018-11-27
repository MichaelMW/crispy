#!/bin/sh


########## work in progress. #########
# use this to annotate output peaks  #
######################################

set -ue

### annotate with ChromHMM

#inPeakSignal="cis_default.peak.bedgraph"

inPeakSignal=$1   ## 4 columns
inAnnoFile1="annotation/HMM.H1.bed"
inAnnoFile2="annotation/HMM.GM12878.bed"

intersectBed -a $inPeakSignal -b $inAnnoFile1 -wa -wb | bedtools groupby -c 8 -o "distinct" > "tmp1"
intersectBed -a $inPeakSignal -b $inAnnoFile2 -wa -wb | bedtools groupby -c 8 -o "distinct" > "tmp2"

header="chrm\tstart\tend\tHMM_H1\tHMM_GM12878"
echo $header
join  <(awk '{print $1"_"$2"_"$3,$4}' tmp1| sort) <(awk '{print $1"_"$2"_"$3,$4}' tmp2| sort) | awk '{split($1,loci,"_");print loci[1],loci[2],loci[3],$2,$3}' | tr " " "\t"

rm tmp1 tmp2

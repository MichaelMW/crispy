#!/bin/sh

sgrnaSignal1="results.SOX2/tiling_GFP.sgRNA.bedgraph"
sgrnaSignal2="results.SOX2/tiling_mCherry.sgRNA.bedgraph"
INREGION="demos/d2.SOX2/regions.bed"
## sgRNA signals -> target region signals
echo "######### 2. sgRNA signals -> target region signals ... #########"

cut -f1-3 $sgrnaSignal1 > tmp
cut -f4 $sgrnaSignal1 | ./bin/val2rank.py -nr | paste -d "\t" tmp -  > tmp.3.1
cut -f1-3 $sgrnaSignal2 > tmp
cut -f4 $sgrnaSignal2 | ./bin/val2rank.py -nr | paste -d "\t" tmp -  > tmp.3.2

cat tmp.3.1 tmp.3.2 > tmp.3
# put sgRNA ranks into region bins
intersectBed -a $INREGION -b tmp.3 -wa -wb | ./bin/collapseBed -c 7 -o "list" -q  > tmp.4
./bin/rra.R --inFile="tmp.4" --outFile="tmp.5" --method="RRA"
# convert and filter pvalues of region bins to signals.
regionSignal="results.SOX2/tiling_merged.region.bedgraph"
tail -n+2 tmp.5 | tr "_" "\t" | awk -v RRACUTOFF=1 -v OFS="\t" '{if($4<RRACUTOFF){$4=-log($4);print}}' | ./bin/mySortBed > $regionSignal


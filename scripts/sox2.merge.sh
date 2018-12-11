#!/bin/sh

#### Sox2: Merging GFP and mCherry channels. 
### Strategy 2: merge rankings of sgRNA pvalues from GPF and from mcherry -> RRA

./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o results.SOX2 \
			-n 1 \
            -p "tiling_GFP" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GpMn,S2.GpMn"\
			-t "png"
			

./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o results.SOX2 \
			-n 1 \
            -p "tiling_mCherry" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GnMp,S2.GnMp"\
			-t "png"

### manaully merge sgRNAs from tiling_GFP and tiling_mCherry

sgrnaSignal1="results.SOX2/tiling_GFP.sgRNA.bedgraph"
sgrnaSignal2="results.SOX2/tiling_mCherry.sgRNA.bedgraph"
INREGION="demos/d2.SOX2/regions.bed"
## sgRNA signals -> target region signals
echo "######### 2. sgRNA signals -> target region signals ... #########"
paste -d"\t" <(cut -f1-3 $sgrnaSignal1) <(cut -f4 $sgrnaSignal1| ./bin/val2rank.py -nr) > tmp.3.1
paste -d"\t" <(cut -f1-3 $sgrnaSignal2) <(cut -f4 $sgrnaSignal2| ./bin/val2rank.py -nr) > tmp.3.2
cat tmp.3.1 tmp.3.2 > tmp.3
# put sgRNA ranks into region bins
intersectBed -a $INREGION -b tmp.3 -wa -wb | ./bin/collapseBed -c 7 -o "list" -q  > tmp.4
./bin/rra.R --inFile="tmp.4" --outFile="tmp.5" --method="RRA"
# convert and filter pvalues of region bins to signals.
regionSignal="results.SOX2/tiling_merged.region.bedgraph"
tail -n+2 tmp.5 | tr "_" "\t" | awk -v RRACUTOFF=1 -v OFS="\t" '{if($4<RRACUTOFF){$4=-log($4);print}}' | ./bin/mySortBed > $regionSignal

## try macs2
subtractBed -a $INREGION -b $regionSignal | awk '{print $0"\t"0}' | cat - $regionSignal | ./bin/mySortBed > tmp; mv tmp $regionSignal
## call peaks
peakSignal="results.SOX2/tiling_merged.peak.bedgraph"
MINLEN=20
MAXGAP=10000
PEAKCUTOFF=1
macs2 bdgpeakcall -i $regionSignal -l $MINLEN -g $MAXGAP -c $PEAKCUTOFF -o /dev/stdout | tail -n+2 | cut -f1-3,5 > $peakSignal
areacutoff=20
awk -v areacutoff=$areacutoff '{if($4>areacutoff) print}' $peakSignal > tmp; mv tmp $peakSignal

## cleanup
rm tmp.*

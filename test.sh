#!/bin/bash


###### notes on setting three key parameters:
## -n is fdr for sgRNA cut off. Yarui Used 0.05 in the nature method paper. Increase this value up to 1 to get more peaks.
## -a is rra peaks. Crispy default is use 1 to keep all (recommended). 
## -c higher is more stringent. 3 means 0.001 peaks filter in macs2. Lower this value down to 0+ to get more peaks. 


## this is the current script for testing new features. 

rm -rf "results.test"

# default CPM filter cutoff
./crispy.sh -i demos/d1.Yarui/data.tsv \
            -r demos/d1.Yarui/regions.bed \
            -s demos/d1.Yarui/oligos.tsv \
            -o "results.test" \
            -p "cis_enriched" \
            -b "ctr1,ctr2" \
            -f "cis1,cis2,cis3,cis4,cis5" \
			-t "png" 

# sox2
./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o "results.test" \
            -p "tiling_4fg" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GnMp,S1.GpMn,S2.GnMp,S2.GpMn" \
			-t "png"

./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o "results.test" \
            -p "tiling_gfp" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GpMn,S2.GpMn" \
			-t "png"

./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o "results.test" \
            -p "tiling_mcherry" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GnMp,S2.GnMp" \
			-t "png"


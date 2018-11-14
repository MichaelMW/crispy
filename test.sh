#!/bin/bash


###### notes on setting three key parameters:
## -n is fdr for sgRNA cut off. Yarui Used 0.05 in the nature method paper. Increase this value up to 1 to get more peaks.
## -a is rra peaks. Crispy default is use 1 to keep all (recommended). 
## -c higher is more stringent. 3 means 0.001 peaks filter in macs2. Lower this value down to 0+ to get more peaks. 


## this is the current script for testing new features. 

rm -rf "results.Yarui"

# loose CPM filter cutoff
./crispy.sh -i demos/d1.Yarui/data.tsv \
            -r demos/d1.Yarui/regions.bed \
            -s demos/d1.Yarui/oligos.tsv \
            -o results.Yarui \
            -p "cis_loose" \
            -b "ctr1,ctr2" \
            -f "cis1,cis2,cis3,cis4,cis5" \
			-q "cis1,cis2,cis3,cis4,cis5;ctr1,ctr2;high1,high2,high3" \
			-u 0 \
			-v 0 \
			-t "png" 

# default CPM filter cutoff
./crispy.sh -i demos/d1.Yarui/data.tsv \
            -r demos/d1.Yarui/regions.bed \
            -s demos/d1.Yarui/oligos.tsv \
            -o results.Yarui \
            -p "cis_default" \
            -b "ctr1,ctr2" \
            -f "cis1,cis2,cis3,cis4,cis5" \
			-q "cis1,cis2,cis3,cis4,cis5;ctr1,ctr2;high1,high2,high3" \
			-t "png" 

# extreme CPM filter cutoff
./crispy.sh -i demos/d1.Yarui/data.tsv \
            -r demos/d1.Yarui/regions.bed \
            -s demos/d1.Yarui/oligos.tsv \
            -o results.Yarui \
            -p "cis_extreme" \
            -b "ctr1,ctr2" \
            -f "cis1,cis2,cis3,cis4,cis5" \
			-q "cis1,cis2,cis3,cis4,cis5;ctr1,ctr2;high1,high2,high3" \
			-u 50 \
			-v 1 \
			-t "png" 

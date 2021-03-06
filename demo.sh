#!/bin/bash


###### notes on setting three key parameters:
## -n is fdr for sgRNA cut off. Yarui Used 0.05 in the nature method paper. Increase this value up to 1 to get more peaks.
## -a is rra peaks. Crispy default is use 1 to keep all (recommended). 
## -c higher is more stringent. 3 means 0.001 peaks filter in macs2. Lower this value down to 0+ to get more peaks. 


rm -rf "results.Yarui"

## demo1, Yarui, 2017 nature method paper. 
# loose pvalue cut off
./crispy.sh -i demos/d1.Yarui/data.tsv \
            -r demos/d1.Yarui/regions.bed \
            -s demos/d1.Yarui/oligos.tsv \
            -o results.demo \
            -p "cis_loose" \
            -b "ctr1,ctr2" \
            -f "cis1,cis2,cis3,cis4,cis5" \
           	-n 0.1 \
           	-c 1\
			-t png 

## demo2
# stringent pvalue cutoff
./crispy.sh -i demos/d1.Yarui/data.tsv \
            -r demos/d1.Yarui/regions.bed \
            -s demos/d1.Yarui/oligos.tsv \
            -o results.demo \
            -p "cis_stringent" \
            -b "ctr1,ctr2" \
            -f "cis1,cis2,cis3,cis4,cis5" \
           	-n 0.05 \
           	-c 3


## demo2.1
# stringent pvalue cutoff, with qnorm
./crispy.sh -i demos/d1.Yarui/data.tsv \
            -r demos/d1.Yarui/regions.bed \
            -s demos/d1.Yarui/oligos.tsv \
            -o results.demo \
            -p "cis_stringent_qnorm" \
            -b "ctr1,ctr2" \
            -f "cis1,cis2,cis3,cis4,cis5" \
			-q 1 \
            -n 0.05 \
            -c 3


## demo3
# example to show unconventional usage. 
# no replicate mode; call only depleted sgRNA; use method=min instead of RRA
./crispy.sh -i demos/d1.Yarui/data.tsv \
            -r demos/d1.Yarui/regions.bed \
            -s demos/d1.Yarui/oligos.tsv \
            -o results.demo \
            -p "cis_noRep_depletion" \
            -b "ctr1" \
            -f "cis1" \
			-d -1 \
			-m "min"


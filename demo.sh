#!/bin/bash


###### notes on setting three key parameters:
## -n is fdr for sgRNA cut off. Yarui Used 0.05 in the nature method paper. Increase this value up to 1 to get more peaks.
## -a is rra peaks. Crispy default is use 1 to keep all (recommended). 
## -c higher is more stringent. 3 means 0.001 peaks filter in macs2. Lower this value down to 0+ to get more peaks. 

## demo1, Yarui, 2017 nature method paper. 
# loose
./crispy.sh -i demos/d1.Yarui/data.tsv \
            -r demos/d1.Yarui/regions.bed \
            -s demos/d1.Yarui/oligos.tsv \
            -o results.Yarui \
            -p "cis_loose" \
            -b "ctr1,ctr2" \
            -f "cis1,cis2,cis3,cis4,cis5" \
           	-n 0.1 \
           	-c 1

## demo2
# stringent
./crispy.sh -i demos/d1.Yarui/data.tsv \
            -r demos/d1.Yarui/regions.bed \
            -s demos/d1.Yarui/oligos.tsv \
            -o results.Yarui \
            -p "cis_stringent" \
            -b "ctr1,ctr2" \
            -f "cis1,cis2,cis3,cis4,cis5" \
           	-n 0.05 \
           	-c 3

## demo3
# no replicate mode
./crispy.sh -i demos/d1.Yarui/data.tsv \
            -r demos/d1.Yarui/regions.bed \
            -s demos/d1.Yarui/oligos.tsv \
            -o results.Yarui \
            -p "cis_stringent_noRep" \
            -b "ctr1" \
            -f "cis1" \
            -n 0.05 \
            -c 3


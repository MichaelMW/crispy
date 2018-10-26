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


## demo2, Sox2 from Xingjie & Xiaoyu
# using tiling 50bp bins as regions input. 
./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o results.SOX2 \
            -p "tiling" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GnMp,S2.GpMn"


# using original DHS, histone peaks targets as regions input. 
./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/targetRegions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o results.SOX2 \
            -p "target" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GnMp,S2.GpMn"


## demo3, FMR1, AFF2 double tag from Xingjie. 
## XR016, 17 contains FMR1 low AFF2 high signals
## XR019, 23 contains FMR1 high AFF2 low signals.
# screen for FMR1-specific peaks
./crispy.sh -i demos/d3.FMR1AFF2/reads.tsv \
            -r demos/d3.FMR1AFF2/regions.bed \
            -s demos/d3.FMR1AFF2/oligos.tsv \
            -o results.FMR1AFF2 \
            -p "FMR1_specific_sgRNA" \
            -b "XR016,XR017" \
            -f "XR020,XR024"

# screen for AFF2-specific peaks
./crispy.sh -i demos/d3.FMR1AFF2/reads.tsv \
            -r demos/d3.FMR1AFF2/regions.bed \
            -s demos/d3.FMR1AFF2/oligos.tsv \
            -o results.FMR1AFF2 \
            -p "AFF2_specific_sgRNA" \
            -b "XR019,XR023" \
            -f "XR020,XR024"


## demo4, SIN3a from Xingjie. Default parameters. 
# default parameters. 
./crispy.sh -i demos/d4.SIN3A/reads.tsv \
            -r demos/d4.SIN3A/regions.bed \
            -s demos/d4.SIN3A/oligos.tsv \
            -o results.SIN3A \
            -p "cis_default" \
            -b "XR029,XR034" \
            -f "XR031,XR036,XR032,XR037"


# extremely loose parameters.
./crispy.sh -i demos/d4.SIN3A/reads.tsv \
            -r demos/d4.SIN3A/regions.bed \
            -s demos/d4.SIN3A/oligos.tsv \
            -o results.SIN3A \
            -p "cis_loose" \
            -b "XR029,XR034" \
            -f "XR031,XR036,XR032,XR037" \
			-n 1 \
			-c 0.1



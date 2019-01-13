#!/bin/bash


###### notes on setting three key parameters:
## -n is fdr for sgRNA cut off. Yarui Used 0.05 in the nature method paper. Increase this value up to 1 to get more peaks.
## -a is rra peaks. Crispy default is use 1 to keep all (recommended). 
## -c higher is more stringent. 3 means 0.001 peaks filter in macs2. Lower this value down to 0+ to get more peaks. 


## this is the current script for testing new features. 

rm -rf "results.test"

## default CPM filter cutoff
#./crispy.sh -i demos/d1.Yarui/data.tsv \
#            -r demos/d1.Yarui/regions.bed \
#            -s demos/d1.Yarui/oligos.tsv \
#            -o "results.test" \
#            -p "cis_enriched" \
#            -b "ctr1,ctr2" \
#            -f "cis1,cis2,cis3,cis4,cis5" \
#			-t "png" 

#sin3a
./crispy.sh -i demos/d4.SIN3A/sin3a_ipsc_reads_3L_3R_12262018.tsv \
            -r demos/d4.SIN3A/sin3a_regions.bed \
            -s demos/d4.SIN3A/sin3a_oligos.tsv \
            -o "results.test" \
            -p "cis_29_34" \
            -b "XR029,XR034" \
			-f "XR032,XR037,XR031,XR036" \
			-q "XR029,XR034;XR030,XR035;XR032,XR037,XR031,XR036" \
			-n 0.025 \
			-c 2 \
            -t "png"


# sox2
# ./crispy.sh -i demos/d2.SOX2/reads.tsv \
#             -r demos/d2.SOX2/regions.bed \
#             -s demos/d2.SOX2/oligos.tsv \
#             -o "results.test" \
#             -p "tiling_4fg" \
#             -b "S1.Unsorted,S2.Unsorted" \
#             -f "S1.GnMp,S1.GpMn,S2.GnMp,S2.GpMn" \
# 			-t "png" \
# 			-n 0.13

# ./crispy.sh -i demos/d2.SOX2/reads.tsv \
#             -r demos/d2.SOX2/regions.bed \
#             -s demos/d2.SOX2/oligos.tsv \
#             -o "results.test" \
#             -p "tiling_gfp" \
#             -b "S1.Unsorted,S2.Unsorted" \
#             -f "S1.GpMn,S2.GpMn" \
# 			-t "png" \
# 			-n 0.12

# ./crispy.sh -i demos/d2.SOX2/reads.tsv \
#             -r demos/d2.SOX2/regions.bed \
#             -s demos/d2.SOX2/oligos.tsv \
#             -o "results.test" \
#             -p "tiling_mcherry" \
#             -b "S1.Unsorted,S2.Unsorted" \
#             -f "S1.GnMp,S2.GnMp" \
# 			-t "png" \
# 			-n 0.12


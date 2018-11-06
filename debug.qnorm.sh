#!/bin/sh


## demo 2.1
# use all 4 single channel results as foreground
./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o results.SOX2 \
            -p "tiling_4fg" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GnMp,S1.GpMn,S2.GnMp,S2.GpMn"

# use qnorm, within fgs and within bgs individually.
./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o results.SOX2 \
            -p "tiling_4fg_qnorm" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GnMp,S1.GpMn,S2.GnMp,S2.GpMn" \
            -q 1

# use qnorm on all 4 fg, but choose gfp channel as fg only.
./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o results.SOX2 \
            -p "tiling_4fg_qnorm4fg_gfp" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GpMn,S2.GpMn" \
            -q "S1.GnMp,S1.GpMn,S2.GnMp,S2.GpMn"


# use qnorm on all 4 fg, but choose mcherry channel as fg only.
./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o results.SOX2 \
            -p "tiling_4fg_qnorm4fg_mcherry" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GnMp,S2.GnMp" \
            -q "S1.GnMp,S1.GpMn,S2.GnMp,S2.GpMn"


#!/bin/sh

# using tiling 50bp bins as regions input.
./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o results.SOX2 \
            -p "tiling.GnMp" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GnMp,S2.GnMp"

./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o results.SOX2 \
            -p "tiling.GpMn" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GpMn,S2.GpMn"

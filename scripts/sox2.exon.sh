#!/bin/sh

### filter reads by exon, check out the pvalues
./crispy.sh -i demos/d2.SOX2/reads.exon.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o results.SOX2 \
            -p "tiling_GFP_exon" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GpMn,S2.GpMn"\
            -t "png"


./crispy.sh -i demos/d2.SOX2/reads.exon.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o results.SOX2 \
            -p "tiling_mCherry_exon" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GnMp,S2.GnMp"\
            -t "png"


./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o results.SOX2 \
            -p "tiling_GFP" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GpMn,S2.GpMn"\
            -t "png"


./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o results.SOX2 \
            -p "tiling_mCherry" \
            -b "S1.Unsorted,S2.Unsorted" \
            -f "S1.GnMp,S2.GnMp"\
            -t "png"


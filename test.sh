#!/bin/sh

## test the new fold change signal. 
## new -e option controls the lower limit for FC plot. default: -0.5
## interestingly, if you set very loose pvalue cutoff (-n 1, which includes all sgRNAs)
## you got much better resolution. (because you included more sgRNAs.)
## this is not too bad, since RRA controls the overall stringency of the test. 
./crispy.sh -i demos/d1.Yarui/data.tsv \
            -r demos/d1.Yarui/regions.bed \
            -s demos/d1.Yarui/oligos.tsv \
            -o results.test \
            -p "yarui_all" \
            -b "ctr1,ctr2" \
            -f "cis1,cis2,cis3,cis4,cis5" \
            -n 1 \
            -c 2\
            -e -0.5\
            -t png

## how about using this n=1 for analysis?
## Dppa2
./crispy.sh -i demos/d6.Dppa2/reads.tsv \
            -r demos/d6.Dppa2/region.bed \
            -s demos/d6.Dppa2/oligos.tsv \
            -o results.test \
            -p "Dppa2_all" \
            -b "U1,U2,Ck1,Ck2" \
            -f "Mn1,Mn2" \
            -n 1 \
            -c 2\
            -t png

## Esrrb
./crispy.sh -i demos/d5.Esrrb/reads.tsv \
            -r demos/d5.Esrrb/regions.tile.bed \
            -s demos/d5.Esrrb/oligos.tsv \
            -o results.test \
            -p "Esrrb_all_M" \
            -b "U1,U2" \
            -f "M1,M2" \
            -n 1 \
            -c 1 \
            -t png

./crispy.sh -i demos/d5.Esrrb/reads.tsv \
            -r demos/d5.Esrrb/regions.tile.bed \
            -s demos/d5.Esrrb/oligos.tsv \
            -o results.test \
            -p "Esrrb_all_G" \
            -b "U1,U2" \
            -f "G1,G2" \
            -n 1 \
            -c 1 \
            -t png

## Sox2
./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o results.test \
            -p "Sox2_all_M" \
            -b "U1,U2" \
            -f "M1,M2" \
            -n 1 \
            -c 2 \
            -t png

./crispy.sh -i demos/d2.SOX2/reads.tsv \
            -r demos/d2.SOX2/regions.bed \
            -s demos/d2.SOX2/oligos.tsv \
            -o results.test \
            -p "Sox2_all_G" \
            -b "U1,U2" \
            -f "G1,G2" \
            -n 1 \
            -c 2 \
            -t png


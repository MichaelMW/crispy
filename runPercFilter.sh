#!/bin/sh

######## for Xiaoyu ########
######## a new script "percFilter.sh" can be used on top of existing crispy run 
######## to get region bins that passed a percentage cutoff
######## there are three output files
######## "regionstat" give you each region nPos, nAll, and nPos/nAll 
######## "goodregion" give you region with nPos/nAll more than specified percCutOff 
######## "goodpeak" give you region RRA signal within "goodregion"

######## here is a demo 


####### Sox2 #######

#### mcherry ####
### crispy run
INREAD="demos/d2.SOX2/reads.tsv"
INREGION="demos/d2.SOX2/regions.bed"
INSGRNA="demos/d2.SOX2/oligos.tsv"
DIR="results.Sox2"
PREFIX="Sox2_all_M"
./crispy.sh -i $INREAD \
            -r $INREGION \
            -s $INSGRNA \
            -o $DIR \
            -p $PREFIX \
            -b "U1,U2" \
            -f "M1,M2" \
            -n 0.05 \
            -c 0

### percentage filter
percCut="0.5" # eg. 0.5 means a bin has to have over 50% good sgRNA
## run script
./percFilter.sh $INREGION $INSGRNA $DIR $PREFIX $percCut
## see results in $DIR.percFilter

#### GFP ####
### crispy run
INREAD="demos/d2.SOX2/reads.tsv"
INREGION="demos/d2.SOX2/regions.bed"
INSGRNA="demos/d2.SOX2/oligos.tsv"
DIR="results.Sox2"
PREFIX="Sox2_all_G"
./crispy.sh -i $INREAD \
            -r $INREGION \
            -s $INSGRNA \
            -o $DIR \
            -p $PREFIX \
            -b "U1,U2" \
            -f "G1,G2" \
            -n 0.05 \
            -c 0

### percentage filter
percCut="0.5" # eg. 0.5 means a bin has to have over 50% good sgRNA
## run script
./percFilter.sh $INREGION $INSGRNA $DIR $PREFIX $percCut
## see results in $DIR.percFilter

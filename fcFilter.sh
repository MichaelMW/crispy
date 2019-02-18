#!/bin/sh


### this script input 
INSGRNA=$1
DIR=$2
PREFIX=$3
fcPvalCut=$4  	# eg. "0.05" # a sgRNA need to have <=$fcPvalCut to be considered for FC calculation. 

### hard coded input arguments. Don't change. 
SGRNAKW="test" 	# keywords indicating sgRNA is the probe of test.
FCLIMIT="-0.5" 	# sgRNA with FC < FCLIMIT will be changed to FC = FCLIMIT. # -0.5 is the default from the previous science paper. 
fcNCut=1 		# eg. "2" # a region need to have >= $fcNCut to be considered for FC calculation. Default is 1 to consider all regions. Make this input parameter later. 

### derived arguments
pvalTsv=$DIR/$PREFIX.sgRNA.tsv
fcSignalFiltered=$DIR/$PREFIX.fc.$fcPvalCut.bedgraph

### original code to get FC
# get FC
### at this point, filter based on pvalu  <----
tail -n+2 $pvalTsv | awk -v kw=$SGRNAKW -v OFS="\t" -v fcLimit=$FCLIMIT -v fcPvalCut=$fcPvalCut '{if($NF==kw && $5<=fcPvalCut){FC=($2>fcLimit?$2:fcLimit); print $1, FC}}' > tmp.1.fc

# map to genome loci
awk '{print $4"\t"$1"_"$2"_"$3}' $INSGRNA > tmp.2.fc
./bin/rp tmp.1.fc tmp.2.fc | tr "_" "\t" > tmp.3.fc
# get the best possible resolution by moving between boundaries of $INSGRNA
tail -n+2 $INSGRNA | ./bin/overlap2cont.py > tmp.4.fc

### at this point, filter based on number of counts <-----
# mean signal in region
intersectBed -a tmp.4.fc -b tmp.3.fc -wa -wb | ./bin/collapseBed -c 7,7 -o "mean,len" -q | awk -v OFS="\t" -v fcNCut=$fcNCut '{if($5>=fcNCut){print $1,$2,$3,$4}}' > $fcSignalFiltered

## clean space
rm tmp.{1,2,3,4}*

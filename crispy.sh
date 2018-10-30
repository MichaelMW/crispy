#!/bin/bash
set -ue

## check for prerequisits
# bedtools
# macs2

command -v macs2 >/dev/null 2>&1 || { echo "macs2 not found. Try this to install: 'pip install -U --no-deps MACS2'"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "bedtools not found. Try this page to install: 'https://bedtools.readthedocs.io/en/latest/content/installation.html'"; exit 1;}

## Default parameters
OUTDIR="crispyOUT"
PREFIX="results"
NBCUTOFF=0.05
DIRECTION=0
RRACUTOFF=1
PEAKCUTOFF=2
MINLEN=50
MAXGAP=200

## Hardcoded paramters
SGRNAKW="test"  ## keywords used to grep from sgRNA for the designed ones. 

## Get arguments
while getopts ":i:r:s:o:p:b:f:n:d:a:c:l:g" opt;
do
	case "$opt" in
		i) INREAD=$OPTARG;; 	# input read number for each sgRNA. (id name should match $INSGRNA)
		r) INREGION=$OPTARG;; 	# input user defined regions to output peaks. 
		s) INSGRNA=$OPTARG;; 	# input sgRNA information. requires chrom, start, end, sgRNAid(match $INREAD file)
		o) OUTDIR=$OPTARG;; 	# name of output dir
		p) PREFIX=$OPTARG;; 	# prefix for experiment output
		b) BG=$OPTARG;; 		# background column names. eg. "unsortedRep1,unsortedRep2"
		f) FG=$OPTARG;; 		# foreground column names. eg. "gpf+mcherry-,gpf-mcherry+"
		n) NBCUTOFF=$OPTARG;;  	# pval cutoff for sgRNA using negative bionomial test. Default=0.05.
		d) DIRECTION=$OPTARG;; 	# set to 1 to get only enriched sgRNA in fg. -1 for only depleted sgRNA. 0 for both. Default=0
		a) RRACUTOFF=$OPTARG;;  # pval cutoff for crest-seq signal in defined region using robust ranking aggregation test. Default=1.
		c) PEAKCUTOFF=$OPTARG;; # pval cutoff for final peaks from macs2. Default=2. (2 means 1e-2, ie, 0.01)
		l) MINLEN=$OPTARG;; 	# min length of final peaks from macs2. Default=50.
		g) MAXGAP=$OPTARG;; 	# max gap between signals to call final peaks from macs2. Default=200.
		\?) echo "Invalid option: -$OPTARG" 1>&2; exit 1;;
	esac
done


## Check arguments


## start run
echo "CRISPY started runing ..."
echo "./crispy $@"

## reads -> sgRNA
echo "######### read counts -> sgRNA signals ... #########"
./bin/call.gRNA.R --inFile=$INREAD --bg="$BG" --fg="$FG" --outDir="$OUTDIR" --prefix="$PREFIX"

## sgRNA -> sgRNA signals in loci
#echo "######### getting sgRNA signal from pvalues ... #########"
pvalTsv="$OUTDIR/$PREFIX.pvalues.tsv"  # from last step
sgrnaSignal="$OUTDIR/$PREFIX.pvalues.bedgraph" # output of this step, good for visualizing all sgRNA signals.
if [ "$DIRECTION" -eq 1 ]; then
	echo "DIRECTION==1, using only enriched sgRNA ..."
	tail -n+2 $pvalTsv | awk -v kw=$SGRNAKW -v nbcutoff=$NBCUTOFF -v OFS="\t" '{if($NF==kw && $5<nbcutoff && $2>0){print $1, -log($5)}}' > tmp1
elif [ "$DIRECTION" -eq -1 ]; then
	echo "DIRECTION==-1, using only depleted sgRNA ..."
	tail -n+2 $pvalTsv | awk -v kw=$SGRNAKW -v nbcutoff=$NBCUTOFF -v OFS="\t" '{if($NF==kw && $5<nbcutoff && $2<0){print $1, log($5)}}' > tmp1
else
	echo "DIRECTION==0, using sgRNAs in both directions ..."
	tail -n+2 $pvalTsv | awk -v kw=$SGRNAKW -v nbcutoff=$NBCUTOFF -v OFS="\t" '{if($NF==kw && $5<nbcutoff){print $1, -log($5) * (($2>0)-0.5)*2}}' > tmp1
fi
awk '{print $4"\t"$1"_"$2"_"$3}' $INSGRNA > tmp2
./bin/rp tmp1 tmp2 | tr "_" "\t" > $sgrnaSignal
rm -rf tmp1 tmp2
echo

## sgRNA signals -> ranks
echo "######### sgRNA signals -> target region signals ... #########"
#echo "######### converting sgRNA signals to ranks ... #########"
rankBed="$OUTDIR/$PREFIX.rank.bed" #
paste -d"\t" <(cut -f1-3 $sgrnaSignal) <(cut -f4 $sgrnaSignal| ./bin/val2rank.py -n) > $rankBed

# ranks -> put in defined regions
#echo "######### intersecting sgRNA ranks with target regions ... #########"
rankBinBed="$OUTDIR/$PREFIX.rankBin.bed"
intersectBed -a $INREGION -b $rankBed -wa -wb | ./bin/collapseBed -c 7 -o "list" -q  > $rankBinBed

# rra in defined regions
rraFile="$OUTDIR/$PREFIX.rra.tsv"
./bin/rra.R --inFile=$rankBinBed --outFile=$rraFile

# convert and filter rankFile to rankBed 
rraBed="$OUTDIR/$PREFIX.rra.bedgraph"
tail -n+2 $rraFile | tr "_" "\t" | awk -v RRACUTOFF=$RRACUTOFF -v OFS="\t" '{if($4<RRACUTOFF){$4=-log($4); print}}' | ./bin/mySortBed > $rraBed
echo 

## macs2 for peak smoothing
echo "######### target region signals -> peak smoothing ... #########"
## fill 0s for macs2
subtractBed -a $INREGION -b $rraBed | awk '{print $0"\t"0}' | cat - $rraBed | ./bin/mySortBed > tmp; mv tmp $rraBed
## call peaks
peakFile="$OUTDIR/$PREFIX.peaks.bedgraph"
macs2 bdgpeakcall -i $rraBed -l $MINLEN -g $MAXGAP -c $PEAKCUTOFF -o /dev/stdout | tail -n+2 | cut -f1-3,10 > $peakFile
echo "Crispy Done!"
echo
echo
echo

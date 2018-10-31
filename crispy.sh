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
DIRECTION=1
METHOD="RRA"
RRACUTOFF=1
PEAKCUTOFF=2
MINLEN=50
MAXGAP=200

## Hardcoded paramters
SGRNAKW="test"  ## keywords used to grep from sgRNA for the designed ones. 

## Get arguments
while getopts ":i:r:s:o:p:b:f:n:d:m:a:c:l:g" opt;
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
		d) DIRECTION=$OPTARG;; 	# set to 1 to get only enriched sgRNA in fg. -1 for only depleted sgRNA. Default=1.
		m) METHOD=$OPTARG;; 	# method used to call aggregate pvalue from pooled sgRNA pvalues. Choose from [RRA,min,geom.mean,median,stuart]. Default=RRA
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
echo "######### 1. read counts -> sgRNA signals ... #########"
./bin/call.gRNA.R --inFile=$INREAD --bg="$BG" --fg="$FG" --outDir="$OUTDIR" --prefix="$PREFIX"

## sgRNA -> sgRNA signals in loci
#echo "######### getting sgRNA signal from pvalues ... #########"
pvalTsv="$OUTDIR/$PREFIX.sgRNA.tsv"  # from last step
sgrnaSignal="$OUTDIR/$PREFIX.sgRNA.bedgraph" # output of this step, good for visualizing all sgRNA signals.
if [ "$DIRECTION" -eq -1 ]; then
	echo "DIRECTION==-1, using only depleted sgRNA ..."
	tail -n+2 $pvalTsv | awk -v kw=$SGRNAKW -v nbcutoff=$NBCUTOFF -v OFS="\t" '{if($NF==kw && $5<nbcutoff && $2<0){print $1, -log($5)}}' > tmp.1
else
	echo "DIRECTION==1, using only enriched sgRNA ..."
	tail -n+2 $pvalTsv | awk -v kw=$SGRNAKW -v nbcutoff=$NBCUTOFF -v OFS="\t" '{if($NF==kw && $5<nbcutoff && $2>0){print $1, -log($5)}}' > tmp.1
fi
awk '{print $4"\t"$1"_"$2"_"$3}' $INSGRNA > tmp.2
./bin/rp tmp.1 tmp.2 | tr "_" "\t" > $sgrnaSignal
echo

## sgRNA signals -> target region signals
echo "######### 2. sgRNA signals -> target region signals ... #########"
paste -d"\t" <(cut -f1-3 $sgrnaSignal) <(cut -f4 $sgrnaSignal| ./bin/val2rank.py -nr) > tmp.3
# put sgRNA ranks into region bins
intersectBed -a $INREGION -b tmp.3 -wa -wb | ./bin/collapseBed -c 7 -o "list" -q  > tmp.4
./bin/rra.R --inFile="tmp.4" --outFile="tmp.5" --method=$METHOD
# convert and filter pvalues of region bins to signals. 
regionSignal="$OUTDIR/$PREFIX.region.bedgraph"
tail -n+2 tmp.5 | tr "_" "\t" | awk -v RRACUTOFF=$RRACUTOFF -v OFS="\t" '{if($4<RRACUTOFF){$4=-log($4); print}}' | ./bin/mySortBed > $regionSignal
echo 

## macs2 for peak smoothing
echo "######### 3. target region signals -> peak smoothing ... #########"
## fill 0s for macs2
subtractBed -a $INREGION -b $regionSignal | awk '{print $0"\t"0}' | cat - $regionSignal | ./bin/mySortBed > tmp; mv tmp $regionSignal
## call peaks
peakSignal="$OUTDIR/$PREFIX.peak.bedgraph"
macs2 bdgpeakcall -i $regionSignal -l $MINLEN -g $MAXGAP -c $PEAKCUTOFF -o /dev/stdout | tail -n+2 | cut -f1-3,10 > $peakSignal
rm tmp.{1,2,3,4,5}
echo -e "Crispy Done!\n\n\n"



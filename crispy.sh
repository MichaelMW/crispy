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
QNORM="0"
MINCPM="5"
MINCPMRATIO="0.5"
PLOTFORMAT="PDF"
NBCUTOFF=0.05
DIRECTION=1
METHOD="RRA"
RRACUTOFF=1
MINSGRNA=2
MINLEN=1000

#PEAKCUTOFF=2
#MINLEN=50
#MAXGAP=200

## Hardcoded paramters
SGRNAKW="test"  ## keywords used to grep from sgRNA for the designed ones. 

## Get arguments
while getopts ":i:r:s:o:p:b:f:q:t:u:v:n:d:m:a:c:l" opt;
do
	case "$opt" in
		i) INREAD=$OPTARG;; 	# input read number for each sgRNA. (id name should match $INSGRNA)
		r) INREGION=$OPTARG;; 	# input user defined regions to output peaks. 
		s) INSGRNA=$OPTARG;; 	# input sgRNA information. requires chrom, start, end, sgRNAid(match $INREAD file)
		o) OUTDIR=$OPTARG;; 	# name of output dir
		p) PREFIX=$OPTARG;; 	# prefix for experiment output
		b) BG=$OPTARG;; 		# background column names. eg. "unsortedRep1,unsortedRep2"
		f) FG=$OPTARG;; 		# foreground column names. eg. "gpf+mcherry-,gpf-mcherry+"
		q) QNORM=$OPTARG;;		# quantile normalization of input reads among fgs and bgs, respectively. 0 or 1. default 0. You can also specify by using colnames, "," and ";". eg. "-q cis1,cis2;ctr1,ctr2;high1,high2"
		t) PLOTFORMAT=$OPTARG;;	# format of qc plot. default="pdf". Use png for faster loading.
		u) MINCPM=$OPTARG;;		# minimal read count cutoff for sgRNA to be deemed as express. Use 0 to disable this filter. Default:5
		v) MINCPMRATIO=$OPTARG;;# minimal ratio of FGs or BGs that have sgRNA read count higher than MINCPM. Use 0 to disable this filter. Default:0.5
		n) NBCUTOFF=$OPTARG;;  	# pval cutoff for sgRNA using negative bionomial test. Default=0.05.
		d) DIRECTION=$OPTARG;; 	# set to 1 to get only enriched sgRNA in fg. -1 for only depleted sgRNA. Default=1.
		m) METHOD=$OPTARG;; 	# method used to call aggregate pvalue from pooled sgRNA pvalues. Choose from [RRA,min,geom.mean,median,stuart]. Default=RRA
		a) RRACUTOFF=$OPTARG;;  # pval cutoff for crest-seq signal to be retained in region and peak. Default=1.
		c) MINSGRNA=$OPTARG;;	# each bin >= $MINSGRNA positive sgRNA to be included in peak calling. Default=2.
		l) MINLEN=$OPTARG;; 	# min length of final peaks. Default=1000.
		\?) echo "Invalid option: -$OPTARG" 1>&2; exit 1;;
	esac
done

## start run
echo "CRISPY started runing ..."
echo "./crispy $@"

## Check arguments
echo "######### 0. checking parameters #########"
## check on negative label:
if grep -q "negative" "$INREAD"; then
	echo "negative label successfully found in reads file: $INREAD"
else
	echo "In order for the sgRNA FDR calculation to work, you need to label 'negative' in your reads file for the control sgRNAs: $INREAD"
	exit
fi

## reads -> sgRNA
echo "######### 1. read counts -> sgRNA signals ... #########"
./bin/call.gRNA.R --inFile=$INREAD --bg="$BG" --fg="$FG" --outDir="$OUTDIR" --prefix="$PREFIX" --plotFormat="$PLOTFORMAT" --qnorm="$QNORM" --min_cpm="$MINCPM" --min_cpm_ratio="$MINCPMRATIO" --pvalCut="$NBCUTOFF" --direction="$DIRECTION"

## sgRNA -> sgRNA signals in loci
#echo "######### getting sgRNA signal from pvalues ... #########"
pvalTsv="$OUTDIR/$PREFIX.sgRNA.tsv"  # from last step
sgrnaSignal="$OUTDIR/$PREFIX.sgRNA.bedgraph" # output of this step, good for visualizing "good" sgRNA signals.
sgrnaSignalAll="$OUTDIR/$PREFIX.sgRNA_all.bedgraph" # output of this step, good for visualizing "all" sgRNA signals.

## all sgRNA
tail -n+2 $pvalTsv | awk -v kw=$SGRNAKW -v OFS="\t" '{if($NF==kw){print $1, -log($5)*(($2>0)-0.5)*2}}' > tmp.1.all
## good sgRNA
if [ "$DIRECTION" -eq -1 ]; then
	echo "DIRECTION==-1, using only depleted sgRNA ..."
	tail -n+2 $pvalTsv | awk -v kw=$SGRNAKW -v nbcutoff=$NBCUTOFF -v OFS="\t" '{if($NF==kw && $5<nbcutoff && $2<0){print $1, -log($5)}}' > tmp.1
else
	echo "DIRECTION==1, using only enriched sgRNA ..."
	tail -n+2 $pvalTsv | awk -v kw=$SGRNAKW -v nbcutoff=$NBCUTOFF -v OFS="\t" '{if($NF==kw && $5<nbcutoff && $2>0){print $1, -log($5)}}' > tmp.1
fi

awk '{print $4"\t"$1"_"$2"_"$3}' $INSGRNA > tmp.2
./bin/rp tmp.1 tmp.2 | tr "_" "\t" > $sgrnaSignal
./bin/rp tmp.1.all tmp.2 | tr "_" "\t" > $sgrnaSignalAll

echo

## sgRNA signals -> target region signals
echo "######### 2. sgRNA signals -> target region signals ... #########"
paste -d"\t" <(cut -f1-3 $sgrnaSignal) <(cut -f4 $sgrnaSignal| ./bin/val2rank.py -nr) > tmp.3
# put sgRNA ranks into region bins
intersectBed -a $INREGION -b tmp.3 -wa -wb | ./bin/collapseBed -c 7 -o "list" -q  > tmp.4
./bin/rra.R --inFile="tmp.4" --outFile="tmp.5" --method=$METHOD --minSgRNA=$MINSGRNA
# convert and filter pvalues of region bins to signals. 
regionSignal="$OUTDIR/$PREFIX.region.bedgraph"
tail -n+2 tmp.5 | tr "_" "\t" | awk -v RRACUTOFF=$RRACUTOFF -v MINSGRNA=$MINSGRNA -v OFS="\t" '{if($4<RRACUTOFF && $5>= MINSGRNA){$4=-log($4); print}}' | ./bin/mySortBed > $regionSignal
echo 

### macs2 for peak smoothing ### remove this step for now. Try something else for peak calling. 
#echo "######### 3. target region signals -> peak smoothing ... #########"
### fill 0s for macs2
#subtractBed -a $INREGION -b $regionSignal | awk '{print $0"\t"0}' | cat - $regionSignal | ./bin/mySortBed > tmp; mv tmp $regionSignal
### call peaks
#peakSignal="$OUTDIR/$PREFIX.peak.bedgraph"
#macs2 bdgpeakcall -i $regionSignal -l $MINLEN -g $MAXGAP -c $PEAKCUTOFF -o /dev/stdout | tail -n+2 | cut -f1-3,5 > $peakSignal

### use inhouse script for merging. curated from CREST-seq code. ###
echo "######### 3. target region signals -> peak smoothing ... #########"
peakSignal="$OUTDIR/$PREFIX.peak.bedgraph"
# merge to MINLEN and sum signals at each peak.
mergeBed -i $regionSignal -c 4 -o max | awk -v MINLEN=$MINLEN -v OFS="\t" '{if($3-$2<=MINLEN){center=int(($2+$3)/2+0.5);ext=int(MINLEN/2+0.5); $2=center-ext; $3=center+ext}; print}' | mergeBed -c 4 -o sum > $peakSignal

## clean up
rm tmp.{1,2,3,4,5}*
echo -e "Crispy Done!\n\n\n"

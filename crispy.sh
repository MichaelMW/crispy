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
NBCUTOFF=1
DIRECTION=1
METHOD="RRA"
RRACUTOFF=1
MINSGRNA=2
MINLEN=1000
FCLIMIT=-0.5

## Hardcoded paramters
SGRNAKW="test"  ## keywords used to grep from sgRNA for the designed ones. 

## Get arguments
while getopts ":i:r:s:o:p:b:f:q:t:u:v:e:n:d:m:a:c:l" opt;
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
		e) FCLIMIT=$OPTARG;; 	# fold change limit for CRISPY signal track. Default: -0.5
		n) NBCUTOFF=$OPTARG;;  	# pval cutoff for sgRNA using negative bionomial test. Default=0.05.
		d) DIRECTION=$OPTARG;; 	# set to 1 to get only enriched sgRNA in fg. -1 for only depleted sgRNA. Default=1.
		m) METHOD=$OPTARG;; 	# method used to call aggregate pvalue from pooled sgRNA pvalues. Choose from [RRA,min,geom.mean,median,stuart]. Default=RRA
		a) RRACUTOFF=$OPTARG;;  # pval cutoff for crest-seq signal to be retained in region and peak. Default=1.
		c) MINSGRNA=$OPTARG;;	# each bin >= $MINSGRNA positive sgRNA to be included in peak calling. Default=2.
		l) MINLEN=$OPTARG;;     # min length of final peaks. Default=1000.
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
	echo "SUCCESS: negative label successfully found in reads file: $INREAD"
else
	echo "FAILED: In order for the sgRNA FDR calculation to work, you need to label 'negative' in your reads file for the control sgRNAs: $INREAD"
	exit
fi
## curate input files for "^M" label. (Rare, but it happens ...)
tr -d '\15\32' <$INREAD > tmp; mv tmp $INREAD
tr -d '\15\32' <$INREGION > tmp; mv tmp $INREGION
tr -d '\15\32' <$INSGRNA > tmp; mv tmp $INSGRNA

## reads -> sgRNA
echo "######### 1. read counts -> sgRNA signals ... #########"
./bin/call.gRNA.R --inFile=$INREAD --bg="$BG" --fg="$FG" --outDir="$OUTDIR" --prefix="$PREFIX" --plotFormat="$PLOTFORMAT" --qnorm="$QNORM" --min_cpm="$MINCPM" --min_cpm_ratio="$MINCPMRATIO" --pvalCut="$NBCUTOFF" --direction="$DIRECTION"

## sgRNA -> sgRNA signals in loci
#echo "######### getting sgRNA signal from pvalues ... #########"
pvalTsv="$OUTDIR/$PREFIX.sgRNA.tsv"  # from last step
fcSignal="$OUTDIR/$PREFIX.fc.bedgraph" # output 1. Track for high resolution ean fold change in input regions.
sgrnaSignalAll="$OUTDIR/$PREFIX.sgRNA_all.bedgraph" # output 2. Track for pvalue of "all" sgRNA signals.
sgrnaSignal="$OUTDIR/$PREFIX.sgRNA.bedgraph" # output 3. Track for pvalue of "good" sgRNA signals.

## mean sgRNA fold change track.
# get FC
tail -n+2 $pvalTsv | awk -v kw=$SGRNAKW -v OFS="\t" -v fcLimit=$FCLIMIT '{if($NF==kw){FC=($2>fcLimit?$2:fcLimit); print $1, FC}}' > tmp.1.fc
# map to genome loci
awk '{print $4"\t"$1"_"$2"_"$3}' $INSGRNA > tmp.2.fc
./bin/rp tmp.1.fc tmp.2.fc | tr "_" "\t" > tmp.3.fc
# get the best possible resolution by moving between boundaries of $INSGRNA
tail -n+2 $INSGRNA | ./bin/overlap2cont.py > tmp.4.fc
# mean signal in region
intersectBed -a tmp.4.fc -b tmp.3.fc -wa -wb | ./bin/collapseBed -c 7 -o "mean" -q  > $fcSignal

## all sgRNA pvalue
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
intersectBed -a tmp.4.fc -b tmp.3 -wa -wb | ./bin/collapseBed -c 7 -o "list" -q  > tmp.4.hres
./bin/rra.R --inFile="tmp.4" --outFile="tmp.5" --method=$METHOD
./bin/rra.R --inFile="tmp.4.hres" --outFile="tmp.5.hres" --method=$METHOD
# convert and filter pvalues of region bins to signals. 
regionSignal="$OUTDIR/$PREFIX.region.bedgraph"
regionSignalHres="$OUTDIR/$PREFIX.region.hres.bedgraph"
tail -n+2 tmp.5 | tr "_" "\t" | awk -v RRACUTOFF=$RRACUTOFF -v MINSGRNA=$MINSGRNA -v OFS="\t" '{if($4<=RRACUTOFF && $5>= MINSGRNA){$4=-log($4); print}}' | ./bin/mySortBed > $regionSignal
tail -n+2 tmp.5.hres | tr "_" "\t" | awk -v RRACUTOFF=$RRACUTOFF -v MINSGRNA=$MINSGRNA -v OFS="\t" '{if($4<=RRACUTOFF && $5>= MINSGRNA){$4=-log($4); print}}' | ./bin/mySortBed > $regionSignalHres
echo 

#echo "######### 3. target region signals -> peak smoothing ... #########"
### use inhouse script for merging. curated from CREST-seq code. ###
echo "######### 3. target region signals -> peak smoothing ... #########"
peakSignal="$OUTDIR/$PREFIX.peak.bedgraph"
#peakSignalHres="$OUTDIR/$PREFIX.peak.hres.bedgraph"
# merge to MINLEN and sum signals at each peak.
mergeBed -i $regionSignal -c 4 -o max | awk -v MINLEN=$MINLEN -v OFS="\t" '{if($3-$2<=MINLEN){center=int(($2+$3)/2+0.5);ext=int(MINLEN/2+0.5); $2=center-ext; $3=center+ext}; print}' | mergeBed -c 4 -o sum > $peakSignal
#mergeBed -i $regionSignalHres -c 4 -o max | awk -v MINLEN=$MINLEN -v OFS="\t" '{if($3-$2<=MINLEN){center=int(($2+$3)/2+0.5);ext=int(MINLEN/2+0.5); $2=center-ext; $3=center+ext}; print}' | mergeBed -c 4 -o sum > $peakSignalHres

## clean up
rm tmp.{1,2,3,4,5}*
echo -e "Crispy Done!\n\n\n"

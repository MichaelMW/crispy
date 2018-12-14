# Crispy
A lightweight versatile pipeline for crispr-screening analysis.

## Demo

* run `./demo.sh` to check if the pipeline works
* run `./test.sh` to test the latest feature

## More Examples
* read `moreExample.sh` for optimal parameters for your experiments. 

## Read Counts
* to generate read count table as input for `crispy.sh`, checkout `counter/counter.demo.sh`
* run `counter/counter.py -h`

## Todos
* discrepency between regions and peaks. (optimize macs2 parameters.)
* automate read counts and QC

## history
* v1.4.1:
	* included an "sgRNA_all.bedgraph" track in the output
* v1.4:
 	* final peak calling using in-house script replicating CREST-seq
* v1.3.8:
 	* Included FDR guideline for "-n" pvalue cutoff choice. 
* v1.3.7:
	* Highlight positive sgRNA based on -n and -d		
* v1.3.6:
	* Now qc plot includes negative and positive sgRNAs as reference to guide user pick pvalue cutoff for '-n $NBCUTOFF'.

* v1.3.5:
	* Qnorm supports batch operator using ";" and ",". eg. -q "cis1,cis2;ctr1,ctr2;high1,high2"
	* More QC option
		* MIN_CPM filter (-u) and MIN_CPM_RATIO filter (-v). 
	* More QC plot
		* Violin + Boxplot for read count distribution
		* Read count density plot

* v1.3:
	* added support for direction (use -d 1 for enriched sgRNA [up-regulator], -d -1 for depleted sgRNA [down-regulator] )
	* added support for aggregation methods (use -m [RRA, min,geom.mean,median,stuart]. Default=RRA)
	* added support for quantile normalization within fgs and bgs. (use -q flag. Defualt=False if not specified.)
	* bug fix
	* removed cluttered output files. Each run now have sgRNA signals, region signals and final peak signals. 

* v1.2:
	* Counter template included. 
	* Support for experiments without replication.


## Licence
* MIT

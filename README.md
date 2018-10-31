# Crispy
A lightweight versatile pipeline for crispr-screening analysis.

## Demo

* run `./demo.sh` to check if the pipeline works

## More Examples
* read `moreExample.sh` for optimal parameters for your experiments. 

## Read Counts
* to generate read count table as input for `crispy.sh`, checkout `counter/counter.demo.sh`
* run `counter/counter.py -h`

## Todos
* discrepency between regions and peaks. (optimize macs2 parameters.)
* Check on depletion vs enrichment (done)
* Quantile normalization
* PCA
* automate read counts (low)


## history
* v1.3:
	* added support for direction (use -d 1 for enriched sgRNA [up-regulator], -d -1 for depleted sgRNA [down-regulator] )
	* added support for aggregation methods (use -m [RRA, min,geom.mean,median,stuart]. Default=RRA)
	* added support for quantile normalization within fgs and bgs. (use -q flag. Defualt=False if not specified.)
	* bug fix

* v1.2:
	* Counter template included. 
	* Support for experiments without replication.
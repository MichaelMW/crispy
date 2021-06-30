#!/bin/sh


## test input
inRegion="input/mmr_ATAC_EXON_regions_unique.bed"
inOligo="input/MMR_unique_oligos_03272019_1.tsv"
inSgRna="input/MMR_CRISPRi_XR058_59(unsorted)_XR066_67(6tg_80)_chr2.sgRNA.tsv"
pvalCut=0.039894

./getStats.sh $inRegion $inOligo $inSgRna $pvalCut > mergedStats.tsv



#!/bin/sh

chrm="chr3"
start=31000000
end=37000000
res=50

awk -v chrm=$chrm -v start=$start -v end=$end -v res=$res -v OFS="\t" 'BEGIN {for (i= start; i <= end - res; i = i + res) print chrm, i, i+res}'


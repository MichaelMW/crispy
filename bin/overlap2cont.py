#!/usr/bin/env python
# encoding: utf-8

## input a bedfile with overlapping regions
## output a bedfile with continous regions divided by incremental boundaries from input

from collections import defaultdict
from sys import stdin
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)


## check input is number
def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

## add boundaries to chrm dictionary
chrm2poss = defaultdict(list)
for line in stdin.readlines():
	ls = line.strip().split()
	chrm, start, end = ls[:3]
	if is_number(start) and is_number(end):
		chrm2poss[chrm].append(int(start))
		chrm2poss[chrm].append(int(end))

## convert boundaries to regions
for chrm, poss in chrm2poss.items():
	poss_sorted = sorted(poss)
	for i, pos in enumerate(poss_sorted):
		if i==0:
			continue
		prepos = poss_sorted[i-1]
		locus = "{}\t{}\t{}".format(chrm, prepos, pos)
		print(locus)



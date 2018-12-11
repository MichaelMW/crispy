#!/usr/bin/env python
# encoding: utf-8

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)


#### read oligo id
oligos = []
inFile = "oligos.exon.tsv"
with open(inFile) as f:
	for l in f.readlines():
		ls = l.strip().split()
		oligoID = ls[3]
		oligos.append(oligoID)

#### filter 
inFile = "reads.tsv"
with open(inFile) as f:
	header = f.readline().rstrip()
	print(header)
	for l in f.readlines():
		ls = l.strip().split()
		oligoID = ls[0]
		status = ls[-1]
		if oligoID in oligos or status != "test":
			print(l.rstrip())
		


#!/usr/bin/env python
# encoding: utf-8

####### todo:
# reverse compliment #


### input fasta file
### input gRNA.bed, design file
### output count number for next step. 

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
from sys import stdin, stderr
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f1', dest='flank1', help="left flank of sgRNAseq, before rc")
parser.add_argument('-f2', dest='flank2', help="right flank of sgRNAseq, before rc")
parser.add_argument('-R', dest='rc', default="0", help="reverse compliment mode, 0: only input; 1: rc; this is useful for sgRNA2 when the input fastq files are paired reads.")
parser.add_argument('-H', dest='hasHeader', default=True, help="oligo file has header")
parser.add_argument('-i', dest='oligo', help="input oligo design file. Use the following fields as default: chrm, start, end, sgRNAid, barcode, set, sgRNA, [sgRNA2 if pair reads]; if user provides a non-conventional oligo design file, minimally it should contain sgRNAid and sgRNAseq, and -G and -B should be used to indicate which columns (the column index starts with 0) sgRNAid and sgRNAseq are in the -i input file. eg. -i 'myOligo.tsv' -G 0 -S 1")
parser.add_argument('-I', dest='sgRNAid', default="3", help="column index of sgRNAid in the input oligo design file, index starts with 0. default to 3.")
parser.add_argument('-E', dest='exact', default=True, help="if pattern resulted in exact read match.")
parser.add_argument('-S', dest='sgRNAseq', default="6", help="column index of sgRNAseq in the input oligo design file, index starts with 0. default to 6 for single reads, specified to 7 for usage of sgRNA2.")

args = parser.parse_args()
flank1 = args.flank1.upper()
flank2 = args.flank2.upper()
inFile = args.oligo
rc = str(args.rc)
idx_sgRNAid = int(args.sgRNAid)
idx_sgRNAseq = int(args.sgRNAseq)
hasHeader = bool(args.hasHeader)
exact = bool(args.exact)

###
# reverse compliment
def func_rc(seq):
	string = seq.upper()
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	return "".join([complement.get(base, base) for base in reversed(seq)])


### 1. build gRNA profile
#inFile = "gRNA.bed"
gid2bar = {}
gids = []
with open(inFile) as f:
	if hasHeader:
		header = f.readline()
	else:
		pass
	for l in f.readlines():
		ls = l.strip().split()
		gid, sgRNAseq = ls[idx_sgRNAid], ls[idx_sgRNAseq]
		gids.append(gid)
		sgRNAseq = sgRNAseq.upper()
		if rc == "1":
			gid2bar[gid] = func_rc(sgRNAseq)
		else:
			gid2bar[gid] = sgRNAseq


sgRNAseqs = set(gid2bar.values())
blens = [len(sgRNAseq) for sgRNAseq in sgRNAseqs]
maxbl, minbl = max(blens), min(blens)

### 2. compile regex
import re

if rc == "0":
	pattern = flank1 + "([ACGT]{" + str(minbl) + "," + str(maxbl) + "})" + flank2 
if rc == "1":
	pattern = func_rc(flank2) + "([ACGT]{" + str(minbl) + "," + str(maxbl) + "})" + func_rc(flank1)
stderr.write("using pattern = " + pattern + "\n")
prog = re.compile(pattern)

### 3. search for sgRNA pattern that matches

hits = []
hitc = 0
readc = 0
for seq in stdin.readlines():
	readc += 1
	seq = seq.upper()
	hit = prog.search(seq)
	# this step that act like a fuzzy match improves runtime
	if hit:
		#stderr.write(hit + "\n")
		match = hit.group()
		matchSeq = match[len(flank1):len(match)-len(flank2)]
		hitc += 1
		# this step checks if the sequence is actually in the original gRNA design
		if exact: # very fast
			if matchSeq in sgRNAseqs:
				if rc == "1":
					matchSeq = func_rc(matchSeq)
				hits.append(matchSeq)
		else: # very slow
			for sgRNAseq in sgRNAseqs:
				if sgRNAseq in matchSeq:
					hits.append(sgRNAseq)
					continue
			#stderr.write(sgRNAseq + "\n")

### count
from collections import Counter
hitCounts = Counter(hits)

#print(hitCounts)

### print
if rc == "1":
	for gid in gids:
		sgRNAseq = gid2bar[gid]
		if func_rc(sgRNAseq) in hitCounts:
			print("\t".join([gid, str(hitCounts[func_rc(sgRNAseq)])]))
elif rc == "0":
	for gid in gids:
		sgRNAseq = gid2bar[gid]
		if sgRNAseq in hitCounts:
			print("\t".join([gid, str(hitCounts[sgRNAseq])]))
else:
	print("rc mode unknown, exit")
	exit


stderr.write("total read count: " + str(readc)+"\n")
stderr.write("total read match (pattern): " + str(hitc)+"\n")
stderr.write("total exact match (sgRNA): " + str(len(hits))+"\n")


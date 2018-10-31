#!/usr/bin/env python


##### convert values to ranks, tie breaker as average #####
from __future__ import division
from sys import argv, stdin, stdout
from signal import signal, SIGPIPE, SIG_DFL
import argparse
signal(SIGPIPE, SIG_DFL)

# parse args
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inVals", help="input values")
parser.add_argument("-d", "--delimiter", default="\n",
                    help="used to split and join input and output")
parser.add_argument("-n", "--normalize", action='store_true',
                    default=False, help="normalize over length of input to [0,1]")
parser.add_argument("-r", "--reverse", action='store_true',
                    default=False, help="high intial value leads to low ranks")
args = parser.parse_args()
sep = args.delimiter

### functions ###

def rank_simple(vector, reverse):
    return sorted(range(len(vector)), key=vector.__getitem__, reverse=reverse)

def rankdata(a, sep, normalize, reverse):
    n = len(a)
    ivec = rank_simple(a, reverse)
    svec = [a[rank] for rank in ivec]
    sumranks = 0
    dupcount = 0
    newarray = [0]*n
    for i in range(n):
        sumranks += i
        dupcount += 1
        if i == n-1 or svec[i] != svec[i+1]:
            averank = sumranks / float(dupcount) + 1
            for j in range(i-dupcount+1, i+1):
                if normalize:
                    newarray[ivec[j]] = str(averank/n)
                else:
                    newarray[ivec[j]] = str(averank)
            sumranks = 0
            dupcount = 0
    return sep.join(newarray)

inData = [float(i) for i in list(stdin.readlines())]
stdout.write(rankdata(inData, sep=sep, normalize=args.normalize, reverse=args.reverse))

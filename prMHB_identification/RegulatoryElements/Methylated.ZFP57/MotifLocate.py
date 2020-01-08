#!/usr/bin/env python

# Time-stamp: <2020-01-08 Zhaowei Yu>
"""Description: Locate the location of the motif in the genome
@author: Zhaowei Yu
@contact: edwardbioyu@gmail.com
"""

import os,sys
import re
import twobitreader

inf = open(sys.argv[1],"rU")
outf = open(sys.argv[2],"w")
motif = sys.argv[3]
genome = twobitreader.TwoBitFile('/mnt/Storage/home/yuzhaowei/data/Genome/mm9/mm9.2bit')

inf_lines = inf.readlines()
for line in inf_lines:
	chrom = line.strip().split('\t')[0]
	start = line.strip().split('\t')[1]
	end = line.strip().split('\t')[2]
	peak = line.strip().split('\t')[3]

	sequence = genome[chrom][int(start):int(end)].upper()
	pattern = re.compile(motif)
	if motif in sequence:
		match = pattern.finditer(sequence)
		for s in match:
			loci = [int(start) + int(s.start()),int(start) + int(s.start()) + len(motif)]
			outf.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n" . format(chrom, loci[0], loci[1], peak, start, end))

inf.close()
outf.close()

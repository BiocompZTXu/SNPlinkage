#!/usr/bin/env python
# coding=utf-8
# Author: Xu Zhongtian
# Mail: 738502908@qq.com
# Created Time: Sun 01 May 2016 12:34:02 AM CST
# Filename: calculating_base_frequency_each_pos.py
#BAM file is 0-based coordinate
#SAM fole, mpileup file is 1-base coordinate

import sys
import warnings

if len(sys.argv)!=2:
	print (" This scripts process mpileup format file to calculate nucleotide freqency at each postion \
			\n So please prepare mipleup file with the following command: \
			\n $ samtools mpileup -f /path/genome.fa /path/sortBamFile.bam \
			\n Usage: python %s mpileup.txt "%sys.argv[0])
	raise SystemExit(1)
inFile = open(sys.argv[1],'r')

print 'Chrom\tPos\tREF\tDepth\tA\tG\tC\tT\tdeleletion\tinsertion\tprimary_allele\tprimary_allele_freq\tsecondary_allele\tsecondary_allele_freq'
for line in inFile:
	line = line.strip().split('\t')
	if len(line) >= 5:
		chrom = line[0]; pos = line[1]; ref = line[2]; depth = line[3]
		types = {'A':0,'G':0,'C':0,'T':0,'deletion':0,'insertion':0}
		
		ATCG = ["A","T","C","G"]
		ATCG.remove(ref)

		bases = line[4].upper()

		types[ref] =  (bases.count(".") + bases.count(","))# / float(depth)
		for nucleotide in ATCG: 
			types[nucleotide] = bases.count(nucleotide)# / float(depth)
		types['deletion'] = bases.count("-") #/ float(depth)
		types['insertion'] = bases.count("+") #/ float(depth)
		count_sorted = types.keys()
		count_sorted.sort(key=types.__getitem__, reverse=True)
		primary_allele = count_sorted[0]
		primary_allele_freq = types[count_sorted[0]] / float(depth)
		secondary_allele = count_sorted[1]
		secondary_allele_freq = types[count_sorted[1]] / float(depth)


#		if (int(depth) != sum(types.values())):
#			warnings.warn("Total count do not add up to coverage. sum of counts: %i coverage: %i\n" % (sum(types.values()), depth))
		
	
	print "\t".join([chrom,pos,ref,depth]) + "\t"+ "\t".join(map(str,[types['A'],types['G'],types['C'],types['T'],types['deletion'],types['insertion'],primary_allele,round(primary_allele_freq,2),secondary_allele,round(secondary_allele_freq,2)]))
		
	
	




#!/usr/bin/env python
# coding=utf-8
# Author: Xu Zhongtian
# Mail: 738502908@qq.com
# Created Time: Fri 29 Apr 2016 12:30:29 PM CST
# Filename: calculating_breadth_depth_according_BAM.py
# Last modifed time: Mon Apr  3 11:05:45 2017
# And BAM.file used be indexed by: samtools index map_sorted.bam

import sys
import os
import re 
import pprint
from subprocess import Popen, PIPE

def get_contigs(bam):
    header, err = Popen(["samtools","view","-H",bam], stdout=PIPE, stderr=PIPE).communicate()
    #if err != "":
    #    raise Exception(err)
    contigs = {}
    for x in re.findall("@SQ\WSN:(?P<chrom>[A-Za-z0-9_\.]*)\WLN:(?P<length>[0-9]+)", header.decode("utf-8")):
        contigs[x[0]] = int(x[1])
    return contigs

def coverage(bam):
    if os.path.isfile(bam) == False:
        raise Exception("Bam file does not exist")
    contigs = get_contigs(bam)
	
	
    coverage_dict = dict()
    for c in contigs.keys():
        command = "samtools depth -r %s %s | awk '{sum+=$3;cnt++}END{print cnt \"\t\" sum}'" % (c, bam)
        coverage_dict[c] = {}
        if  Popen(command, stdout=PIPE, shell = True).communicate()[0].decode("utf-8").strip().split("\t") == ['']:
            coverage_dict[c]["Bases_Mapped"], coverage_dict[c]["Sum_of_Depths"] = 0,0
            coverage_dict[c]["Breadth_of_Coverage"] = coverage_dict[c]["Bases_Mapped"] / float(contigs[c])
            coverage_dict[c]["Depth_of_Coverage"] = coverage_dict[c]["Sum_of_Depths"] / float(contigs[c])
            coverage_dict[c]["Length"] = int(contigs[c])
        else:
            coverage_dict[c]["Bases_Mapped"], coverage_dict[c]["Sum_of_Depths"] = map(int,Popen(command, stdout=PIPE, shell = True).communicate()[0].decode("utf-8").strip().split("\t"))
            coverage_dict[c]["Breadth_of_Coverage"] = coverage_dict[c]["Bases_Mapped"] / float(contigs[c])
            coverage_dict[c]["Depth_of_Coverage"] = coverage_dict[c]["Sum_of_Depths"] / float(contigs[c])
            coverage_dict[c]["Length"] = int(contigs[c])


#Calculate Genome Wide Breadth of Coverage and Depth of Coverage
    genome_length = float(sum(contigs.values()))
    print(genome_length)
    coverage_dict_genome = dict()
    coverage_dict_genome["genome"] = {}
    coverage_dict_genome["genome"]["Length"] = int(genome_length)
    coverage_dict_genome["genome"]["Bases_Mapped"] = sum([x["Bases_Mapped"] for k, x in coverage_dict.items() if k != "genome"])
    coverage_dict_genome["genome"]["Sum_of_Depths"] = sum([x["Sum_of_Depths"] for k, x in coverage_dict.items() if k != "genome"])
    coverage_dict_genome["genome"]["Breadth_of_Coverage"] = sum([x["Bases_Mapped"] for k, x in coverage_dict.items() if k != "genome"]) / float(genome_length)
    coverage_dict_genome["genome"]["Depth_of_Coverage"] = sum([x["Sum_of_Depths"] for k, x in coverage_dict.items() if k != "genome"]) / float(genome_length)

    return coverage_dict,coverage_dict_genome


#Main part 
if __name__ == "__main__":
    for i in range(1,len(sys.argv)):
        print("\t" + "\t".join(['length','covered_base_count','Breadth_of_coverage','total_base_mapped','average_Depth']))
        bam = sys.argv[i]
        coverage_dict, coverage_dict_genome = coverage(bam)
        for chrom,stats in sorted(coverage_dict.items(), key=lambda x: int(re.findall('\d+',x[0])[0])):
            print(chrom+"\t"+'\t'.join(map(str,[stats['Length'],stats["Bases_Mapped"],round(stats['Breadth_of_Coverage'],4),stats["Sum_of_Depths"],round(stats["Depth_of_Coverage"],2)])))
        for chrom,stats in sorted(coverage_dict_genome.items()):
            print(chrom+"\t"+'\t'.join(map(str,[stats['Length'],stats["Bases_Mapped"],round(stats['Breadth_of_Coverage'],4),stats["Sum_of_Depths"],round(stats["Depth_of_Coverage"],2)])))

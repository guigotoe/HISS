#!/usr/bin/python
# -*- coding: utf-8 -*-

# Guillermo Torres MSc,
# ggtorrese@unal.edu.co, guigotoe@gmail.com
# August 2012
# This script make random data set.
# ********************* How to use ***************************
# 
# 
# ************************************************************
###################
     
import random, sys, re
import numpy as np 
from Bio import SeqIO

def main():
    infile =  sys.argv[1] #"rsstest.txt"#
    outfile = open("ave_%s"%infile,'w')
    probenumber = 0
    genenumber = 4 # number of random genes; default 20 random genes
    Probes = {}
    GENES = {}
    
    for line in open("%s"%infile,'r'):
        line = line.strip("\n").split()
        ids = re.split("(_\d+)$",line[0])
        name = ids[0]
        geneid = int(ids[1].strip("_"))
        ssvalue = float(line[1])       
        try:
            if GENES[name]:
                pass
        except KeyError:
            pscore = aveprobe(GENES)
            for k in GENES.keys():
                write(pscore,k,outfile)
            Probes = {}
            GENES = {}       
        try:
            Probes[geneid].append(ssvalue)
        except KeyError:
            Probes[geneid] = []
            Probes[geneid].append(ssvalue)  
        GENES[name] = Probes   
    pscore = aveprobe(GENES)
    for k in GENES.keys():
                write(pscore,k,outfile)
        
def write(pscore,name,outf):
    for k in range(len(pscore)):
        outf.write("%s\t%s\n"%(name, round(float(pscore[k]),4)))
    
def aveprobe(gene):
    if gene == {}:
        return
    else:
        for k in gene.keys():
            matrix = np.matrix(gene[k].values())
            avePscore = matrix.mean(axis=0).tolist()[0]
            return avePscore

if __name__ == "__main__": main()
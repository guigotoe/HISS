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
    infile = sys.argv[1] # "rNmax.txt"#
    outfile = open("ave_%s"%infile,'w')
    SEED={}
    for line in open("%s"%infile,'r'):
        line=line.strip("\n").split()
        name=re.split("(_\d+)$",line[0])[0]
        try :
            SEED[name][0].append(line[1])
            SEED[name][1].append(line[2])
        except KeyError:
            SEED[name]=[[],[]]
            SEED[name][0].append(line[1])
            SEED[name][1].append(line[2])
        
    for k in SEED.keys():
        outfile.write("%s\t%s\t%s\n"%(k,round(np.array(SEED[k][0],float).mean(),2),round(np.array(SEED[k][1],float).mean(),2)))


if __name__ == "__main__": main()
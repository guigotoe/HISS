#!/usr/bin/python

# Guillermo Torres MSc,
# ggtorrese@unal.edu.co, guigotoe@gmail.com
# August 2012
# This script make random data set.
# ********************* How to use ***************************
# make_random_dataset Number_Of_Sequences Sequence_Length 
# 
# ************************************************************
###################
     
import random, sys
from Bio import SeqIO

def main():
    rdataset = open("randomDS.fa","w")
    for record in SeqIO.parse("%s"%(sys.argv[1]), "fasta"):
        i =1
        while i <=int(sys.argv[2]):
            randomseq = "".join(random.sample(record.seq,len(record.seq)))
            rdataset.write(">%s_%s\n%s\n"%(record.id,i,randomseq))
            i+=1
            
  

if __name__ == "__main__": main()
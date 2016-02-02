#!/usr/bin/python
# Guillermo Torres MSc,
# ggtorrese@unal.edu.co, guigotoe@gmail.com
# Bioinformatic Research Group
# Biotechnology Institute of National University of Colombia
# 17 August 2011
# update November 2013
#* last update September 2014
#* trimming Ambioguosity "N", only the unambiougosity sequence that represent bigged than 70% of
#* whole sequence is retained, else is not considered for future analysis.
# This is written as part of IN SILCO MICROARRAY pipeline, but could be split to serve to different purposes.
#** How to use ****
#** fieDeb.py clustersimilarity(0-1) infilename
#** example: python fileDeb.py --infa path/nirs.fasta --clusID 0.85 --out genes.fa 
#** 
#** END ***

import os, time, sys, shutil, re
from tempfile import mkstemp
import subprocess as sub
from optparse import OptionParser
from Bio import SeqIO

## PARSING... ##
 #* Options *#
parser = OptionParser()
usage = """\nfileDeb exploits the sequence similarities, clustering them, in order to estimate a biological sequence from consensus sequence derived from a multiple alignment using UCLUST\n 
REQUIREMENTS:\n** Usearch v 5.2\n
\t%prog --infa IN_fileNAME --culsID 0.85 --out outFileName\n\t%prog {-h --help}\n
Part of in silico Hybridization system by Guillermo G. Torres (ggtorrese@unal.edu.co) - Bioinformatics
Group - IBUN - National University of Colombia"""
parser = OptionParser(usage=usage,version="%prog 1.0")
parser.add_option("--infa","--input_fasta",type="string",action="store",dest="multifa",#default='nirs.fasta',#metavar=" infile_NAME",#default='nirs.fasta', 
                      help="Input multifasta file with raw genes information (from Data Bases)")
parser.add_option("--out","--out_file",type="string",action="store",dest="out_file_name",default='genes.fa',#metavar=" out_file_NAME" 
                      help="Out multifasta file name of consensus gene sequences - default name = genes.fa")
parser.add_option("--clusID","--cluster_id",type="float",action="store",dest="clus_id",default=0.85,#metavar=" Genes clustering identity % (float)",
                      help="Denoising identity percentage - default % value = 0.85")
parser.add_option("-S","--Seed", dest="seed",action="store_true",default=False,
                      help="Extract the seeds sequences from cluster")
(o,args) = parser.parse_args()
if len(args)!=0:
    parser = OptionParser(usage="\n\t%prog --infa IN_fileNAME --culsID 0.85 --out outFileName\n\t%prog {-h --help}")
    parser.error("incorrect number of arguments\n")
path = os.getcwd()   
def main():
    '''
    Check the unambiguosity of the sequences and then send them to cluster process
    '''
    inf = o.multifa
    clID = o.clus_id
    out = o.out_file_name
    
    ## Removing ambiguosities...
    fd,entry = mkstemp(suffix='.tmp',prefix='_',dir=path)    # making temporary file tuple
    for seq_record in SeqIO.parse(inf, "fasta"):
        fragments = {}
        lenghts = []
        for fragment in seq_record.seq.split('N'):
            if len(fragment) > 0:
                lenghts.append(len(fragment))
                fragments[len(fragment)] = fragment
        if len(lenghts) > 0:
            longest = max(lenghts)
            fragment_seq = fragments[longest]
            rate = len(fragment_seq)/float(len(seq_record))
            if rate >= 0.7:
                os.write(fd,">%s\n%s\n"%(seq_record.description,fragment_seq))     # filling the temporary file with gene sequence in fasta format

    ## Clustering
    clustering(entry,clID,o.seed)
    os.close(fd)
    os.remove(entry)
    print ("Clustering Done!")
    
      
def clustering(inf,clID,out):
    '''
    Use Usearch to clusterize the input sequences
    '''
    outf = o.out_file_name
    sorting = "usearch --sort %s --output %s.sorted" % (inf, outf)
    consenso = "usearch --cluster %s.sorted --id %s --uc %s.uc --consout %s.fa --log %s.log" % (outf, clID, outf, outf, outf)
    seeds = "usearch --cluster %s.sorted --id %s --uc %s.uc --seedsout %s.fa --log %s.log" % (outf,clID,outf,outf,outf)
    
    sort = sub.Popen(sorting, shell = True, close_fds = True, stdout=sub.PIPE)
    time.sleep(1)
    if sort.wait() == 0:
        if out == 1:
            cluster = sub.Popen(seeds, shell = True, close_fds = True, stdout=sub.PIPE)
        else:
            cluster = sub.Popen(consenso, shell = True, close_fds = True, stdout=sub.PIPE)
        if cluster.wait() == 0:
            os.system("rm -r %s.sorted"%(outf))

if __name__ == "__main__": main()

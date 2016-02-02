#!/usr/bin/python
# Guillermo Torres MSc,
# ggtorrese@unal.edu.co, guigotoe@gmail.com
# Bioinformatic Research Group
# Biotechnology Institute of National University of Colombia
# last update: November 2013
# **HISS integrates: Comparison + Section processes
# This is written as part of IN SILCO MICROARRAY pipeline, but could be
# splitted to serve different porpuoses.
#/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Requirements:-Biopython >= v1.52; -New Blast installed (BLAST+)
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

import re, sys, os, math 
from tempfile import mkstemp
from collections import Counter
import subprocess as sub
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIStandalone
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastnCommandline


## PARSING... ##
 #* Options *#
parser = OptionParser()
usage = """\nThis scripts is able to make a in silico hibridization between probes from CrSo and reads from metagenome sequencing. 
Here 2 basic steps coming to scene; 1) Comparison: Execute BLAST to make alignments between probes and reads. 
2) Selection: retrieve the homology reads trought evaluation the lenght and mismatch of the alignments. 
The Results are presented in tabular file.\n REQUIREMENTS:\n** Biopython >= v 1.52\n** BLAST \n
\t%prog --prb PROBES_fileNAME --refDB path/DB_NAME [options] arg \n\t%prog {-h --help}\n
Part of in silico Hybridization system by Guillermo G. Torres (ggtorrese@unal.edu.co) - Bioinformatics
Group - IBUN - National University of Colombia"""
parser = OptionParser(usage=usage,version="%prog 1.0")
parser.add_option("-p","--probe_fasta",type="string",action="store",dest="multifa",metavar=" file_NAME",#default='./Proof/probes0.fa', 
                      help="Input multifasta file with probes of genes to map")
parser.add_option("-o","--out_file",type="string",action="store",dest="out_file_name",metavar=" out_file_NAME",default='hissOut', 
                      help="Out tabular files name - default name = hissOut")
parser.add_option("-r","--Illumina_reads",type="string",action="store",dest="reads_db",metavar=" path/DB_NAME",#default='./Proof/nirsdb',#L1_1XDgR',#nirsdb', 
                      help="Reference data base to evaluate the specificity of gene regions from which will be designed the probes")
parser.add_option("-n","--processors_num",type="int",action="store",dest="proc",metavar=" Number of processor used by BLAST (int)",default=1,
                      help="Number of processor used by BLAST - default value = 1")
parser.add_option("-v", action="store_true",dest="verbose",default=False,
                      help = "Verbose - default False")
(o,args) = parser.parse_args()
if len(args)!=0 or (o.multifa == None and o.reads_db == None):
    parser = OptionParser(usage="\n\t%prog -p path/IN_ProbefileNAME -r path/DB_NAME [options] arg \n\t%prog {-h --help}")
    parser.error("incorrect number of arguments\n")
   
## temp raw probes file ##
path = os.getcwd()


def main():
    outf = open('%s_counts.txt'%o.out_file_name,'w') 
    outf.write ("ProbeID\tNoGenesFounded\n")
    outfull = open('%s_full.txt'%o.out_file_name,'w')
    outfull.write("ProbeID\tSeqHibritedID\tAl_Len\tPerc_Identity\tMismatch\tAl_Probe_Pos\tAl_SeqHib_Pos\n") # Writing the header of the full file
  
    
    handle = open(o.multifa, "rU")
    fl = int(open(o.multifa,"rU").read().count(">")) # length of file
    lc=0 # line counter
    if o.verbose==True:sys.stderr.write("""
Your parameters:\ninput_fasta = %s\nreads_DB = %s\n
Your out file name = %s\nHibridization...\n"""%(o.multifa,o.reads_db,o.out_file_name))
    try:
        for record in SeqIO.parse(handle, "fasta"):
            if re.search("[NRYSWKMBDHV]",str(record.seq)) != None:pass
            else:
                fd,entry = mkstemp(suffix='.tmp',prefix='_',dir=path)    # making temporary file tuple
                os.write(fd,">%s\n%s"%(record.id,record.seq)) # filling the temporary file with probe sequence in fasta format)
                try: HISS(entry)     # HISS function execution
                except:sys.stderr.write('can not hibridate the sequence with ID: %s\n'%record.id)#;continue ## If BLAST cant run
                os.close(fd)
                os.remove(entry)
                
            ## Progress counter for verbose mode ##
            lc+=1
            if o.verbose==True:
                sys.stderr.write('\r'+'' *0)
                sys.stderr.write(str(int(lc*100/fl))+'%')
                sys.stdout.flush()
        sys.stderr.write('Hibridization was done!')
    except TypeError:sys.stderr.write('Error!... in probe file to Hibridization') 
    
def HISS(infa): 
    '''Input: Name file of the probe sequence in fasta format
    The function Execute BLAST (comparation process), then
    read the output file and execute the Selector function 
    which performs the selection process'''
         
    ## BLAST ## 

    blastn_cline = 'blastall -p blastn -i %s -d %s -m 0 -o blast.out -v 500 -F F -a %s'%(infa,o.reads_db,o.proc)
    blast = sub.Popen(str(blastn_cline),shell=True,close_fds=True,stdout=sub.PIPE,stderr=sub.PIPE )     # BAST execution
    proc_out, proc_err = blast.communicate()
    if proc_err != '':               ## Is it able to execute this version?
        blastn_cline = NcbiblastnCommandline(query=infa, db=o.reads_db,evalue=10,out='blast.out',dust="no",word_size=11,outfmt=0)
        blast = sub.Popen(str(blastn_cline),shell=True,close_fds=True,stdout=sub.PIPE,stderr=sub.PIPE )     # BAST execution
        proc_out, proc_err = blast.communicate()
          
    ## Geting Hits ##
    
    if blast.wait() == 0:
        Hits = getHits('blast.out') ## Run Selection process and function to write the results
   
def getHits(infile):
    ''' BLAST parser using Biopython
    Input: name of blast out file in standard ouput format
    Outputs: 2 files 
    ''' 
    inf = open(infile)
    parser = NCBIStandalone.BlastParser()
    error_parser = NCBIStandalone.BlastErrorParser(inf)
    iterator = NCBIStandalone.Iterator(inf, parser)
    err_iterator = NCBIStandalone.Iterator(inf, error_parser)
    
    ## *** Parsing *** ##
    QUERY,HITS = {},{}
    q,hits=0,0
    for record in iterator:
        hits = 0
        if record.alignments is []:
            hits = 0
            HITS[record.query] = hits
        else:
            q+=1
            QUERY[q] = [record.query,record.query_letters,record.database]
            for alignment in record.alignments:
                for hsp in alignment.hsps:                                  
                    if Selector(hsp,record.query_letters):  #-->## ** Selection Process **## 
                        hits+=1
                        FullFile(record,alignment,hsp)
                    else:pass
        HITS[record.query] = hits
    
    outf = open('%s_counts.txt'%o.out_file_name,'a') 
    for k in HITS.keys():
        outf.write ("%s\t%s\n"%(k,HITS[k]))         
 
def Selector(hsp,pLen):
    '''Selection Process:
    Input: HSP and probe length
    Output: Boolean -> True if the hsp is a hit or
    False if the hsp is not a hit
    '''
    ## Variables -> %ID, strand, Alignment start and end points for probe and read.
    identity = round(100*hsp.identities[0]/hsp.identities[1])
    mismatch = hsp.identities[1] - hsp.identities[0]
    
    #*** minimum length of the alignment according to probe lenght (pLen) 
       
    if int(pLen) > 25: minL = 35 
    else:minL = 21
  
    ## ** Selection Process throught selection function **##

    #err_rule = math.fabs(math.ceil(int(9.2729*(math.log(float(hsp.identities[0]))) - 29.428))) 
    #err_rule = math.fabs(math.ceil(int(9.6413*(math.log(float(hsp.identities[0]))) - 30.483))) 
    
    Poly = math.fabs(math.floor(0.0006*math.pow(float(hsp.identities[0]),2) + 0.1631*float(hsp.identities[0]) -4.0082))
    Log = math.fabs(math.ceil(int(9.9581*(math.log(float(hsp.identities[0]))) - 32.067))) 
    
    PolyRule = round((100 - (Poly*100/float(hsp.identities[0]))),2)
    LogRule = round((100 - (Log*100/float(hsp.identities[0]))),2)
    
    if int(hsp.identities[1]) >= minL and float(identity) >= LogRule: return True
        #print 'pLen:',pLen,'minL:',minL,'*Alg len:',int(hsp.identities[1]),'Log:',Log,'MM:',mismatch,'LogRule:',LogRule,'#%id:',identity
    else: return False 

def FullFile(record,alignment,hsp):
    '''This Funtion writes the output file
    Input: record of the query and hsp, info resulted by Comparison process
    Output: 1 File of hits with probes and the hibrited sequences with all 
    the alignment info
    '''
    identity = round(100*hsp.identities[0]/hsp.identities[1])
    mismatch = hsp.identities[1] - hsp.identities[0]
    if hsp.strand[0] == 'Plus': 
        query_end = int(hsp.query_start) + hsp.identities[1]
    else: 
        query_end = int(hsp.query_start) - hsp.identities[1]
    if hsp.strand[1] == 'Plus': 
        subjct_end = int(hsp.sbjct_start) + hsp.identities[1]
    else: 
        subjct_end = int(hsp.sbjct_start) - hsp.identities[1]
   
    PPos = "%s...%s"%(hsp.query_start,query_end)
    SPos = "%s...%s"%(hsp.sbjct_start,subjct_end)
    outfull = open('%s_full.txt'%o.out_file_name,'a')
    outfull.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(record.query,alignment.title.lstrip("> "),hsp.identities[1],identity,mismatch,PPos,SPos))
    #print record.query,alignment.title.rstrip("> "),hsp.identities[1],identity,mismatch,PPos,SPos#,'qend:',query_end,'send:',subjct_end

   
    
if __name__ == "__main__": main()
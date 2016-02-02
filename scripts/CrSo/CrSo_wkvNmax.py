#!/usr/bin/python
# Guillermo Torres MSc,
# ggtorrese@unal.edu.co, guigotoe@gmail.com
# Bioinformatic Research Group
# Biotechnology Institute of National University of Colombia
# last update: 24 - Dic 2014
#   * Working file; retrieve info about PIS and SS to analyse their behavior 
#   * Nmax of PIS corrected
# update: July 201
#   * Raw_file writing was corrected
#   * BLAST execution exception was corrected
#   * Function getScores were divided by 2 : getNucScore and getOligoScore
#   * In Function getProbes an ambiguosity of the oligo sequence check operation was implemented
#   * New recursive function for overlaping check operation were impelented
#   * More comments were placed
# * update: November 2013
#
# **CrSo integrates: 3 esential process: getScore + getProbes + probeDenoising
# This is written as part of IN SILCO MICROARRAY pipeline, but could be
# splitted to serve different porpuoses.
#/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Requirements:-Biopython >= v1.52; -Blast installed (BLAST+), Usearch V 5.2
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

import sys,re,os,time,itertools,errno
import numpy as np
from tempfile import mkstemp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIStandalone
import subprocess as sub
from collections import defaultdict
from collections import OrderedDict
from optparse import OptionParser

## PARSING... ##
 #* Options *#
parser = OptionParser()
usage = """\nThis scripts is able to design probes from multifasta gene file (denoised with 
filedebugger.py). Here 3 basic steps coming to scene; 1) getScore: calculate a position homology score 
for each gene. 2) getProbes: retrieve the best probes from each gene using a global score 
calculated with position homology score. 3) ProbeDenoising: verify the unique attribute of the probes. 
The probes are presented in fasta file.\n REQUIREMENTS:\n** Biopython >= v 1.52\n** Usearch v 5.2\n
\t%prog --infa IN_fileNAME --refDB path/DB_NAME [options] arg \n\t%prog {-h --help}\n
Part of in silico Hybridization system by Guillermo G. Torres (ggtorrese@unal.edu.co) - Bioinformatics
Group - IBUN - National University of Colombia"""
parser = OptionParser(usage=usage,version="%prog 1.0")
parser.add_option("-i","--infa","--input_fasta",type="string",action="store",dest="multifa",metavar=" file_NAME",default='nirs.fasta', 
                      help="Input multifasta file with genes for which will design probes")
parser.add_option("-o","--out","--out_file",type="string",action="store",dest="out_file_name",metavar=" out_file_NAME",default='probes.fa', 
                      help="Out multifasta file name with the probes of target genes - default name = probes.fa")
parser.add_option("-r","--refDB","--reference_DB",type="string",action="store",dest="ref_db",metavar=" path/DB_NAME",default='nirsdb', 
                      help="Reference data base to evaluate the specificity of gene regions from which will be designed the probes")
parser.add_option("-l","--pl","--probe_length",type="int",action="store",dest="probe_len",metavar=" probe length (int)",default=100,
                      help="Number of nucleotides that compound the probe - default value = 100nt")
parser.add_option("-p","--pn","--probe_number",type="int",action="store",dest="probe_num",metavar=" number of probe (int)",default=1,
                      help="Number o prober to map one gene - default value = 1")
parser.add_option("-s","--ov","--overlap",type="int",action="store",dest="p_overlap",metavar=" probe overlapping length (int)",default=0,
                     help="Length of probe overlapping, in nucleotides - default value = 0")
parser.add_option("-m","--minL","--minimum_length",type="int",action="store",dest="min_len",metavar=" Minimum alignment length (int)",default=7,
                      help="Alignment length Lower threshold for specificity evaluation - default value = 21nt")
parser.add_option("-x","--maxL","--maximum_length",type="float",action="store",dest="max_len",metavar=" Minimum alignment length (float)",default=0.85,
                      help="Alignment length highest threshold for specificity evaluation, measured by a percentage of gene aligned - default % value = 0.85")  
parser.add_option("-y","--maxS","--maximum_similarity",type="float",action="store",dest="max_sim",metavar=" Minimum alignment identity % (float)",default=0.85, 
                      help="Alignment Identity percentage highest threshold for specificity evaluation - default % value = 0.85")
parser.add_option("-c","--cov","--coverage",type="float",action="store",dest="tot_cov",metavar=" % of target gene coverage (float)",default=0.99,
                      help="Alignment length of homologous gene with target gene, measured by a percentage of gene aligned - default % value = 0.99")
parser.add_option("-t","--clusID","--cluster_id",type="float",action="store",dest="clus_id",metavar=" Probes clustering identity % (float)",default=0.85,
                      help="Denoising identity percentage - default % value = 0.85")
parser.add_option("-n","--proc","--processors_num",type="int",action="store",dest="proc",metavar=" Number of processor used by BLAST (int)",default=1,
                      help="Number of processor used by BLAST - default value = 1")
parser.add_option("-v", action="store_true",dest="verbose",default=False,
                      help = "Verbose - default False")
(o,args) = parser.parse_args()
if len(args)!=0 or (o.multifa == None and o.ref_db == None):
    parser = OptionParser(usage="\n\t%prog --infa path/IN_fileNAME --refDB path/DB_NAME [options] arg \n\t%prog {-h --help}")
    parser.error("incorrect number of arguments\n")
   
## for temp files ##
path = os.getcwd()

def main():
    rp = "raw_probes.fa"
    open(rp,'w')            ## Opening new raw_probe file
    handle = open(o.multifa, "rU")
    fl = int(open(o.multifa,"rU").read().count(">")) # length of file
    lc=0 # line counter
    if o.verbose==True:sys.stderr.write("""
Your parameters:\ninput_fasta = %s\nreference_DB = %s\n
Your out file name = %s\nGetting Probes...\n"""%(o.multifa,o.ref_db,o.out_file_name))
    for record in SeqIO.parse(handle, "fasta"):
        if len(str(record.seq)) < 100: 
            if o.verbose==True:sys.stderr.write("%s: Gene Sequence too short for probe design; min lenght:100nt"%record.id)
            pass
        else:
            fd,entry = mkstemp(suffix='.tmp',prefix='_',dir=path)    # making temporary file tuple
            os.write(fd,">%s\n%s"%(record.id,record.seq))     # filling the temporary file with gene sequence in fasta format
            try:CrSo(entry,record.id,record.seq,rp)     # Execution function CrSo
            except:sys.stderr.write('can not read sequence with ID %s'%record.id)
            os.close(fd)
            os.remove(entry)
            
    ## Progres counter for verbose mode ##
        lc+=1
        if o.verbose==True:
            sys.stderr.write('\r'+'' *0)
            sys.stderr.write(str(int(lc*100/fl))+'%')
            sys.stdout.flush()

    statinfo = os.stat(rp)
    if statinfo.st_size!=0:     ## checking if probes were designed, in order to execute the denoising process
        probeDenoising(rp)
    else:sys.stderr.write('There were no probes designed!, check your Input file')

def CrSo(infa,ID,seq,rp):
    '''Generates probes using as 
    input: gene sequence in fasta file, Identifier, genes equence and raw output file
    '''
    Probes = {} 
    ## BLAST ##
    
    blastn_cline = 'blastall -p blastn -i %s -d %s -m 0 -o blast.out -v 500 -F F -a %s'%(infa,o.ref_db,o.proc)
    blast = sub.Popen(str(blastn_cline),shell=True,close_fds=True,stdout=sub.PIPE,stderr=sub.PIPE )
    proc_out, proc_err = blast.communicate()
    if proc_err != '':
        blastn_cline = NcbiblastnCommandline(query=infa, db=o.ref_db,evalue=10,out='blast.out',dust="no",word_size=11,outfmt=0)
        blast = sub.Popen(str(blastn_cline),shell=True,close_fds=True,stdout=sub.PIPE,stderr=sub.PIPE )
        proc_out, proc_err = blast.communicate()
        
    ## Geting Score and Probes ##
    if blast.wait() == 0:
        nmax_file = open("Nmax.txt",'a')
        nscore = getNucScore('blast.out')  ## This function generate the score for each nucleotide of the gene 
        nmax_file.write("%s\t%s\n"%(ID,nscore))
        #oscore = getOligoScore(nscore)  ## This function uses the nucleotide score matrix to calculate the oligo specificity scores
        #Probes = getProbes(oscore,seq,o.probe_num,Probes)    ## This function select from the posible oligos the definitive probes.
        
        rp = open(rp,'a')
        for k in Probes.keys():
            rp.write(">%s_%s %s\n%s\n"%(ID,k,Probes[k][0],str(Probes[k][1])))   # Write:ProbeID_NumProbe [StartNuc EndNuc ProbeScore]
    
def getNucScore(blastfile):
    
    ## Parsing the Blast File
    inf = open("%s/%s"%(path,blastfile))
    parser = NCBIStandalone.BlastParser()
    error_parser = NCBIStandalone.BlastErrorParser(inf)
    iterator = NCBIStandalone.Iterator(inf, parser)
    err_iterator = NCBIStandalone.Iterator(inf, error_parser)
    for record in iterator:
        GenMatrix=dict.fromkeys(range(1,record.query_letters+1))
        PosMatrix=defaultdict(list)         ## this dic has the form: GenePos:[al1(1/0),...,ali(1/0)]; It determine the position specificity
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                if hsp.identities[1] > o.min_len:
                    if hsp.identities[1] > (o.tot_cov*record.query_letters) and round(float(hsp.identities[0])/float(hsp.identities[1]),3) >= o.max_sim:pass
                    else:
                        if hsp.identities[1] > (o.max_len*record.query_letters) and round(float(hsp.identities[0])/float(hsp.identities[1]),3) >= 0.99:pass
                        else:
                            rel_pos=0
                            for index, i in enumerate(hsp.query):     ## Running throught alignment
                                if i == ' ' or i == '-':pass        ## Jumping the gaps
                                else:
                                    pos=rel_pos+hsp.query_start
                                    if hsp.match[index] == '|':PosMatrix[pos].append(1)     ## If there is a match 1 --> not specificity
                                    elif hsp.match[index] == ' ':PosMatrix[pos].append(0)   ## If there is a mismatch 0 --> nucleotide specific
                                    rel_pos+=1
    hits=[]
    for key in PosMatrix.keys(): 
        hits.append(sum(PosMatrix[key]))
    return "%s\t%s"%(max(hits),np.mean(hits))
    #print GenMatrix
    #for key in GenMatrix.keys():
        ##try:GenMatrix[key]=round(float(sum(PosMatrix[key]))/float(max(hits)),4) # PIS for each position
    #    PIS = round(float(sum(PosMatrix[key]))/max(hits),4)
    #    GenMatrix[key] = PIS
    #return GenMatrix    ## Matrix for each position of gene secuence with it Position un-specificity Score (0 - 1) where 0 is specific and 1 not specific
    
def getOligoScore(scores):
    '''This function calculate the Oligo specificity score SS and it orders oligos according to SS
    Input: Matrix with positional or nucleotide specificity Score (PIS) (0 - 1) where 0 is specific and 1 not specific; PIS for each nucleotide of the gene sequence
    Output: Dictionary with each posible oligo with its SS socre and position from gene sequence
    '''
    #print scores
    #pos = range(1,len(scores)+1)
    Oligos ={}
    for i in xrange(len(scores)-o.probe_len+1):     ## Running each position of the gene sequence in order to calculate the probe score of specificity SS (0-100) 100 = specific, 0 = no specific 
        SSscore = []                                ## List of scores (SSscore) of each nucleotide that belong to the probe
        for j in range(i,i+o.probe_len): SSscore.append(scores[j+1])    ## Appending nucleotide score to SSscore list
        SS = round (100*(1 - sum(SSscore)/o.probe_len),3)               ## Calculus of Probe specificity score SS using each nucleotide scores that compound the probe  
        Oligos[i] = [i,i+o.probe_len,SS]                                ## Puting each probe with its info.[start_pos,end_pos,specificity_score(SS)]
    
    Oligos = OrderedDict(sorted(Oligos.items(), key=lambda t: t[1][2],reverse=True))    # ordering probes by SS --> higher to lower specificity
    return Oligos

def getProbes(Oligos,seq,np,Probes):
    ''' This function verify the unambiguosity of the sequence and return the bests X probes for the next step, the Denoising.
    Input: All posible Oligos, sequence of the gen, number of probes needed, definitive Pobres
    Output: Dictionary with the X best probes unambious and with the overlapping allowed
    '''
    #print np,Probes
    if Oligos == {}:return Probes
    
    oligo=Oligos[Oligos.keys()[0]]
    oseq = seq[oligo[0]:oligo[1]]
    del Oligos[Oligos[Oligos.keys()[0]][0]]
    
    if np == 0:return Probes
    else:
        if re.search("[NRYSWKMBDHV]",str(oseq)) != None:return getProbes(Oligos,seq,np,Probes)     ## Evaluating ambiuosity of the oligo sequence; If it is an ambiguos seq then check the other oligo
        else:
            if overlaping(oligo,Probes,o.probe_num-np):         ## Evaluating oligo overlaping                
                Probes[o.probe_num-(np-1)] = [oligo,oseq]       ## One more probe  
                return getProbes(Oligos,seq,np-1,Probes)     ## Recursion if more probes are needed. --> according to number of probes parameter
            else:return getProbes(Oligos,seq,np,Probes)      ## Next probe, because oligo overlaping evaluation was failed

def overlaping(oligo,refProb,n):
    ''' Overlaping evaluation of the probes.
    Input: oligo and List of Probes for the gene
    The function evaluate if there is a bigger overlaping than the allowed
    Return True if the overlaping between the oligo and the probes are lower than the threshold otherwise return False
    '''
    if n == 0: return True
    
    ## Evaluate for each probe if absolute value of (OligoEnd - Probe end) --> Non overlaping region -- is bigger than the
    ## threshold, calculated with (ProbeLenght - OverlappingAllowed); then return True.
    else:  
        if abs(oligo[1]-refProb[n][0][1]) >= int(o.probe_len)-int(o.p_overlap): return overlaping(oligo,refProb,n-1) ## recursion to evaluate the oligo with all probes
        else:return False

def probeDenoising(inf):
    print inf
    fx,uc = mkstemp(dir=path)
    seeds = "usearch --cluster %s --id %s --uc %s.uc --log %s.log" % (inf,o.clus_id,inf,inf) #--seedsout probes.fa
    uc2fasta = "usearch --uc2fasta %s --input %s --output %s" % (uc, inf,o.out_file_name)
    cluster = sub.Popen(seeds, shell = True, close_fds = True, stdout=sub.PIPE)
    lc=0
    if cluster.wait() == 0:
        fl = int(open("%s.uc"%inf,'r').read().count("\n")) # length of file
        if o.verbose==True:sys.stderr.write("Denoising...\n")      
        for i in open("%s.uc"%inf,'r'):
            if re.search("^C",i):
                if int(i.split('\t')[2]) == 1:os.write(fx,i)
            lc+=1
            if o.verbose==True:
                sys.stderr.write('\r'+'' *0)
                sys.stderr.write(str(int(lc*100/fl))+'%')
                sys.stdout.flush()       
        os.close(fx)
        uc2fa = sub.Popen(uc2fasta, shell = True, close_fds = True, stdout=sub.PIPE) 
        if uc2fa.wait() == 0:
            os.remove(uc)
            os.system("rm -r %s.log %s"%(inf,inf))
            sys.stderr.write('All done!; Probes designed in file %s\n'%o.out_file_name)
            

if __name__ == "__main__": main()

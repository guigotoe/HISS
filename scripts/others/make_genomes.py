#!/usr/bin/python3

# Guillermo Torres MSc,
# ggtorrese@unal.edu.co, guigotoe@gmail.com
# August 2012
# This script make random genome with one gene of our interest.
# ********************* How to use ***************************
# make_genomes fastainfile fastaoutfile 
# 
# ************************************************************
###################
     
import random, screed, math

def main():
    random.seed(1)                  # make reproducible 
    
    infname = "NCDScomplete_p.fasta" #infname = sys.argv[1]
    outname = "NCDScomplete_p_G" #outname = sys.argv[2]
    
    outfile='%s.fasta' % outname
    outf=open(outfile,"w+")
    logfile='%s.log' % outname
    outlog=open(logfile,"w+")
    i=0
    
    for g in screed.open("%s" % infname):
        name = g.name
        gen = g.sequence.lower()
        power = int(g.description)
    
        GENOMELEN = int(1e6)
        GENLEN = len(gen)
    
        genome = generate_genome(GENOMELEN)
        pos = random.randint(1, GENOMELEN) - 1
        genome = genome[:pos] + gen + genome[pos+GENLEN:]
        
        outf.write(">genome%d %d\n%s\n" % (i, power, genome))
        outlog.write("genome%d for sequence => %s >> %d to %d\n" % (i, name, pos, (pos+GENLEN)))
        i+=1

def generate_genome(length):
    x = ["A"] + ["G"] + ["C"] + ["T"]
    y = int(length/4)
    x = x*y
    
    random.shuffle(x)
    return "".join(x)       

if __name__ == "__main__": main()
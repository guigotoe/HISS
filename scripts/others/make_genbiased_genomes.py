#!/usr/bin/env python

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
    outname = "NCDScomplete_p_10gxG" #outname = sys.argv[2]
    
    outfile='%s.fasta' % outname
    outf=open(outfile,"w+")
    logfile='%s.log' % outname
    outlog=open(logfile,"w+")
    
    N_GENES = 10
    N_GENOMES = 100
    GENOMELEN = int(1e6)
    POWER_RANGE = 3
    
    seqs = {}
    indices =[]
    index = 0
    
    for g in screed.open("%s" % infname):
        name = g.name
        gen = g.sequence.lower()
        GENLEN = len(gen)
        power = g.description
        seqs[index] = [gen, name, power]
        indices += [index]
        index += 1
        
        G_id = 0
        
    for i in range(N_GENOMES):
        count = N_GENES #count = random.randint(1, N_GENES)
        genome = generate_genome(GENOMELEN)
        stp = []
        for i in range(count):
            hit = random.choice(indices)
            genseq = seqs[hit][0]
            genlen = len(genseq)
            genname = seqs[hit][1]
            genpower = seqs[hit][2]
            genomepower = random.choice(range(1, POWER_RANGE + 1))
            
            pos = ubication(stp, genlen, GENOMELEN)
            genome = genome[:pos] + genseq + genome[pos+GENLEN:]
            stp.append(pos)
            outlog.write("genome%d for sequence => %s >> %d to %d\n" % (G_id, genname, pos, (pos+genlen)))

        outf.write(">genome%d_%d %d\n%s\n" % (G_id, count, genomepower, genome))
        
            
        G_id += 1

def generate_genome(length):
    x = ["A"] + ["G"] + ["C"] + ["T"]
    y = int(length/4)
    x = x*y
    
    random.shuffle(x)
    return "".join(x)

def ubication (stp, GENLEN, GENOMELEN):
    pos = random.randint(1, GENOMELEN - GENLEN) - 1
    for p in stp:
        if p<pos<(p+GENLEN) or p<(pos+GENLEN)<(p+GENLEN):
            pos = ubication(stp, GENLEN, GENOMELEN)
    return pos

if __name__ == "__main__": main()
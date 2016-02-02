#!/usr/bin/env python
# Guillermo Torres MSc,
# ggtorrese@unal.edu.co, guigotoe@gmail.com
# August 2012
# This script make random reads from data set. This generates homogeneous data 
# in order to sequence coverage
# ********************* How to use ***************************
# make_reads  dataset.fasta
# 
# ************************************************************
###################

import screed, sys, random

def main ():
    random.seed(1)                  # make this reproducible.

    infname = "NCDScomplete_p_10gxG.fasta" #infname = sys.argv[1]
    outname = "NCDScomplete_p_10gxG_10X" #outname = sys.argv[2]

    outfile='%s.fasta' % outname
    outf=open(outfile,"w+")
    logfile='%s.log' % outname
    outlog=open(logfile,"w+")

    COVERAGE=10
    READLEN=100
    ERROR_RATE=100

    #record = iter(screed.open(sys.argv[1])).next()
    #record = iter(screed.open('NCDScomplete_p.fasta')).next()
    record = screed.open(infname)
    id = 0
    total_reads = 0
    total_mutated = 0
    total_mutations = 0

    for g in record:
        name = "n%s" % g.name # this is for noise # g.name 
        genome = g.sequence
        len_genome = len(genome)
    
        n_reads = int(len_genome*COVERAGE / float(READLEN))
        reads_mut = 0
        total_mut = 0
    
        for i in range(n_reads):
            if len_genome < 100:
                pass
            else:
                start = random.randint(0, len_genome - READLEN)
                read = genome[start:start + READLEN].upper()

    # reverse complement?
            if random.choice([0, 1]) == 0:
                read = screed.rc(read)

    # error?
            was_mut = False
            for _ in range(READLEN):
                while random.randint(1, ERROR_RATE) == 1:
                    pos = random.randint(1, READLEN) - 1
                    read = read[:pos] + random.choice(['a', 'c', 'g', 't']) + read[pos+1:]
                    was_mut = True
                    total_mut += 1

            if was_mut:
                reads_mut += 1
            
            outf.write('>n%dread%d\n%s\n' % (id,i, read))

        #print >>sys.stderr, "%d of %d reads mutated; %d total mutations from sequence %s" % \
        #(reads_mut, n_reads, total_mut, name)
        outlog.write("%d of %d reads mutated; %d total mutations from sequence %d => %s\n" % \
                   (reads_mut, n_reads, total_mut, id, name))
        total_reads += n_reads
        total_mutated += reads_mut
        total_mutations += total_mut
        id += 1

    outlog.write("TOTAL: %d of %d reads mutated; %d total mutations" % \
                 (total_mutated, total_reads, total_mutations))
    print >>sys.stderr, "Done!!"
    
if __name__ == "__main__": main()
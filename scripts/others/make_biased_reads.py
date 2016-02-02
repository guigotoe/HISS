#!/usr/bin/env python
# Based on diginorm script by Titus Brown 
# Guillermo Torres MSc,
# ggtorrese@unal.edu.co, guigotoe@gmail.com
# August 2012
# This script make random reads from data set.
# ********************* How to use ***************************
# make_reads  dataset.fasta
# 
# ************************************************************
###################

import screed, sys, random, math

#random.seed(1)                  # this generates the same result every time.

outfile='reads3.fasta'
outf=open(outfile,"w")

#COVERAGE=20
#N_READS = int(len_genome*COVERAGE / float(READLEN))
N_READS = int(1e6)
READLEN=100
ERROR_RATE=100

indices = []
seqs = []
powers = {}
z = []

index = 0
#for r in screed.open(sys.argv[1]):
for r in screed.open('dataset1.fasta'):
    power = int(r.description)
    count = int(math.pow(10, power))
    indices += [index] * count
    seqs.append(r.sequence)
    powers[index] = power

    index += 1

n_reads = N_READS
reads_mut = 0
total_mut = 0

z = []
for i in range(n_reads):
    index = random.choice(indices)
    sequence = seqs[index]
    start = random.randint(0, len(sequence) - READLEN)
    read = sequence[start:start + READLEN].upper()

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

    #print ('>read%d %d\n%s' % (index, i, read))
    outf.write('>read%d %d\n%s\n' % (index, i, read))
    z.append(index)

y = []
for i in set(z):
    y.append((i, z.count(i)))       # information about read from (by index) and abundance of it.
y.sort()
print >>sys.stderr, y

print >>sys.stderr, "%d of %d reads mutated; %d total mutations" % \
    (reads_mut, n_reads, total_mut)
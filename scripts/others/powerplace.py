#!/usr/bin/env python
# Guillermo Torres MSc,
# ggtorrese@unal.edu.co, guigotoe@gmail.com
# August 2012
# This script put a power descriptor in fasta comment line.
# ********************* How to use ***************************
# make_reads  dataset.fasta
# 
# ************************************************************
###################


import sys, re, random


random.seed(1)                  # make this reproducible, please.

outfile='NCDScomplete_p.fasta'
outf=open(outfile,"w")

POWER_RANGE = 3

for g in open('NCDScomplete.fasta'):
    g = g.strip("\n")
    if re.match('>', g):
        power = random.choice(range(1, POWER_RANGE + 1))
        outf.write("%s %d\n" % (g.strip("\n"), power))
    else:
        if re.search('r', g): g=re.sub('r', random.choice(['a', 'g']), g) 
        if re.search('y', g): g=re.sub('y', random.choice(['c', 't']), g)
        if re.search('m', g): g=re.sub('m', random.choice(['a', 'c']), g)
        if re.search('k', g): g=re.sub('k', random.choice(['g', 't']), g)
        if re.search('s', g): g=re.sub('s', random.choice(['c', 'g']), g)
        if re.search('w', g): g=re.sub('w', random.choice(['a', 't']), g)
        if re.search('h', g): g=re.sub('h', random.choice(['a', 'c', 't']), g)
        if re.search('b', g): g=re.sub('b', random.choice(['c', 'g', 't']), g)
        if re.search('v', g): g=re.sub('v', random.choice(['a', 'c', 'g']), g)
        if re.search('d', g): g=re.sub('d', random.choice(['a', 'g', 't']), g)
        if re.search('n', g): g=re.sub('n', random.choice(['a', 'c', 'g', 't']), g)
        outf.write("%s\n" % (g.strip("\n")))
print ('Done!!')

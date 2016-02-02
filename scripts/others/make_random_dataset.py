#!/usr/bin/python3

# Guillermo Torres MSc,
# ggtorrese@unal.edu.co, guigotoe@gmail.com
# August 2012
# This script make random data set.
# ********************* How to use ***************************
# make_random_dataset Number_Of_Sequences Sequence_Length 
# 
# ************************************************************
###################
     
import random

def main():
          
    random.seed(1)                  # make random sequences
    
    N_SEQUENCES = 20000
    #SEQUENCES_LENGTH = 500
    POWER_RANGE = 3                   # for subsequent abundance calculation
    outfile='noise.fasta'
    outf=open(outfile,"w+")
    
    for i in range(N_SEQUENCES):
        SEQUENCES_LENGTH = random.randint(500, 1000)
        transcript = generate_sequence(SEQUENCES_LENGTH)
        power = random.choice(range(1, POWER_RANGE + 1))
        outf.write(">noise%d %d\n%s\n" % (i, power, transcript))
    
def generate_sequence(length):
    x = ["A"] + ["G"] + ["C"] + ["T"]
    y = int(length/4)
    x = x*y

    random.shuffle(x)
    return "".join(x)        

if __name__ == "__main__": main()
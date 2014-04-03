# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 10:58:43 2013

@author: Daniel Struck
"""
from Bio import SeqIO
from itertools import permutations
import random

ref = "Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
alphabet = ['A','C','G','T']

sequences = {}
with open("HIV1_REF_2010_genome_DNA_PURE.fasta") as f:
    for record in SeqIO.parse(f, "fasta") :
        sequences[record.id] = record.seq
    
noise_levels = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
window_sizes = [50,100,200,300,400]
stepping_sizes = [50,100,200,300,400]

hxb2_start = 790 # gag start
hxb2_end = 9417 # nef end

# introduce noise
def introduce_noise(input_seq):
    
    seq = list(input_seq)
    noise = random.choice(noise_levels)
    
    for i in range(int(len(str(input_seq).replace("-","")) * noise)):
        while True:
            pos_mut = random.randrange(len(seq))
            orig_nuc = seq[pos_mut]
            if orig_nuc != '-': break # do not mutate gaps
        seq[pos_mut] = random.choice([x for x in alphabet if x != orig_nuc])
    return seq


count = 0
with open("2013-12-20_recomb.fasta","w") as w:
    
    for perm in permutations(sequences.keys(), 2):

        stepping_size = random.choice(stepping_sizes)
        window_size = random.choice(window_sizes)

        background, insert = perm[0], perm[1]

        if background.split(".")[1] == insert.split(".")[1]:
            continue
        
        seq_background = introduce_noise(sequences[background]) # introduce random amount of noise 
        seq_insert = introduce_noise(sequences[insert]) # introduce random amount of noise 

        for i in xrange(hxb2_start, hxb2_end - window_size, stepping_size):

            # randomly skip 10%, else the file gets too big!
            if random.random() > 0.10:
                continue

            pos = 1
            new_seq = []
            for j in xrange(len(sequences[ref])):
                if pos >= i and pos < i + window_size:
                    new_seq.append(seq_insert[j])
                elif pos >= hxb2_start and pos <= hxb2_end:
                    new_seq.append(seq_background[j])
                if sequences[ref][j] != "-": pos += 1
                
            w.write(">{},URF,{},{},{},{}\n".format(count, background.split(".")[1], insert.split(".")[1], window_size, stepping_size))
            w.write("".join(new_seq).replace("-","") + "\n") # remove the gap symbol, else the file gets too big
            count += 1
    
print count

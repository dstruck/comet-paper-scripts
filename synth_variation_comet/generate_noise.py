# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:17:38 2013

@author: Daniel Struck

generate recombined sequence to be tested by subtyping tools

"""
from __future__ import division
from Bio import SeqIO
import random

ref = "Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
alphabet = ['A','C','G','T']
sequences = {}
with open("HIV1_REF_2010_genome_DNA.fasta") as f:
    for record in SeqIO.parse(f, "fasta") :
        sequences[record.id] = str(record.seq)
    
parameters = [(100,50),(200,100),(400,200),(600,300),(800,400),(1200,600),(1600,800)] # window_size, stepping_size
noise_levels = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20]
subtypes_to_exclude = ["O","N","P","CPZ"]

# introduce noise
dataset = {}
print "start noise generation"
for noise in noise_levels:
    temp_sequences = {}
    for name in sequences:
        if name.split(".")[1] in subtypes_to_exclude: continue
        seq = [x for x in sequences[name]] # put the nucleotides of a sequence in an array
        for i in range(int(len(sequences[name].replace("-","")) * noise)):
            while True:
                pos_mut = random.randrange(len(sequences[name]))
                orig_nuc = sequences[name][pos_mut]
                if orig_nuc != '-': break  # do not mutate gaps
            new_nuc = random.choice([x for x in alphabet if x != orig_nuc])
            seq[pos_mut] = new_nuc
        temp_sequences[name] = "".join(seq)
    dataset[noise] = temp_sequences
print "finished noise generation"
    
hxb2_start = 790 # gag start
hxb2_end = 9417 # nef end

count = 0
print "start writing down the sequences"
with open("dataset/2014-01-09_noise.fasta","w") as w:
    for noise in noise_levels:
        print "noise",noise
        for param in parameters:
            window_size = param[0]
            stepping_size = param[1]
            for name in dataset[noise]:
                for start_point in xrange(hxb2_start, hxb2_end - window_size, stepping_size):
                    pos = 1
                    new_seq = []
                    gap_count = 0
                    for (ref_nuc, query_nuc) in zip(dataset[noise][ref], dataset[noise][name]):
                        if pos == start_point:
                            # no larger gaps then the stepping size / 2
                            if len(new_seq) >= window_size or gap_count > (stepping_size / 2):
                                break
                            elif query_nuc != '-':
                                new_seq.append(query_nuc)
                            else:
                                gap_count += 1
                        elif ref_nuc != "-": pos += 1
                    if window_size - len(new_seq) <= 3: # maximum 3 NUC missing
                        w.write(">" + name + "_n" + str(noise) + "_w" + str(window_size) 
                        + "_s" + str(stepping_size)  + "_s" + str(start_point) +  "\n")
                        w.write("".join(new_seq) + "\n")
                        count += 1
        

print "finished"
print "count",count

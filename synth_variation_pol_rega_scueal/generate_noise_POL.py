# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:07:46 2013

@author: daniel
"""
from __future__ import division
from Bio import SeqIO
import random

ref = "Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
alphabet = ['A','C','G','T']

subtypes_to_remove = []
# SCUEAL has no reference sequences for those subtypes
subtypes_to_remove.append("26_AU")
subtypes_to_remove.append("46_BF")
subtypes_to_remove.append("42_BF")
subtypes_to_remove.append("45_cpx")
subtypes_to_remove.append("38_BF1")
subtypes_to_remove.append("44_BF")
subtypes_to_remove.append("49_cpx")
subtypes_to_remove.append("P")

# REGAv2 does not support thoses subtypes
# source: http://bioafrica.mrc.ac.za/CRFs/CRFs.php
subtypes_to_remove.append("15_01B")
subtypes_to_remove.append("16_A2D")
subtypes_to_remove.append("17_BF")
subtypes_to_remove.append("21_A2D")
subtypes_to_remove.append("22_01A1")
subtypes_to_remove.append("23_BG")
subtypes_to_remove.append("26_AU")
subtypes_to_remove.append("28_BF")
subtypes_to_remove.append("32_06A1")
subtypes_to_remove.append("33_01B")
subtypes_to_remove.append("34_01B")
subtypes_to_remove.append("36_cpx")
subtypes_to_remove.append("44_BF")
subtypes_to_remove.append("45_cpx")
subtypes_to_remove.append("46_BF")
subtypes_to_remove.append("48_01B")
subtypes_to_remove.append("CPZ")

# generally remove those subtypes
subtypes_to_remove.append("N")
subtypes_to_remove.append("O")

# load the LANL subtype reference sequences
subtypes = set()
sequences = {}
with open("HIV1_REF_2010_genome_DNA.fasta") as f:
    for record in SeqIO.parse(f, "fasta"):
        subtype = record.id.split(".")[1]
        if subtype in subtypes_to_remove: continue
        if subtype not in subtypes: # take only the first reference sequence to make the dataset smaller
            sequences[record.id] = str(record.seq)
            subtypes.add(subtype)
    
parameters = [(200,200),(400,200),(800,400),(1200,600),(1600,800)] # window_size, stepping_size
noise_levels = [0.0, 0.025, 0.05, 0.075, 0.1]

# introduce noise
print "start noise generation"
dataset = {}
for noise in noise_levels:
    temp_sequences = {}
    for name in sequences:
        seq = [x for x in sequences[name]]
        for i in range(int(len(sequences[name].replace("-","")) * noise)):
            while True:
                pos_mut = random.randrange(len(sequences[name]))
                orig_nuc = sequences[name][pos_mut]
                if orig_nuc != '-': break # do not mutate gaps
            seq[pos_mut] = random.choice([x for x in alphabet if x != orig_nuc])
        temp_sequences[name] = "".join(seq)
    dataset[noise] = temp_sequences
print "finished noise generation"
    
hxb2_start = 2085 # pol start
hxb2_end = 5096  # pol end

count = 0
print "start writing down the sequences"
with open("dataset/2013-08-29_noise_POL.fasta","w") as w, open("dataset/2013-08-29_noise_POL_names.csv","w") as n:
    n.write("\t".join(("id","name","noise","window_size","stepping_size","start_pos"))+"\n")
    for noise in noise_levels:
        print "noise",noise
        for param in parameters:
            window_size = param[0]
            stepping_size = param[1]
            for name in dataset[noise]:
                for i in xrange(hxb2_start, hxb2_end - window_size, stepping_size):
                    pos = 1
                    new_seq = []
                    for (ref_nuc, query_nuc) in zip(dataset[noise][ref], dataset[noise][name]):
                        if pos >= i and pos < i + window_size:
                            new_seq.append(query_nuc)
                        if ref_nuc != "-":
                            pos += 1

                    w.write(">" + str(count) +  "\n")
                    w.write("".join(new_seq) + "\n")

                    n.write("\t".join((str(count),name,str(noise),str(window_size), 
                    str(stepping_size),str(i)))  +  "\n")
                    count += 1


print "finished"
print "count",count
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:07:46 2013

@author: Daniel Struck
"""
from __future__ import division
from Bio import SeqIO
import random

ref = "Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
alphabet = ['A','C','G','T']

subtypes_to_remove = []
subtypes_to_remove.append("P")
subtypes_to_remove.append("CPZ")
subtypes_to_remove.append("N")
subtypes_to_remove.append("O")

# load the LANL subtype reference sequences
#subtypes = set()
sequences = {}
with open("HIV1_REF_2010_genome_DNA.fasta") as f:
    for record in SeqIO.parse(f, "fasta"):
        subtype = record.id.split(".")[1]
        if subtype in subtypes_to_remove: continue
        sequences[record.id] = str(record.seq)
    
noise_levels = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]

# introduce noise
print "start noise generation"
dataset = {}
for noise in noise_levels:
    for j in range(10):
        temp_sequences = {}
        for name in sequences:
            seq = [x for x in sequences[name]]
            for i in range(int(len(sequences[name].replace("-","")) * noise)):
                while True:
                    pos_mut = random.randrange(len(sequences[name]))
                    orig_nuc = sequences[name][pos_mut]
                    if orig_nuc != '-': break # do not mutate gaps
                new_nuc = random.choice([x for x in alphabet if x != orig_nuc])
                seq[pos_mut] = new_nuc
            temp_sequences[name] = "".join(seq)
        dataset["{}-{}".format(noise,j)] = temp_sequences
print "finished noise generation"
    
hxb2_start = 790 # gag start, 790
hxb2_end = 9417 # nef end

count = 0
print "start writing down the sequences"
with open("2013-12-18_noise_param.fasta","w") as w:
    for k in dataset.keys():
        print "key",k
        for name in dataset[k]:
            pos = 1
            new_seq = []
            for (ref_nuc, query_nuc) in zip(dataset[k][ref], dataset[k][name]):
                if pos >= hxb2_start and pos <= hxb2_end:
                    new_seq.append(query_nuc)
                if ref_nuc != "-":
                    pos += 1
            w.write(">" + ",".join((str(count),name,str(k))) +  "\n")
            w.write("".join(new_seq) + "\n")
            count += 1


print "finished"
print "count",count

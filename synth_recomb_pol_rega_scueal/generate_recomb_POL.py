# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 10:06:57 2013

@author: Daniel Struck
"""
from Bio import SeqIO
from itertools import permutations

ref = "Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455"

sequences = {}
subtypes = set()
with open("HIV1_REF_2010_genome_DNA_PURE.fasta") as f:
    for record in SeqIO.parse(f, "fasta") :
        subtype = record.id.split(".")[1]
        if subtype not in subtypes: # take only the first ref sequence to keep the dataset smaller for REGA and SCUEAL
            sequences[record.id] = str(record.seq)
            subtypes.add(subtype)
    
parameters = [(100,100),(200,200),(300,200),(400,200),(600,200),(800,200)] # window_size, stepping_size

hxb2_start = 2085 # pol start
hxb2_end = 5096  # pol end

count = 0
with open("dataset/2013-10-18_recomb_POL.fasta","w") as w, open("dataset/2013-10-18_recomb_POL_names.csv","w") as n:
    n.write("\t".join(("id","background","insert","window_size","stepping_size","start_pos"))+"\n")
    for param in parameters:
        window_size = param[0]
        stepping_size = param[1]
        for perm in permutations(sequences.keys(), 2):
            background = perm[0]
            insert = perm[1]
            if background.split(".")[1] == insert.split(".")[1]:
                continue
            for i in xrange(hxb2_start, hxb2_end - window_size, stepping_size):
                pos = 1
                new_seq = []
                for (ref_nuc, back_nuc, ins_nuc) in zip(sequences[ref], sequences[background], sequences[insert]):
                    if pos >= i and pos < i + window_size:
                        new_seq.append(ins_nuc)
                    elif pos >= hxb2_start and pos <= hxb2_end:
                        new_seq.append(back_nuc)
                    if ref_nuc != "-": pos += 1
                    
                w.write(">" + str(count) + "\n")
                w.write("".join(new_seq).replace("-","") + "\n") # remove the gap symbol, else the file gets too big
                n.write("\t".join((str(count), background, insert, str(window_size), str(stepping_size), str(i))) + "\n")
                count += 1
    
print count

# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:17:38 2013

@author: Daniel Struck

generate recombined sequence to be tested by subtyping tools

"""
from __future__ import division
from Bio import SeqIO
from itertools import combinations
import numpy as np

ref = "Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455"

sequences = {}
with open("HIV1_REF_2010_genome_DNA.fasta") as f:
    for record in SeqIO.parse(f, "fasta") :
        sequences[record.id] = record.seq
    
hxb2_start = 790 # gag start
hxb2_end = 9417  # nef end

subtypes_to_exclude = ["O","N","P","CPZ"]

result = []
for comb in combinations(sequences.keys(), 2):

    ref1, ref2 = comb[0], comb[1]
    subtype1, subtype2 = ref1.split(".")[1], ref2.split(".")[1]
    
    if subtype1 in subtypes_to_exclude or subtype2 in subtypes_to_exclude: continue
    
    # do not include CRFs
    if "_" in subtype1 or "_" in subtype2: continue

    # differences within subtypes    
    if subtype1 != subtype2:
        continue

    pos = 1
    count = 0
    count_diff = 0
    for (ref_nuc, back_nuc, ins_nuc) in zip(sequences[ref], sequences[ref1], sequences[ref2]):
        if pos >= hxb2_start and pos <= hxb2_end:
            if back_nuc != ins_nuc:
                count_diff += 1
            count += 1
        if ref_nuc != "-": pos += 1
    diff_pct =  count_diff / count * 100
    print subtype1,subtype2,diff_pct
    result.append(diff_pct)

print "---"
print "min:",min(result)
print "max:",max(result)
print "avg:",sum(result)/len(result)
print "std:",np.std(result)

# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 12:30:26 2013

@author: Daniel Struck

create blast nucleotide database from the COMET training set with the following command:
makeblastdb  -in training_subtype_03.fasta -dbtype nucl

run blast with the following command:
(blast version: 2.2.29)
blastn -num_threads 12 -db training_subtype_03.fasta -query 2014-01-09_noise.fasta -num_descriptions 1 -num_alignments 0 > blast_results.txt

"""
from __future__ import division
from collections import Counter

n = 0
window_n = Counter()
noise_n = Counter()
pos = 0
window_pos = Counter()
noise_pos = Counter()
query = ""

with open("blast_results.txt") as f:
    for line in f:
        if line.startswith("Query"):
            
            query =  line.split(".")[1]
            token = line.split("_")
            noise = float(token[-4][1:])
            window_size = int(token[-3][1:])
            
            window_n[window_size] += 1
            noise_n[noise] += 1
            n += 1
            
            for _ in range(6): line = next(f) # forward to result

            result = line.split(".")[0].lstrip()
            if query == result:
                pos += 1
                window_pos[window_size] += 1
                noise_pos[noise] += 1

            
print "\tdetected\tn\tsensitivity"           
print "total:\t{}\t{}\t{}".format(pos, n, pos / n * 100)

print

for noise in sorted(noise_n.keys()):
    print "noise {}:\t{}\t{}\t{}".format(noise, noise_pos[noise], noise_n[noise], noise_pos[noise] / noise_n[noise])
 
print
   
for window in sorted(window_n.keys()):
    print "window {}:\t{}\t{}\t{}".format(window, window_pos[window], window_n[window], window_pos[window] / window_n[window])    
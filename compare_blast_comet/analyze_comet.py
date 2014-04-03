# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 12:30:26 2013

@author: Daniel Struck
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

with open("2014-01-09_noise.fasta.csv") as f:
    next(f) # header
    for line in f:
        token = line.split("\t")
        query = token[0].split(".")[1]
        result = token[-2].rstrip()

        token2 = token[0].split("_")
        noise = float(token2[-4][1:])
        window_size = int(token2[-3][1:])

        n += 1
        window_n[window_size] += 1
        noise_n[noise] += 1
        
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
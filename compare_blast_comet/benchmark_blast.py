# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 09:53:51 2014

@author: Daniel Struck
"""
import os, time

for i in range(11):
    start_time = time.time()
    os.system(".blastn -num_threads 12 -db training_subtype_03.fasta -query /home/daniel/2014-01-09_noise.fasta -num_descriptions 1 -num_alignments 0 > blast_results.txt");
    print "{}. elapsed {}".format(i, time.time() - start_time)

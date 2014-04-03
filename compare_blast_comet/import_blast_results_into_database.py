# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 12:30:26 2013

@author: Daniel Struck

create blast results for the COMET <> BLAST LANL comparison:
blastn -num_threads 12 -db training_subtype_03.fasta -query ../lanl_comparison/dataset/2013-08-30_pol_lanl.fasta -num_descriptions 1 -num_alignments 0 > blast_lanl_pol.txt 

import BLAST results into the database

for the comparison see ../lanl_comparison

"""
from __future__ import division
import psycopg2

conn = psycopg2.connect(database='comet')
cur = conn.cursor()

with open("blast_lanl_pol.txt") as f:
    for line in f:
        if line.startswith("Query"):
            
            db_id =  int(line.split(" ")[-1])
            for _ in range(6): line = next(f)
            subtype = line.strip().split(" ")[0].split(".")[0]
            cur.execute("UPDATE sequences SET subtype_blast = %s WHERE id = %s",
                    (subtype, db_id))
                    
conn.commit()
cur.close()
conn.close
print "finished"

            



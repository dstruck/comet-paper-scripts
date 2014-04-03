# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 12:47:37 2013

@author: Daniel Struck
"""
from __future__ import division
import psycopg2 

conn = psycopg2.connect(database='comet')
cur = conn.cursor()

subtype = "11_CPX"

ids = []
with open("2013-08-30_pol_lanl.csv") as f:
    next(f)
    for line in f:
        token = line.rstrip().split("\t")
        if token[2].upper() == subtype:
            ids.append(token[0])


cur.execute("SELECT id, subtype_comet, subtype_scueal, subtype_rega, sequence FROM sequences WHERE id IN %s",
                (tuple(ids),))
for row in cur:
    if row[1].upper() == subtype:
        print "{}\t{}\t{}\t{}".format(row[0],row[1].replace(",","_").replace(";","_"),row[2].replace(",","_"),row[3].replace("HIV-1 Subtype ",""))


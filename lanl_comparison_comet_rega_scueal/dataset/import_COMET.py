# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 11:37:05 2013

@author: Daniel Struck

the dump of the complete database can be found in: database_dump.sql.bz2

"""
import psycopg2

conn = psycopg2.connect(database='comet')
cur = conn.cursor()

#with open('2013-08-30_pol_lanl.fasta.major.csv') as f:
with open('2013-08-30_pol_lanl.fasta.csv') as f:
    next(f) #header
    for line in f:
        token = line.rstrip().split("\t")
        name = token[0]
        subtype = token[2]
#        cur.execute("UPDATE sequences SET subtype_comet_major = %s WHERE id = %s",
        cur.execute("UPDATE sequences SET subtype_comet = %s WHERE id = %s",
                    (subtype, name))
        if cur.rowcount != 1:
            print "should not happen"
                        
conn.commit()
cur.close()
conn.close
print "finished"

            



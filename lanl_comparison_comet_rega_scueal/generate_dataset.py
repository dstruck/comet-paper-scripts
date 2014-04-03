# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 14:49:27 2013

@author: Daniel Struck
"""
import psycopg2 

conn = psycopg2.connect(database="comet")

count = 0    
with conn.cursor() as cur:
    cur.execute("""SELECT subtype_lanl,count(id) as anzahl 
                    FROM sequences 
                    WHERE gene = 'POL' 
                    AND LENGTH(sequence) > 800
                    AND subtype_LANL NOT IN ('-','A','U','A1C','02G','A1D','BF','BG','CD','A1G','02A1',
                                             '19B','0206','AD','O','BF1','BC','AG','01A1','F','01B')
                    GROUP BY subtype_lanl
                    HAVING count(id) >= 50
                    ORDER BY subtype_lanl; """)
    subtypes = []
    for row in cur:
        subtypes.append(row[0])
    
    from collections import Counter
    zaehlen = Counter()
    with open("dataset/2013-08-30_pol_lanl.fasta","w") as w, open("dataset/2013-08-30_pol_lanl.csv","w") as n:
        n.write("id\tname\tsubtype_lanl\n")
        for subtype in subtypes:
            cur.execute("""SELECT id, sequence, name, subtype_lanl 
                            FROM sequences 
                            WHERE gene = 'POL' 
                            AND subtype_LANL = %s
                            AND LENGTH(sequence) > 800
                            ORDER BY RANDOM()
                            LIMIT 1000; """,(subtype,))
                            
            for row in cur:
                zaehlen[row[3]] += 1
                w.write(">" + str(row[0]) + "\n")
                w.write(row[1] + "\n")
                n.write("\t".join((str(row[0]),row[2],row[3])) + "\n")
                count += 1
                
print "count:",count            
print "finished"

for entry in zaehlen:
    print entry, zaehlen[entry]
        
'''
comet=# SELECT min(length(sequence)),MAX(length(sequence))  from sequences where gene = 'POL';
 min | max  
-----+------
 289 | 4993
(1 row)
'''

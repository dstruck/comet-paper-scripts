# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:17:38 2013

@author: Daniel Struck

generate recombined sequence to be tested by subtyping tools

"""
from __future__ import division
from collections import Counter, defaultdict

simplify = {"A1":"A", "A2":"A", "F1":"F", "F2":"F"}

scueal_translation = {
'CRF02':'02_AG','CRF03':'03_AB','CRF04':'04_cpx','CRF05':'05_DF','CRF06':'06_cpx'
,'CRF07':'07_BC','CRF08':'08_BC','CRF09':'09_cpx','CRF10':'10_CD','CRF11':'11_CPX','CRF12':'12_BF'
,'CRF13':'13_cpx','CRF14':'14_BG','CRF15':'15_01B','CRF16':'16_A2D','CRF17':'17_BF'
,'CRF18':'18_cpx','CRF19':'19_cpx','CRF20':'20_BG','CRF21':'21_A2D','CRF22':'22_01A1'
,'CRF23':'23_BG','CRF24':'24_BG','CRF25':'25_cpx','CRF27':'27_cpx','CRF28':'28_BF'
,'CRF29':'29_BF','CRF31':'31_BC','CRF32':'32_06A1','CRF33':'33_01B','CRF34':'34_01B'
,'CRF35':'35_AD','CRF36':'36_cpx','CRF37':'37_cpx','CRF39':'39_BF','CRF40':'40_BF'
,'CRF43':'43_02G','CRF47':'47_BF','AE':'01_AE'}

'''
uncomment the tool and only one tool to analyze its result
'''

#filename_dataset = "dataset/SCUEAL/2013-10-18_recomb_POL.csv"
#def retrieve_query_subtype(line):
#    token = line.split("\t")
#    query, id = token[3], token[1]
#    result = False
#    subtypes_found = []
#    if query == 'U' or "recombinant" in query or "CRF" in query or "Complex" in query:
#        for _ in token[2].split(" ")[0].split(","):
#            if _ in scueal_translation: _ = scueal_translation[_]
#            if "_" in _ and "cpx" not in _:
#                [subtypes_found.append(x) for x in _.split("_")[-1]]
#            else:
#                subtypes_found.append(_)
#        result = True
#    return id, result, [(simplify[x] if x in simplify else x) for x in subtypes_found]

#filename_dataset = "dataset/COMET/2013-10-18_recomb_POL.fasta.csv"
#def retrieve_query_subtype(line):
#    token = line.split("\t")
#    query, id = token[2], token[0]
#    result = False
#    subtypes_found = []
#    if "unassigned" in query:
#        result = True
#        subtypes_found = query.replace(" ","").split(";")[-1].split(",")
#        for _ in subtypes_found:
#            if "_" in _ and "cpx" not in _:
#                [subtypes_found.append(x) for x in _.split("_")[-1]]
#    return id, result, [(simplify[x] if x in simplify else x) for x in subtypes_found]

filename_dataset = "dataset/REGAv2/all.csv"
def retrieve_query_subtype(line):
    token = line.rstrip().replace('"','').split(",")
    query, id = token[2], token[0]
    result = False
    subtypes_found = []
    if "Check" in query or "CRF" in query:
        result = True
        subtypes_found = token[15].replace(" ","-").split("-")
    return id, result, [(simplify[x] if x in simplify else x) for x in subtypes_found]


window_sizes = set()
start_positions = defaultdict(set) 

count_window_size = Counter()   
count_query_window_size = Counter()
count_query_window_size_correct = Counter()

count_subtype = Counter()   
count_query_subtype = Counter()

count_insert = Counter()   
count_query_insert = Counter()

count_sp_ws = Counter()   
count_query_sp_ws = Counter()

data = {}
#id      background      insert  window_size     stepping_size   start_pos
with open("dataset/2013-10-18_recomb_POL_names.csv") as f:
    next(f) #header
    for line in f:
        token = line.rstrip().split("\t")
        data[token[0]] = token[1:]

done = set()
with open(filename_dataset) as f:
    
    next(f) # header
    
    for line in f:

        id, query_test, subtypes_found = retrieve_query_subtype(line)
        if id in done: continue
        else: done.add(id) # had to resubmit batches to REGA several times, so some sequences have been subtyped by REGA several times
        
        windows_size = int(data[id][2])        
        background = data[id][0].split(".")[1]
        insert = data[id][1].split(".")[1]        
        start_pos_insert = int(data[id][4])
        
        start_positions[windows_size].add(start_pos_insert)
        window_sizes.add(windows_size)

        if query_test:
            count_query_window_size[windows_size] += 1
            count_query_subtype[background,insert] += 1
            count_query_insert[insert] += 1
            count_query_sp_ws[start_pos_insert, windows_size] += 1

            #did the tool find the right composition
            if background in simplify: background = simplify[background]
            if insert in simplify: insert = simplify[insert]
            if background in subtypes_found and insert in subtypes_found:
                count_query_window_size_correct[windows_size] += 1
            
        count_window_size[windows_size] += 1
        count_subtype[background, insert] += 1
        count_insert[insert] += 1
        count_sp_ws[start_pos_insert, windows_size] += 1


print "window_size\tn\tfound\tquery"
for e in sorted(count_window_size):
    print "\t".join([str(e)
                     ,str(count_window_size[e])
                     ,str(count_query_window_size[e])
                     ,str(count_query_window_size[e] / count_window_size[e])
                     ])
print        
print "window_size\tn\tfound\tcorrectly identified"
for e in sorted(count_window_size):
    print "\t".join([str(e)
                     ,str(count_window_size[e])
                     ,str(count_query_window_size_correct[e])
                     ,str(count_query_window_size_correct[e] / count_window_size[e])
                     ])
print        
print "background, insert\tn\tquery"
for e in sorted(count_subtype):
    print "\t".join([str(e)
                     ,str(count_subtype[e])
                     ,str(count_query_subtype[e] / count_subtype[e])
                     ])
print    
print "insert\tn\tquery"
for e in sorted(count_insert):
    print "\t".join([str(e)
                     ,str(count_insert[e])
                     ,str(count_query_insert[e] / count_insert[e])
                     ])
print
print "n\tinsert start position\twindow size\tquery"
for ws in sorted(window_sizes):
    for sp in sorted(start_positions[ws]):
        if (sp,ws) in count_sp_ws:
            print "\t".join([str(count_sp_ws[sp,ws])
                             ,str(sp)
                             ,str(ws)
                             ,str(count_query_sp_ws[sp,ws] / count_sp_ws[sp,ws])
                             ])
print
for ws in sorted(window_sizes):
    print "window_size\tn\t" + "\t".join([str(x) for x in sorted(start_positions[ws])])
    print str(ws) + "\t" + str(count_window_size[ws]) + "\t" + "\t".join([str(count_query_sp_ws[x, ws]/count_sp_ws[x, ws]) for x in sorted(start_positions[ws])])


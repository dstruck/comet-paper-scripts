# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:17:38 2013

@author: Daniel Struck

generate recombined sequences to be tested by COMET

"""
from __future__ import division
from collections import Counter, defaultdict

simplify = {"A1":"A", "A2":"A", "F1":"F", "F2":"F"}

filename_dataset = "dataset/2014-01-09_recomb.fasta.csv"
def retrieve_query_subtype(line):
    token = line.split("\t")
    query, id = token[2], token[0]
    result = False
    subtypes_found = []
    
    if "unassigned" in query:
        result = True
        subtypes_found = query.replace(" ","").split(";")[-1].split(",")
        for _ in subtypes_found:
            if "_" in _ and "cpx" not in _:
                [subtypes_found.append(x) for x in _.split("_")[-1]]
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
with open("/mnt/data/bioinfo/comet/2014-01-09_recomb_names.csv") as f:
    next(f) #header
    for line in f:
        token = line.rstrip().split("\t")
        data[token[0]] = token[1:]

done = set()
with open(filename_dataset) as f:
    
    next(f) # header
    
    for line in f:

        id, query_test, subtypes_found = retrieve_query_subtype(line)
        
        windows_size = int(data[id][2])        
        background = data[id][0].split(".")[1]
        insert = data[id][1].split(".")[1]        
        start_pos_insert = int(data[id][4])
        
        start_positions[windows_size].add(start_pos_insert)
        window_sizes.add(windows_size)
        
        if query_test:

            #did COMET find the right composition
            if background in simplify: background = simplify[background]
            if insert in simplify: insert = simplify[insert]
            if background in subtypes_found and insert in subtypes_found:
                count_query_window_size_correct[windows_size] += 1
            else:
#                print background, insert, subtypes_found
                pass

            count_query_window_size[windows_size] += 1
            count_query_subtype[background,insert] += 1
            count_query_insert[insert] += 1
            count_query_sp_ws[start_pos_insert, windows_size] += 1
            
        count_window_size[windows_size] += 1
        count_subtype[background, insert] += 1
        count_insert[insert] += 1
        count_sp_ws[start_pos_insert, windows_size] += 1


print "window_size\tn\tfound\tsensitivity"
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
print "background, insert\tn\tsensitivity"
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
print "n\tinsert start position\twindow size\tsensitivity"
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

# generate plot
'''
import matplotlib.pyplot as plt

plots, labels = [], []

for ws in sorted(window_sizes):
    labels.append(ws)
    X = [x for x in sorted(start_positions[ws])]
    Y = [(count_query_sp_ws[x, ws]/count_sp_ws[x, ws]) for x in sorted(start_positions[ws])]
    plot, = plt.plot(X,Y,"o-")
    plots.append(plot)

plt.legend(plots, labels)
plt.title("synthetic recombination near full length: sens. / position")
plt.show()
'''
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:55:49 2013

@author: daniel
"""
from __future__ import division

pure = ['A1','A2','B','C','D','F1','F2','G','H','J','K']

def retrieve_query_subtype_comet(line):
    token = line.rstrip().split("\t")
    query_temp = token[2]
    query = []
    if "check" in query_temp:
        query.append(query_temp.split(" ")[0])
        query.append(query_temp.split(" ")[3][:-1])
    else:
        query.append(query_temp.upper())
    return query
    
    
noises, window_sizes, stepping_sizes, start_pos, subtypes, subtypes_pure, subtypes_crf = set(), set(), set(), set(), set(), set(), set()
def prepare_dataset(filename_dataset, retrieve_query_subtype):
    dataset = []
    with open(filename_dataset) as f:
        next(f) #header
        for line in f:
            ref_subtype = line.split("\t")[0].split(".")[1].upper()
            subtypes.add(ref_subtype)
            if ref_subtype in pure:
                subtypes_pure.add(ref_subtype)
            else:
                subtypes_crf.add(ref_subtype)

            query = retrieve_query_subtype(line)
            token = line.split("\t")[0].split("_")
            noise = float(token[-4][1:])
            noises.add(noise)
            ws = int(token[-3][1:])
            window_sizes.add(ws)
            ss = int(token[-2][1:])
            stepping_sizes.add(ss)
            sp = int(token[-1][1:])
            start_pos.add(sp)
            # ref_subtype, noise, window_size, stepping_size, start_pos, tool subtype
            dataset.append((ref_subtype,noise,ws,ss,sp,query))
    return dataset

dataset_comet = prepare_dataset("/mnt/data/bioinfo/comet/2014-01-09_noise.fasta.csv", retrieve_query_subtype_comet)

def sens_spec(data, subtypes_to_analyze):
    mean_sens = []
    mean_spec = []
    for subtype in sorted(subtypes_to_analyze):
        tp, fp, tn, fn = 0,0,0,0
        for e in data:
            if e[0] == subtype:
                if subtype in e[5]:
                    tp += 1
                else:
                    fn += 1
            else:
                if subtype not in e[5]:
                    tn += 1
                else:
                    fp += 1
        sensitivity = tp / (tp + fn)
        mean_sens.append(sensitivity)
        specificity = tn / (tn + fp)
        mean_spec.append(specificity)
    return sum(mean_sens) / len(mean_sens), sum(mean_spec) / len(mean_spec)
#        print "subtype {}\t{}\t{}\t{}\t{}\t{}\t{}".format(subtype,tp, fn, tn, fp, sensitivity, specificity)

# look at results per noise level
print "PURE sensitivities"
print "fragment size\tn\t{}".format("\t".join([str(x) for x in sorted(noises)]))
for window in sorted(window_sizes):
    print "{}\t{}\t".format(window,len([x for x in dataset_comet if x[2] == window])),
    for noise in sorted(noises):
        subset =  [x for x in dataset_comet if x[1] == noise and x[2] == window]
        print "{}\t".format(sens_spec(subset, subtypes_pure)[0]),
    print

print "PURE specificities"
print "fragment size\tn\t{}".format("\t".join([str(x) for x in sorted(noises)]))
for window in sorted(window_sizes):
    print "{}\t{}\t".format(window,len([x for x in dataset_comet if x[2] == window])),
    for noise in sorted(noises):
        subset =  [x for x in dataset_comet if x[1] == noise and x[2] == window]
        print "{}\t".format(sens_spec(subset, subtypes_pure)[1]),
    print

print "CRF sensitivities"
print "fragment size\tn\t{}".format("\t".join([str(x) for x in sorted(noises)]))
for window in sorted(window_sizes):
    print "{}\t{}\t".format(window,len([x for x in dataset_comet if x[2] == window])),
    for noise in sorted(noises):
        subset =  [x for x in dataset_comet if x[1] == noise and x[2] == window]
        print "{}\t".format(sens_spec(subset, subtypes_crf)[0]),
    print

print "CRF specificities"
print "fragment size\tn\t{}".format("\t".join([str(x) for x in sorted(noises)]))
for window in sorted(window_sizes):
    print "{}\t{}\t".format(window,len([x for x in dataset_comet if x[2] == window])),
    for noise in sorted(noises):
        subset =  [x for x in dataset_comet if x[1] == noise and x[2] == window]
        print "{}\t".format(sens_spec(subset, subtypes_crf)[1]),
    print
        
        
        
        
        
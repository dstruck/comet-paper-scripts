# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:55:49 2013

@author: Daniel Struck
"""
from __future__ import division

scueal_translation = {
'CRF02':'02_AG','CRF03':'03_AB','CRF04':'04_cpx','CRF05':'05_DF','CRF06':'06_cpx'
,'CRF07':'07_BC','CRF08':'08_BC','CRF09':'09_cpx','CRF10':'10_CD','CRF12':'12_BF'
,'CRF13':'13_cpx','CRF14':'14_BG','CRF15':'15_01B','CRF16':'16_A2D','CRF17':'17_BF'
,'CRF18':'18_cpx','CRF19':'19_cpx','CRF20':'20_BG','CRF21':'21_A2D','CRF22':'22_01A1'
,'CRF23':'23_BG','CRF24':'24_BG','CRF25':'25_cpx','CRF27':'27_cpx','CRF28':'28_BF'
,'CRF29':'29_BF','CRF31':'31_BC','CRF32':'32_06A1','CRF33':'33_01B','CRF34':'34_01B'
,'CRF35':'35_AD','CRF36':'36_cpx','CRF37':'37_cpx','CRF39':'39_BF','CRF40':'40_BF'
,'CRF43':'43_02G','CRF47':'47_BF','AE':'01_AE'}


def retrieve_query_subtype_comet(line):
    token = line.rstrip().split("\t")
    query_temp, id = token[2], token[0]
    query = []
    if "check" in query_temp:
        query.append(query_temp.split(" ")[0])
        query.append(query_temp.split(" ")[3][:-1])
    else:
        query.append(query_temp.upper())
    return query, id

def retrieve_query_subtype_scueal(line):
    token = line.rstrip().split("\t")
    query_temp, id = token[3].replace("-like",""), token[1]
    if query_temp in scueal_translation: # normalize subtype names
        query_temp = scueal_translation[query_temp]
    query = []
    query.append(query_temp)
    return query, id

def retrieve_query_subtype_rega2(line):
    token = line.rstrip().replace('"','').split(",")
    query_temp, id = token[2].upper().split(" "), token[0]
    query = []
    if len(query_temp) > 2:
        query.append(query_temp[2])

    if len(query_temp) > 3:
        query.append(query_temp[3][1:-1])
    return query, id
    
    
noises, window_sizes, stepping_sizes, start_pos, subtypes = set(), set(), set(), set(), set()
data = {}
with open("dataset/2013-08-29_noise_POL_names.csv") as f:
    next(f) #header
    for line in f:
        token = line.rstrip().split("\t")
        data[token[0]] = token[1:]
        subtypes.add(token[1].split(".")[1].upper())
        noises.add(float(token[2]))
        window_sizes.add(float(token[3]))
        stepping_sizes.add(float(token[5]))
        start_pos.add(float(token[5]))

def prepare_dataset(filename_dataset, retrieve_query_subtype):
    dataset = []
    with open(filename_dataset) as f:
        next(f) #header
        for line in f:
            query, id = retrieve_query_subtype(line)
            # ref_subtype, noise, window_size, stepping_size, start_pos, tool subtype
            dataset.append((data[id][0].split(".")[1].upper(), float(data[id][1]),
                            float(data[id][2]),float(data[id][3]),float(data[id][4])
                            , query))
    return dataset

dataset_comet = prepare_dataset("dataset/COMET/2013-08-29_noise_POL.fasta.csv", retrieve_query_subtype_comet)
dataset_scueal = prepare_dataset("dataset/SCUEAL/2013-08-29_noise_POL_SCUEAL.csv", retrieve_query_subtype_scueal)
dataset_rega2 = prepare_dataset("dataset/REGAv2/all.csv", retrieve_query_subtype_rega2)

def sens_spec(data):
    for subtype in sorted(subtypes):
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
        specificity = tn / (tn + fp)
        print "subtype {}\t{}\t{}\t{}\t{}\t{}\t{}".format(subtype,tp, fn, tn, fp, sensitivity, specificity)

# look at results per noise level
for tool in [(dataset_comet,"comet"),(dataset_scueal,"scueal"),(dataset_rega2,"rega2")]:
    print tool[1]
    for noise in sorted(noises):
        for window in sorted(window_sizes, reverse = True):
            print "noise {} window {}\ttp\tfn\ttn\tfp\tsens.\tspec.".format(noise, window)
            subset =  [x for x in tool[0] if x[1] == noise and x[2] == window]
            sens_spec(subset)
            print 
        
        
        
        
        
        
        
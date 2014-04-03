# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 12:47:37 2013

@author: Daniel Struck
"""
from __future__ import division
import psycopg2 

conn = psycopg2.connect(database='comet')
cur = conn.cursor()

pure = ['A1','A2','B','C','D','F1','F2','G','H','J','K']

subtypes = set()
ids = []
with open("dataset/2013-08-30_pol_lanl.csv") as f:
    next(f)
    for line in f:
        token = line.rstrip().split("\t")
        subtype = token[2].upper()
        if "_" in subtype and  int(subtype.split("_")[0]) > 14:
            continue
        subtypes.add(token[2].upper())
        ids.append(token[0])

subtypes = sorted(subtypes)

scueal_translation = {
'CRF02':'02_AG'
,'CRF03':'03_AB'
,'CRF04':'04_cpx'
,'CRF05':'05_DF'
,'CRF06':'06_cpx'
,'CRF07':'07_BC'
,'CRF08':'08_BC'
,'CRF09':'09_cpx'
,'CRF10':'10_CD'
,'CRF11':'11_cpx'
,'CRF12':'12_BF'
,'CRF13':'13_cpx'
,'CRF14':'14_BG'
,'CRF15':'15_01B'
,'CRF16':'16_A2D'
,'CRF17':'17_BF'
,'CRF18':'18_cpx'
,'CRF19':'19_cpx'
,'CRF20':'20_BG'
,'CRF21':'21_A2D'
,'CRF22':'22_01A1'
,'CRF23':'23_BG'
,'CRF24':'24_BG'
,'CRF25':'25_cpx'
,'CRF27':'27_cpx'
,'CRF28':'28_BF'
,'CRF29':'29_BF'
,'CRF31':'31_BC'
,'CRF32':'32_06A1'
,'CRF33':'33_01B'
,'CRF34':'34_01B'
,'CRF35':'35_AD'
,'CRF36':'36_cpx'
,'CRF37':'37_cpx'
,'CRF39':'39_BF'
,'CRF40':'40_BF'
,'CRF43':'43_02G'
,'CRF47':'47_BF'
,'AE':'01_AE'
}


def convert_comet(subtype):
    result = []
    if "check" in subtype:
        result.append(subtype.split(" ")[0].upper())
        result.append(subtype.split(" ")[-1][:-1].upper())
    elif "unassigned" in subtype:
        result.append("RECOMBINANT")
    else:
        result.append(subtype.upper())
    return result

def convert_scueal(subtype):
    subtype = subtype.replace("-like","")
    if subtype in scueal_translation:
        subtype = scueal_translation[subtype]
    result = []
    if "recombinant" in subtype or "Complex" in subtype:
        result.append("RECOMBINANT")
    else:
        result.append(subtype.upper())
    return result

def convert_rega2(subtype):
    result = []
    if "bootscan" in subtype:
        result.append("RECOMBINANT")
    else:
        subtype = subtype.upper().split(" ")
        if len(subtype) > 2:
            result.append(subtype[2])
        if len(subtype) > 3:
            result.append(subtype[3][1:-1])
    return result
    
def convert_blast(subtype):
    result = []
    result.append(subtype.upper())
    return result
##################################################################

gene = "POL"

print "subtype\tn\tTP\tTN\tFP\tFN\tsensitivity\tspecificity"

for subtype in subtypes:
    
    cur.execute("SELECT subtype_lanl, subtype_comet, subtype_rega, subtype_scueal, subtype_blast, subtype_comet_major, name FROM sequences WHERE id IN %s",
                (tuple(ids),))

    count_n = 0
    count_TP, count_FN, count_FP,count_TN = 0, 0, 0, 0

    for row in cur:

        name = row[6]
        subtype_lanl = row[0].upper()
        subtype_comet = convert_comet(row[1])
        subtype_rega = convert_rega2(row[2])
        subtype_scueal = convert_scueal(row[3])
        subtype_blast = convert_blast(row[4])
        subtype_comet_major = convert_comet(row[5])
        subtype_tool = subtype_scueal # select the tool to analyze
        
        if subtype_lanl == subtype:
            count_n += 1
            if subtype in subtype_tool:
                count_TP += 1
            else:
                count_FN += 1
        else:
            if subtype in subtype_tool and subtype_lanl not in subtype_tool:
                count_FP += 1
            else:
                count_TN += 1
    sensitivity, specificity = 0.0, 0.0
    if count_n > 0:
        sensitivity = count_TP / (count_TP + count_FN)            
        specificity = count_TN / (count_TN + count_FP)            
    print "\t".join((subtype, str(count_n), str(count_TP), str(count_TN), str(count_FP), str(count_FN),str(sensitivity),str(specificity)))

    
agree_all = 0
agree_all_pure = 0
agree_all_crf = 0
disagree_all = 0
disagree_all_pure = 0
disagree_all_crf = 0
agree_comet_rega = 0
agree_comet_rega_pure = 0
agree_comet_rega_crf = 0
agree_comet_scueal = 0
agree_comet_scueal_pure = 0
agree_comet_scueal_crf = 0
agree_rega_scueal = 0
agree_rega_scueal_pure = 0
agree_rega_scueal_crf = 0
n = 0

cur.execute("SELECT subtype_lanl, subtype_comet, subtype_rega, subtype_scueal FROM sequences WHERE id IN %s",
            (tuple(ids),))
for row in cur:
    subtype_lanl = row[0].upper()
    subtype_comet = set(convert_comet(row[1]))
    subtype_rega = set(convert_rega2(row[2]))
    subtype_scueal = set(convert_scueal(row[3]))
    
    if subtype_comet.intersection(subtype_rega) and subtype_comet.intersection(subtype_scueal):
        agree_all += 1
        if subtype_lanl in pure: agree_all_pure += 1
        else: agree_all_crf += 1
    elif subtype_comet.intersection(subtype_rega):
        agree_comet_rega += 1
        if subtype_lanl in pure: agree_comet_rega_pure += 1
        else: agree_comet_rega_crf += 1
    elif subtype_comet.intersection(subtype_scueal):
        agree_comet_scueal += 1
        if subtype_lanl in pure: agree_comet_scueal_pure += 1
        else: agree_comet_scueal_crf += 1
    elif subtype_rega.intersection(subtype_scueal):
        agree_rega_scueal += 1
        if subtype_lanl in pure: agree_rega_scueal_pure += 1
        else: agree_rega_scueal_crf += 1
    else:
        disagree_all += 1
        if subtype_lanl in pure: disagree_all_pure += 1
        else: disagree_all_crf += 1
    n += 1

print "n =\t", n
print "all 3 tools agree\t", agree_all
print "all 3 tools agree, PURE\t", agree_all_pure
print "all 3 tools agree, CRF\t", agree_all_crf
print
print "COMET & REGA agree\t", agree_comet_rega
print "COMET & REGA agree, PURE\t", agree_comet_rega_pure
print "COMET & REGA agree, CRF\t", agree_comet_rega_crf
print
print "COMET & SCUEAL agree\t", agree_comet_scueal
print "COMET & SCUEAL agree, PURE\t", agree_comet_scueal_pure
print "COMET & SCUEAL agree, CRF\t", agree_comet_scueal_crf
print
print "REGA & SCUEAL agree\t", agree_rega_scueal
print "REGA & SCUEAL agree, PURE\t", agree_rega_scueal_pure
print "REGA & SCUEAL agree, CRF\t", agree_rega_scueal_crf
print
print "all 3 tools disagree\t", disagree_all
print "all 3 tools disagree, PURE\t", disagree_all_pure
print "all 3 tools disagree, CRF\t", disagree_all_crf

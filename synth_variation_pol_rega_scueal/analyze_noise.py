# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:55:49 2013

@author: daniel
"""
from __future__ import division
from collections import Counter, defaultdict

scueal_translation = {
'CRF02':'02_AG','CRF03':'03_AB','CRF04':'04_cpx','CRF05':'05_DF','CRF06':'06_cpx'
,'CRF07':'07_BC','CRF08':'08_BC','CRF09':'09_cpx','CRF10':'10_CD','CRF11':'11_CPX','CRF12':'12_BF'
,'CRF13':'13_cpx','CRF14':'14_BG','CRF15':'15_01B','CRF16':'16_A2D','CRF17':'17_BF'
,'CRF18':'18_cpx','CRF19':'19_cpx','CRF20':'20_BG','CRF21':'21_A2D','CRF22':'22_01A1'
,'CRF23':'23_BG','CRF24':'24_BG','CRF25':'25_cpx','CRF27':'27_cpx','CRF28':'28_BF'
,'CRF29':'29_BF','CRF31':'31_BC','CRF32':'32_06A1','CRF33':'33_01B','CRF34':'34_01B'
,'CRF35':'35_AD','CRF36':'36_cpx','CRF37':'37_cpx','CRF39':'39_BF','CRF40':'40_BF'
,'CRF43':'43_02G','CRF47':'47_BF','AE':'01_AE'}

''' COMET '''
filename_dataset = "dataset/COMET/2013-08-29_noise_POL.fasta.csv"
def retrieve_query_subtype(line):
    token = line.upper().rstrip().split("\t")
    query_temp, id = token[2], token[0]
    query = []
    if "CHECK" in query_temp:
        query.append(query_temp.split(" ")[0])
        query.append(query_temp.split(" ")[3][:-1])
    else:
        query.append(query_temp)
    return query, id

''' SCUEAL '''
#filename_dataset = "dataset/SCUEAL/2013-08-29_noise_POL_SCUEAL.csv"
#def retrieve_query_subtype(line):
#    token = line.upper().rstrip().split("\t")
#    query_temp, id = token[3].replace("-LIKE",""), token[1]
#    # normalize subtype names
#    if query_temp in scueal_translation:
#        query_temp = scueal_translation[query_temp]
#    query = []
#    query.append(query_temp.upper())
#    return query, id

''' REGAv2 '''
#filename_dataset = "dataset/REGAv2/all.csv"
#def retrieve_query_subtype(line):
#    token = line.upper().rstrip().replace('"','').split(",")
#    query_temp, id = token[2].upper().split(" "), token[0]
#    query = []
#    if len(query_temp) > 2:
#        query.append(query_temp[2])
#
#    if len(query_temp) > 3:
#        query.append(query_temp[3][1:-1])
#    return query, id





pure = ['A1','A2','B','C','D','F1','F2','G','H','J','K']
other = ['O','N','P','CPZ']

noise_levels = set()
fragment_size = set() 
start_positions = defaultdict(set) 

count_fs = Counter()
count_query_fs = Counter()
count_noise_fs = Counter()
count_query_noise_fs = Counter()

count_pure_fs = Counter()
count_pure_query_fs = Counter()
count_pure_noise_fs = Counter()
count_pure_query_noise_fs = Counter()

count_crf_fs = Counter()
count_crf_query_fs = Counter()
count_crf_noise_fs = Counter()
count_crf_query_noise_fs = Counter()

count_noise_fs_sp = Counter()
count_query_noise_fs_sp = Counter()

count_query_noise_subtype = Counter()
count_noise_subtype = Counter()

count_query_fs_noise_subtype = Counter()
count_fs_noise_subtype = Counter()

data = {}
#id      name    noise   window_size     stepping_size   start_pos
with open("dataset/2013-08-29_noise_POL_names.csv") as f:
    next(f) #header
    for line in f:
        token = line.rstrip().split("\t")
        data[token[0]] = token[1:]


with open(filename_dataset) as f:
    next(f) #header
    for line in f:
        query, id = retrieve_query_subtype(line)

        windows_size = int(data[id][2])        
        background = data[id][0].split(".")[1]
        insert = data[id][1].split(".")[1]        
        start_pos_insert = int(data[id][4])

        ref = data[id][0].split(".")[1].upper()
        fs = int(data[id][2])        
        start_pos = int(data[id][4])        
        noise = float(data[id][1])        

        fragment_size.add(fs)
        noise_levels.add(noise)
        start_positions[fs].add(start_pos)
        
        # need to filter out certain CRF for REGAv2 as REGAv2 only supports till 14_BG
        if "_" in ref:
            crf_number = int(ref.split("_")[0])
            if crf_number not in [1,2,3,4,5,6,7,8,9,10,11,12,13,14]:
                continue
        
        if ref in query:
            count_query_fs[fs] += 1
            count_query_noise_fs[noise, fs] += 1
            count_query_noise_fs_sp[noise, fs, start_pos] += 1
            
            if ref in pure:
                count_pure_query_fs[fs] += 1
                count_pure_query_noise_fs[noise, fs] += 1
            else:
                count_crf_query_fs[fs] += 1
                count_crf_query_noise_fs[noise, fs] += 1
                
            count_query_noise_subtype[noise, ref] += 1
            count_query_fs_noise_subtype[fs, noise, ref] += 1

        
        if ref in pure:
            count_pure_fs[fs] += 1
            count_pure_noise_fs[noise, fs] += 1
        else:
            count_crf_fs[fs] += 1
            count_crf_noise_fs[noise, fs] += 1
            
        count_fs[fs] += 1
        count_noise_fs[noise, fs] += 1
        count_noise_subtype[noise, ref] += 1
        count_fs_noise_subtype[fs, noise, ref] += 1
        count_noise_fs_sp[noise, fs, start_pos] += 1


import operator

print "fragment_size\tn\t" + "\t".join([str(x) for x in sorted(noise_levels)])
for fs in sorted(fragment_size):
    print str(fs) + "\t" + str(count_fs[fs]) + "\t" + "\t".join([str(count_query_noise_fs[x, fs]/count_noise_fs[x, fs]) for x in sorted(noise_levels)])
    
print                     
   
print "PURE\tn\t" + "\t".join([str(x) for x in sorted(noise_levels)])
for fs in sorted(fragment_size):
    print str(fs) + "\t" + str(count_pure_fs[fs]) + "\t" + "\t".join([str(count_pure_query_noise_fs[x, fs]/count_pure_noise_fs[x, fs]) for x in sorted(noise_levels)])

print                     

print "CRF\tn\t" + "\t".join([str(x) for x in sorted(noise_levels)])
for fs in sorted(fragment_size):
    print str(fs) + "\t" + str(count_crf_fs[fs]) + "\t" + "\t".join([str(count_crf_query_noise_fs[x, fs]/count_crf_noise_fs[x, fs]) for x in sorted(noise_levels)])

                  
print

print "noise\tsubtype\tn\tquery"
results = []
for e in count_noise_subtype:
    results.append([e[0], e[1], count_noise_subtype[e], count_query_noise_subtype[e] / count_noise_subtype[e]])

results.sort(key = operator.itemgetter(3), reverse = True)
results.sort(key = operator.itemgetter(0))
for row in results:
    print "\t".join([str(x) for x in row])


import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *
pdf = PdfPages('2013-07-23_subtype_noise.pdf')

x = []
y = []
c = []
for noise in sorted(noise_levels):
    for e in count_noise_subtype:
        if e[0] == noise:
            x.append(noise)
            y.append(count_query_noise_subtype[e]/count_noise_subtype[e])
            if e[1] in pure:
                c.append('b')
            elif e[1] in other:
                c.append('g')
            else:
                c.append('r')

plt.scatter(x,y,c=c,s=30.0,alpha=0.5,edgecolors=None)
plt.ylabel('sensitivity')
plt.xlabel('noise level')
plt.title('pure (blue), crf (red) & other (green): sens. / noise')
savefig(pdf, format='pdf') # note the format='pdf' argument!
close()

for fs in sorted(fragment_size, reverse = True):
    x = []
    y = []
    c = []
    for noise in sorted(noise_levels):
        for e in count_fs_noise_subtype:
            if e[1] == noise and e[0]==fs:
                x.append(noise)
                y.append(count_query_fs_noise_subtype[e]/count_fs_noise_subtype[e])
            if e[2] in pure:
                c.append('b')
            elif e[2] in other:
                c.append('g')
            else:
                c.append('r')
    plt.scatter(x,y,c=c,s=30.0,alpha=0.5,edgecolors=None)
    plt.ylabel('sensitivity')
    plt.xlabel('noise level')
    plt.title('fragme_size: ' + str(fs))
    savefig(pdf, format='pdf') # note the format='pdf' argument!
    close()
plt.show()
#pdf.close()


print

for fs in sorted(fragment_size, reverse=True):
    print "fragment size:",fs
    print "n\tnoise\t" + "\t".join([str(x) for x in sorted(start_positions[fs])])
    for noise in sorted(noise_levels):
        print str(count_noise_fs[noise,fs]) + "\t" + str(noise) + "\t" + "\t".join([str(count_query_noise_fs_sp[noise,fs,x]/count_noise_fs_sp[noise,fs,x]) for x in sorted(start_positions[fs])])
    print
    




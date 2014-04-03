# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:55:49 2013

@author: daniel
"""
from __future__ import division
from collections import Counter, defaultdict

''' COMET '''
filename_dataset = "dataset/2014-01-09_noise.fasta.csv"
def retrieve_query_subtype(line):
    token = line.rstrip().split("\t")
    query_temp, _ = token[2], token[0]
    query = []
    if "check" in query_temp:
        query.append(query_temp.split(" ")[0])
        query.append(query_temp.split(" ")[3][:-1])
    else:
        query.append(query_temp)
    return query



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

with open(filename_dataset) as f:
    next(f) #header
    for line in f:
        query = retrieve_query_subtype(line)
        ref_part = line.split("\t")[0]
        param_token = line.split("\t")[0].split("_")

        ref = ref_part.split(".")[1]
        fs = int(param_token[-3][1:])        
        start_pos = int(param_token[-1][1:])        
        noise = float(param_token[-4][1:])        

        fragment_size.add(fs)
        noise_levels.add(noise)
        start_positions[fs].add(start_pos)
        
        if ref not in query:
#            print ref, query
#            print line
#            print
            pass
        else:
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

        count_fs[fs] += 1
        count_noise_fs[noise, fs] += 1
        if ref in pure:
            count_pure_fs[fs] += 1
            count_pure_noise_fs[noise, fs] += 1
        else:
            count_crf_fs[fs] += 1
            count_crf_noise_fs[noise, fs] += 1
            
        count_noise_subtype[noise, ref] += 1
        count_fs_noise_subtype[fs, noise, ref] += 1
        count_noise_fs_sp[noise, fs, start_pos] += 1

#
# FORMAT THE RESULTS
#

#import operator

#print "fragment_size\tn\t" + "\t".join([str(x) for x in sorted(noise_levels)])
#for fs in sorted(fragment_size):
#    print str(fs) + "\t" + str(count_fs[fs]) + "\t" + "\t".join([str(count_query_noise_fs[x, fs]/count_noise_fs[x, fs]) for x in sorted(noise_levels)])
#
#print                     
   
print "PURE\tn\t" + "\t".join([str(x) for x in sorted(noise_levels)])
for fs in sorted(fragment_size):
    print str(fs) + "\t" + str(count_pure_fs[fs]) + "\t" + "\t".join([str(count_pure_query_noise_fs[x, fs]/count_pure_noise_fs[x, fs]) for x in sorted(noise_levels)])

print                     

print "CRF\tn\t" + "\t".join([str(x) for x in sorted(noise_levels)])
for fs in sorted(fragment_size):
    print str(fs) + "\t" + str(count_crf_fs[fs]) + "\t" + "\t".join([str(count_crf_query_noise_fs[x, fs]/count_crf_noise_fs[x, fs]) for x in sorted(noise_levels)])
                  
print

#print "noise\tsubtype\tn\tsensitivity"
#results = []
#for e in count_noise_subtype:
#    results.append([e[0], e[1], count_noise_subtype[e], count_query_noise_subtype[e] / count_noise_subtype[e]])
#
#results.sort(key = operator.itemgetter(3), reverse = True)
#results.sort(key = operator.itemgetter(0))
#for row in results:
#    print "\t".join([str(x) for x in row])


import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *
pdf = PdfPages('2014-01-10_subtype_noise.pdf')

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
plt.scatter(x, y, c=c, s=30.0, alpha=0.5, edgecolors=None)
plt.ylabel('sensitivity')
plt.xlabel('noise level')
plt.title('pure (blue), CRF (red) & other (green)')
savefig(pdf, format='pdf') # note the format='pdf' argument!
close()

for fs in sorted(fragment_size, reverse = True):

    x = []
    y = []
    c = []
    for noise in sorted(noise_levels):
        for e in count_fs_noise_subtype:
            
            if e[1] == noise and e[0] == fs:
                
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
pdf.close()


print

for fs in sorted(fragment_size, reverse=True):
    print "fragment size:{}\t\tstart position".format(fs)
    print "n\tnoise\t" + "\t".join([str(x) for x in sorted(start_positions[fs])])
    for noise in sorted(noise_levels):
        print str(count_noise_fs[noise,fs]) + "\t" + str(noise) + "\t" + "\t".join([str(count_query_noise_fs_sp[noise,fs,x]/count_noise_fs_sp[noise,fs,x]) for x in sorted(start_positions[fs])])
    print
    




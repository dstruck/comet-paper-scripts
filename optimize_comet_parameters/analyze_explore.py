# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 10:58:54 2013

@author: Daniel Struck

Notes:

- COMET output files not included in the script folder (reason: too many files)

- The output of this analysis script are stored under:

explore_window_size_thresholds/output_kmax_7.csv
explore_window_size_thresholds/output_kmax_8.csv

"""
from __future__ import division
import numpy as np
import os

pure = ['A1','A2','B','C','D','F1','F2','G','H','J','K']
directory = "/mnt/data/bioinfo/comet/param/explore/"

print "window_size\tthreshold_pure\tthreshold_crf\tsens. URF\tsens. PURE/CRF"


def retrieve_subtypes(subtype):
    results = []
    if "check" in subtype:
        results.append(subtype.split(" ")[0])
        results.append(subtype.split(" ")[3][:-1])
    else:
        results.append(subtype)
    return results


for fn in os.listdir(directory):
    token = fn[:-4].split("_")
    window_size = token[0]
    threshold_pure = token[1]
    threshold_crf = token[2]
    
    sensitivity_urf = []
    sensitivity_pure_crf = []
    with open(directory + fn) as f:
        next(f)
        for line in f:
            token = line.split("\t")
            reference = token[0].split(",")[1]
            comet_subtype = token[1].rstrip()

            if reference == "URF":            
                if comet_subtype.startswith("unassigned"):
                    sensitivity_urf.append(1)
                else:
                    sensitivity_urf.append(0)
            else:
                reference_subtype = reference.split(".")[1]
                comet_subtypes = retrieve_subtypes(comet_subtype)
                if reference_subtype in comet_subtypes:
                    sensitivity_pure_crf.append(1)
                else:
                    sensitivity_pure_crf.append(0)
                    
    print "{}\t{}\t{}\t{}\t{}".format(window_size, threshold_pure, threshold_crf, \
    np.mean(sensitivity_urf),np.mean(sensitivity_pure_crf))

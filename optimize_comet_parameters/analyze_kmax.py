# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 10:58:54 2013

@author: Daniel Struck

"""
from __future__ import division
import numpy as np

print "kmax\tsensitivity\tmean log likelihood\tmedian log_likelihood\tAIC"

for kmax in range(1,16):
    sensitivity = []
    log_likelihood = []
    with open("kmax/kmax_{}.csv".format(kmax)) as f:
        next(f)
        for line in f:
            token = line.split("\t")
            reference_subtype = token[0].split(",")[1].split(".")[1]
            comet_subtype = token[1]
            log_likelihood.append(float(token[2]))
            if reference_subtype == comet_subtype:
                sensitivity.append(1)
            else:
                sensitivity.append(0)
                
    parameters = 4**kmax * 3
    aic = 2 * parameters - 2  * sum(log_likelihood) # ln(a) = 2.303 * log(a)
                
    print "{}\t{}\t{}\t{}\t{}".format(kmax, np.mean(sensitivity), np.mean(log_likelihood), np.median(log_likelihood), aic)

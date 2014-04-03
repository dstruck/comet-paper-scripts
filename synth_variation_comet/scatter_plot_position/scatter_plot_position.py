#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 11:21:00 2014

@author: Daniel Struck
"""
import matplotlib.pyplot as plt

noise_level = "0.15"

plots, labels = [], []
with open("synthetic_variation_full_length.csv") as f:
    for line in f:
        token = line.rstrip().split(",")
        if "fragment size" in token[0]:
            fs = int(token[0].split(":")[-1])
            labels.append(fs)
        if token[1] == "noise":
            X = [int(x) for x in token[2:] if x]
        elif token[1] == noise_level:
            Y = [float(x) for x in token[2:] if x]

            plot, = plt.plot(X,Y,"o-")
            plots.append(plot)

plt.legend(plots, labels)
plt.title("noise level:{}".format(noise_level))

            
plt.show()

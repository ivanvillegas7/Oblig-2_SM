# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 12:25:09 2022

@author: Iván
"""

import matplotlib.pyplot as plt

r=[]

C=[]

with open("C(r).txt") as infile:
    infile.readline()
    
    lines= infile.readlines()
    
    for line in lines:
        
        vals = line.split()
        
        r.append(float(vals[0]))
        
        vals[1]=vals[1].replace("(","").replace(")","")
        
        C.append(tuple(map(float, vals[1].split(',')))[0])

plt.plot(r, C)
plt.grid(True)
plt.xlabel(r'$r$ []')
plt.ylabel(r'$C(r))$ []')
plt.title(r'Correlation function of a spin chain.')
#plt.savefig('Correlation function.pdf')

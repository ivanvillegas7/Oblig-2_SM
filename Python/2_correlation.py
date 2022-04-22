# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 12:25:09 2022

@author: Iván
"""

import matplotlib.pyplot as plt

r=[]

C1=[]

C2=[]

with open("C(r)_T_0.25.txt") as infile:
    infile.readline()
    
    lines= infile.readlines()
    
    for line in lines:
        
        vals = line.split()
        
        r.append(float(vals[0]))
        
        C1.append(float(vals[1]))
        
with open("C(r)_T_0.5.txt") as infile:
    infile.readline()
    
    lines= infile.readlines()
    
    for line in lines:
        
        vals = line.split()
        
        C2.append(float(vals[1]))

plt.plot(r, C1, label='T/J=0.25')
plt.plot(r, C2, label='T/J=0.50')
plt.legend()
plt.grid(True)
plt.xlabel(r'$r$ []')
plt.ylabel(r'$C(r))$ []')
plt.title(r'Correlation function of a spin chain.')
#plt.savefig('Correlation function.pdf')

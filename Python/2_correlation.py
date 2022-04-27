# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 12:25:09 2022

@author: Iván
"""

import numpy as np

import matplotlib.pyplot as plt

r=[]

C1=[]

c1=[]

C2=[]

c2=[]

l1=np.exp(1/0.25)-1

l2=np.exp(1/0.25)+2

ratio=l2/l1

l_1=np.exp(1/0.5)-1

l_2=np.exp(1/0.5)+2

ratio_=l_2/l_1

N=16

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
        
for i in range(0, len(r)):
    
    c1.append((l2**N*(1/ratio)**r[i]+l1**N*((ratio)**r[i]+1))/(l2**N+2*l1**N))
    
    c2.append((l_2**N*(1/ratio_)**r[i]+l_1**N*((ratio_)**r[i]+1))/(l_2**N+2*l_1**N))

plt.plot(r, C1, label='Simulated T/J=0.25')
plt.plot(r, C2, label='Simulated T/J=0.50')
plt.plot(r, c1, label='Theoretical T/J=0.25')
plt.plot(r, c2, label='Theoretical T/J=0.50')
plt.legend()
plt.grid(True)
plt.xlabel(r'$r$ [position in the chain]')
plt.ylabel(r'$C(r)$ [$1$]')
plt.title(r'Correlation function of a spin chain.')
plt.savefig('Correlation function.pdf')

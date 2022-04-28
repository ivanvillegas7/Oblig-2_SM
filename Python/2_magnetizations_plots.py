# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 18:23:16 2022

@author: Iván Villegas Pérez
"""

import matplotlib.pyplot as plt


T=[]

m_av=[]

m_2=[]

m_2_16=[]

m_4_16=[]

gamma_16=[]

m_2_8=[]

m_4_8=[]

gamma_8=[]

m_2_32=[]

m_4_32=[]

gamma_32=[]

TJ=[]


with open("m_2D.txt") as infile:
    
    infile.readline()
    
    lines= infile.readlines()
    
    
    for line in lines:
        
        vals = line.split()
        
        
        T.append(float(vals[0]))
        
        m_av.append(float(vals[1]))
        
        m_2.append(float(vals[2]))


with open("m_2D_16.txt") as infile:
    
    infile.readline()
    
    lines= infile.readlines()
    
    
    for line in lines:
        
        vals = line.split()
        
        
        TJ.append(float(vals[0]))
        
        m_2_16.append(float(vals[2]))
        
        m_4_16.append(float(vals[3]))


with open("m_2D_8.txt") as infile:
    
    infile.readline()
    
    lines= infile.readlines()
    
    
    for line in lines:
        
        vals = line.split()
        
        
        m_2_8.append(float(vals[2]))
        
        m_4_8.append(float(vals[3]))


with open("m_2D_32.txt") as infile:
    
    infile.readline()
    
    lines= infile.readlines()
    
    
    for line in lines:
        
        vals = line.split()
        
                
        m_2_32.append(float(vals[2]))
        
        m_4_32.append(float(vals[3]))


for i in range(0, 1000):
    
    gamma_16.append((m_4_16[i])/(16*m_2_16[i]**2))

    gamma_8.append((m_4_8[i]**2)/(8*m_2_8[i]**2))
    
    gamma_32.append((m_4_32[i]**2)/(32*m_2_32[i]**2))


plt.figure()
plt.plot(T, m_av)
plt.grid(True)
plt.xlabel(r'$T/J$ [$1$]')
plt.ylabel(r'$\langle m\rangle$ [$1$]')
plt.title(r'Real part of the average magnetization vs temperature.')
plt.savefig('average_magnetization.pdf')


print(r'The average magnetization when T=0 is ')
print(m_av[0])
print(r'The average magnetization when T goes to infinity is arround')
print(m_av[-1])

plt.figure()
plt.plot(T, m_2)
plt.grid(True)
plt.xlabel(r'$T/J$ [$1$]')
plt.ylabel(r'$\langle |m|^2\rangle$ [$1$]')
plt.title(r'Squared average magnetization vs temperature.')
plt.savefig('average_magnetization_squared.pdf')


print(r'The average squared magnetization when T=0 is ')
print(m_2[0])
print(r'The average squared magnetization when T goes to infinity is arround')
print(m_2[-1])


plt.figure()
plt.plot(TJ, gamma_8, label='L=8')
plt.plot(TJ, gamma_16, label='L=16')
plt.plot(TJ, gamma_32, label='L=32')
plt.legend()
plt.grid(True)
plt.xlabel(r'$T/J$ [1]')
plt.ylabel(r'$\Gamma(L=16)$ [$1$]')
plt.title(r'$\Gamma$ vs temperature for different $L$s.')
plt.savefig('gamma_function.pdf')


for i in range(0, len(gamma_16)):
    
    if gamma_16 == gamma_32 == gamma_8:
        
        print(r'$T_c=$')
        print(TJ[i])

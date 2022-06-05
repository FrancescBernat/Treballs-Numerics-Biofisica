# -*- coding: utf-8 -*-
"""
Created on Fri May 27 09:19:53 2022

@author: Francesc Bernat Bieri Tauler
"""

import numpy as np
from Funcions import *
from functools import partial
import matplotlib.pyplot as plt

t0 = 1; tf = 1000

t = np.arange(t0, tf)

E = Euler(2)

dN = lambda t, N, P, b, a = 1: N*(a-b*P)
dP = lambda  t, N, P, c, d = 1: P*(c*N-d)

N = np.zeros(len(t)+1); P = N.copy()  

N[0] = 2;   P[0] = 2;

h = 0.05

for b, c in zip([0.01, 1, 0.75, 1.25], [0.01, 1, 1.25, 0.75]):
    
    dN = partial(dN, b=b)
    dP = partial(dP, c=c)
    
    for i in t:        
        N[i], P[i] = RK2(dN, dP, i, N[i-1], P[i-1], 0.05)
        
        
        # P[i] = E(P[i-1], N[i-1], P[i-1], 0.05, dP)
        # N[i] = E(N[i-1], N[i-1], P[i], 0.05, dN)
        
    fig, ax = plt.subplots(figsize=(8, 6), dpi=400)
    ax.plot(P, label='Depredador', color='seagreen')
    ax.plot(N, label='Pressa', color='skyblue')
    ax.set(ylabel=r"Nº d'individus", xlabel="Temps", title=f"Modelo depredador presa\n\n b = {b}, c = {c}")
    ax.legend()
    
    ax.yaxis.get_ticklocs(minor=True)
    ax.yaxis.set_ticks_position('both'); 
    ax.xaxis.set_ticks_position('both')
    ax.minorticks_on()
    plt.tight_layout()
    
    fig, ax = plt.subplots(figsize=(8, 6), dpi=400)
    ax.plot(N, P, color='darkblue')
    
    ax.set(ylabel=r"N", xlabel="P", 
           title=f"Modelo depredador presa\n \n b = {b}, c = {c}")
    
    ax.yaxis.get_ticklocs(minor=True)
    ax.yaxis.set_ticks_position('both'); 
    ax.xaxis.set_ticks_position('both')
    ax.minorticks_on()
    plt.tight_layout()
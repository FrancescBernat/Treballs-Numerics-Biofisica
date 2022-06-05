# -*- coding: utf-8 -*-
"""
Created on Fri May 27 09:19:53 2022

@author: Francesc Bernat Bieri Tauler
"""

import numpy as np
from Funcions import * 
from functools import partial
import matplotlib.pyplot as plt

# Si volem guardar les gràfiques resultants, posar la variable guardar en true
guardar = False 

# Definim temps inicial, final i pas de temps
t0 = 1; tf = 50; h = 0.005

t = np.arange(t0, tf, h)

# Definim les equacions de Lotka-Volterra
dN = lambda t, N, P, b, a = 1: N*(a-b*P)
dP = lambda  t, N, P, c, d = 1: P*(c*N-d)

N = np.zeros(len(t)); P = N.copy()  

# Condicions inicials
N0 = 2;         P0 = 2;

for b, c in zip([0.01, 1, 0.75, 1.25], [0.01, 1, 1.25, 0.75]):
    
    dN = partial(dN, b=b) 
    dP = partial(dP, c=c)
    
    for i in range(len(t)):
        if i == 0:
        
            N[i] = N0 
            P[i] = P0
            
        else:
            # Usam el runge kutta definit a l'arxiu Funcions.py
            N[i], P[i] = RK2(dN, dP, i, N[i-1], P[i-1], h)
        
    # Grafiques de evolució temporal
    
    fig, ax = plt.subplots(figsize=(8, 6), dpi=400)
    
    ax.plot(P, label='Depredador', color='seagreen')
    ax.plot(N, label='Pressa', color='skyblue')
    
    ax.set(ylabel=r"Nº d'individus", xlabel="Temps", 
           title=f"Modelo depredador presa\n\n b = {b}, c = {c}")
    ax.legend(loc='upper left')
    
    ax.yaxis.get_ticklocs(minor=True)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    
    ax.minorticks_on()
    plt.tight_layout()
    plt.show()
    
    if guardar:
        # Ens asseguram de l'existencià de la subcarpeta
        path = NovaSubcarpeta("Resultats_Treball_3")
        
        # Definim el nom per guardar la gràfica i guardam
        nom = f"ModelDepPres_{b}_{c}".replace('.', '_')
        fig.savefig( path + '\\' + nom + ".pdf", dpi = 400)
        
    
    ###########################################################################
    
        
    # Grafiques on representam la evolució de les poblacions entre si
    fig, ax = plt.subplots(figsize=(8, 6), dpi=400)
    
    ax.plot(N, P, color='darkblue')
    
    ax.set(ylabel=r"N", xlabel="P", 
           title=f"Modelo depredador presa\n \n b = {b}, c = {c}")
    
    ax.yaxis.get_ticklocs(minor=True)
    ax.yaxis.set_ticks_position('both'); 
    ax.xaxis.set_ticks_position('both')
    
    ax.minorticks_on()
    plt.tight_layout()
    plt.show()
    
    if guardar:
        nom += "_PreVsDep"
        fig.savefig( path + '\\' + nom + ".pdf", dpi = 400)
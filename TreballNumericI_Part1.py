# -*- coding: utf-8 -*-
"""
Created on Fri May 27 09:19:53 2022

@author: Francesc Bernat Bieri Tauler
"""

import numpy as np
from Funcions import *
from functools import partial
import matplotlib.pyplot as plt

guardar = False

# Definim el nom de la subcarpeta
NomPath = "Resultats_Treball_1"

# Definim les equacions
dv = lambda t, v, u, I : 0.04*v**2 + 5*v + 140 - u + I
du = lambda t, v, u, a, b : a*(b*v-u)


if __name__ == "__main__": 
    
    # Definim l'interval de temps
    t0 = 0;         tf = 180;       h = 1e-2 # h és el pas d'integració
    t = np.arange(t0, int(tf), h)
    
    I0 = 15
    
    a = 0.02; 
    B = [0.2, 0.2, 0.25]
    C = [-50, -55, -65]
    D = [2, 4, 2]
    
    # Definim les array's on guardarem el resultat
    u = np.zeros(len(t)); v = u.copy()
    
    # Condicions inicials
    u[0] = 2; v[0] = 0
    
    # Assignam un valor de I constant
    dv = partial(dv, I = I0) 
    
    # Preparam un colormap per obtenir gràfiques de colors diferents
    colors = plt.cm.jet(np.linspace(0,1,6))
    
    l = 0 # Per obtenir un color diferent
    
    for b, c, d in zip(B, C, D): # Iteram per cada tipus diferent de neurona
        
        # Assignam valors constants per a la funció per a aquesta integració
        du = partial(du, a=a, b=b) 
        
        for i in range(len(t)):
            
            # Resolem a cada passa pel runge-Kutta definit a Funcions.py
            v[i], u[i] = RK2(dv, du, i, v[i-1], u[i-1], h)
            
            # Condicions de dispar de la neurona
            if v[i] > 30:
                v[i] = c
                u[i] += d
                
        # Representam els resultats
        Grafic(t, v, color=colors[l], 
               title=f"Potencial de membrana \n b={b} c={c} d={d}", 
               path=NomPath, guardar=guardar);   l += 1
        
        Grafic(t, u, color=colors[l],
               title=f"Variable de recuperacio \n b={b} c={c} d={d}", 
               path=NomPath, guardar=guardar);  l += 1
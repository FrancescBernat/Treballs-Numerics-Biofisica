# -*- coding: utf-8 -*-
"""
Created on Fri May 27 09:19:53 2022

@author: Francesc Bernat Bieri Tauler
"""

import numpy as np
from Funcions import *
from functools import partial
import matplotlib.pyplot as plt

NomPath = "Resultats_Treball_1"

dv = lambda t, v, u, I : 0.04*v**2 + 5*v + 140 - u + I
du = lambda t, v, u, a, b : a*(b*v-u)

if __name__ == "__main__":
    
    # Definim l'interval de temps
    t0 = 0;         tf = 200;       h = 1e-2 # h és el pas d'integració
    t = np.arange(t0, int(tf), h)
    
    I0 = 15
    
    a = 0.02; 
    B = [0.2, 0.2, 0.25]
    C = [-50, -55, -65]
    D = [2, 4, 2]
    
    u = np.zeros(len(t)); v = u.copy()
    
    u[0] = 2; v[0] = 0
        
    dv = partial(dv, I = I0)
    
    colors = plt.cm.jet(np.linspace(0,1,6))
    
    l = 0
    for b, c, d in zip(B, C, D):
        
        du = partial(du, a=a, b=b)
        
        for i in range(len(t)):
            
            v[i], u[i] = RK2(dv, du, i, v[i-1], u[i-1], h)
            
            if v[i] > 30:
                v[i] = c
                u[i] = u[i] + d
                
        
        Grafic(t, v, color=colors[l], title=f"Plot \n {b} {c} {d}");   l += 1
        
        Grafic(t, u, color=colors[l], ylabel="I (mA)");  l += 1

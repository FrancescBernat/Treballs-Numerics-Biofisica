# -*- coding: utf-8 -*-
"""
Created on Fri May 27 09:19:53 2022

@author: Francesc Bernat Bieri Tauler
"""

import numpy as np
from numba import jit
from Funcions import *
from functools import partial
import matplotlib.pyplot as plt

# Definim l'interval de temps
t0 = 0;         tf = 30000/1.6
t = np.arange(t0, int(tf))


I0 = 15

a = 0.02; 

B = [0.2, 0.2, 0.25]
C = [-50, -55, -65]
D = [2, 4, 2]

dv = lambda t, v, u, I = I0 : 0.04*v**2 + 5*v + 140 - u + I
du = lambda t, v, u, a, b : a*(b*v-u)

u = np.zeros(len(t)); v = u.copy()

u[0] = 2; v[0] = 0

h = 1e-2
    
colors = plt.cm.jet(np.linspace(0,1,6))

l = 0
for b, c, d in zip(B, C, D):
    
    du = partial(du, a=a, b=b)
    
    for i in t:
        
        v[i], u[i] = RK2(dv, du, i, v[i-1], u[i-1], h)
        
        if v[i] > 30:
            v[i] = c
            u[i] = u[i] + d
            
    
    Grafic(t, v, color=colors[l], title=f"Plot \n {b} {c} {d}");   l += 1
    
    Grafic(t, u, color=colors[l], ylabel="I (mA)");  l += 1

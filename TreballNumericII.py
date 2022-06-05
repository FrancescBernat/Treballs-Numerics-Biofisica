# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 12:41:13 2022

@author: Francesc Bernat Bieri Tauler
"""

import numpy as np
from tqdm import tqdm
from numba import jit
from Funcions import *
from functools import partial
import matplotlib.pyplot as plt

# Definim l'interval de temps
t0 = 1;         tf = 250;         h = 1e-2
t = np.arange(t0, int(tf), h)

# Definim les condicions inicials

I0 = 15

a = 0.02;   b = 0.2;    c = -65;    d = 8

dr = lambda i, r, tau = 5.6 : -r/tau
Ic = lambda r, v, E, gmax = 0.1 :  gmax * r * (v - E) # Intensitat conectora

dv = lambda t, v, u, Ic, I = I0 : 0.04*v**2 + 5*v + 140 - u + I + Ic
du = lambda t, v, u, a, b : a*(b*v-u)

# Assignam les array on guardarem les dades
u1 = np.zeros(len(t)+1);
v1, u2, v2, r1, r2 = map(np.copy, [u1]*5)

u1[0] = 0; v1[0] = 0

du = partial(du, a=a, b=b)

for i in tqdm(range(len(t))):
    
    r1[i] = RK(dr, i, r1[i-1], h, r1[i-1])
    r2[i] = RK(dr, i, r2[i-1], h, r2[i-1])
    
    
    dv1 = partial(dv, Ic = Ic(r2[i-1], v2[i-1], 0))
    dv2 = partial(dv, Ic = Ic(r1[i-1], v1[i-1], -53))
    
    v1[i], u1[i] = RK2(dv1, du, i, v1[i-1], u1[i-1], h)
    v2[i], u2[i] = RK2(dv2, du, i, v2[i-1], u2[i-1], h)
    
    if v1[i] >= 30:
        v1[i] = c
        u1[i] += d
        
    if v2[i] >= 30:
        v2[i] = c
        u2[i] += d
        
        
plt.figure(dpi=300)
plt.plot(v1, '--', label='Neurona 1')
plt.plot(v2, '-.', label='Neurona 2')
plt.show()
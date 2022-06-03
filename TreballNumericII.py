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
t0 = 1;         tf = 1000;         h = 1e-2
t = np.arange(t0, int(tf), h)

# Definim les condicions inicials

I0 = 15

a = 0.02;   b = 0.2;    c = -65;    d = 8

dv = lambda t, v, u, I = I0 : 0.04*v**2 + 5*v + 140 - u + I
du = lambda t, v, u, a, b : a*(b*v-u)

u = np.zeros(len(t)+1); v = u.copy()

u[0] = 2; v[0] = 0

du = partial(du, a=a, b=b)

for i in tqdm(range(len(t))):
    
    v[i], u[i] = RK2(dv, du, i, v[i-1], u[i-1], h)
    
    if v[i] > 30:
        v[i] = c
        u[i] = u[i] + d
        
plt.figure(dpi=300)
plt.plot(v)
plt.show()
# -*- coding: utf-8 -*-
"""
Created on Fri May 27 09:19:53 2022

@author: Francesc Bernat Bieri Tauler
"""

import numpy as np
from Funcions import *
from tqdm import tqdm
from multiprocessing import Process
import scipy.signal as sgn
from functools import partial
import matplotlib.pyplot as plt

from TreballNumericI_Part1 import dv, du, NomPath

# Definim l'interval de temps
t0, tf = 1, 1000
h = 5e-3 # tf = 1000 i h = 1e-2 donen un resultat que esta be

t = np.arange(t0, int(tf), h)


Imin, Imax = 15, 30
I = np.linspace(Imin, Imax, 150)

a, b, c, d = 0.02, -0.1, -55, 6

# dv = lambda t, v, u, I: 0.04*v**2 + 5*v + 140 - u + I
# du = lambda t, v, u, a=a, b=b : a*(b*v-u)
du = partial(du, a=a, b=b)

u = np.zeros(len(t)+1); v = u.copy()
Freq = np.zeros(len(I))

u[0] = 0;           v[0] = -60

Euler = Euler()

n = 0
find_peaks = sgn.find_peaks

import time; t0 = time.time()

def Time(u, v, t):
    for i in range(len(t)):
        
        u[i] = Euler(u[i-1], v[i-1], u[i-1], h, du)
        v[i] = Euler(v[i-1], v[i-1], u[i-1], h, dv)
        
        # v[i], u[i] = RK2(dv, du, i, v[i-1], u[i-1], h)
        
        if v[i] >= 30:
            v[i] = c
            u[i] += d
    return u, v

for I0 in tqdm(I):
    
    dv = partial(dv, I=I0)
    
    u, v = Time(u, v, t)
            
    peaks, _ = find_peaks(v[2500:], height=15)
    # peaks = sgn.argrelextrema(v[2500:], np.greater)[0]
    if len(peaks) > 2:
        Freq[n] = 1 /  ( np.unique(np.diff(peaks))[0] )
    n += 1
    
print(f"\n Han passat {time.time()-t0} segons!")
peaks, _ = find_peaks(v, height=20)

plt.plot(v)
plt.plot(peaks, v[peaks], "x")

plt.show()

Grafic(I, Freq)
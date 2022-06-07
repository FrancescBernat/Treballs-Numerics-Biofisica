# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 12:41:13 2022

@author: Francesc Bernat Bieri Tauler
"""

import numpy as np
from tqdm import tqdm
from Funcions import *
from functools import partial
import matplotlib.pyplot as plt

# Definim l'interval de temps
t0 = 0;         tf = 200;         h = 5e-3
t = np.arange(t0, int(tf), h)

# Definim les condicions inicials

I0 = 15

a = 0.02;   b = 0.2;    c = -65;    d = 8

# dr = lambda i, r, tau = 5.6 : -r/tau

# Definim les funcions que anirem a usar
def dr(i, r, v, tau=5.6):
    if v > 30:  return -r/tau + 0.5 
    
    else:   return -r/tau

Ic = lambda r, v, E, gmax = 0.1 :  gmax * r * (v - E) # Intensitat conectora

dv = lambda t, v, u, Ic, I = I0 : 0.04*v**2 + 5*v + 140 - u + I + Ic
du = lambda t, v, u, a, b : a*(b*v-u)

# Assignam les array on guardarem les dades
u1 = np.zeros(len(t));
v1, u2, v2, r12, r21 = map(np.copy, [u1]*5)

# Condicions inicials dels resultats
u1[0] = 0.1*-20; v1[0] = -70
u2[0] = 0.1*-10; v2[0] = -70

du = partial(du, a=a, b=b)

for i in tqdm(range(len(t))):
    
    if i == 0: # En cas de estar a l'inici, no fer res
        continue
    
    else:
        
        # Condició de dispar neurona 1
        if v1[i-1] > 30:
            dr1 = partial(dr, v = v1[i-1])
            v1[i-1] = c
            u1[i-1] += d
            
        else:
            dr1 = partial(dr, v = v1[i-1])
            
        # Condició de dispar neeurona 2 
        if v2[i-1] > 30:
            dr2 = partial(dr, v = v2[i-1])
            v2[i-1] = c
            u2[i-1] += d
            
        else:
            dr2 = partial(dr, v = v2[i-1])
        
        r12[i] = RK(dr1, i, r12[i-1], h, r12[i-1])
        r21[i] = RK(dr2, i, r21[i-1], h, r21[i-1])
        
        
        dv1 = partial(dv, Ic = 0)#Ic(r21[i-1], v1[i-1], -53))
        dv2 = partial(dv, Ic = Ic(r12[i-1], v2[i-1], -53))
        
        v1[i], u1[i] = RK2(dv1, du, i, v1[i-1], u1[i-1], h)
        v2[i], u2[i] = RK2(dv2, du, i, v2[i-1], u2[i-1], h)
    
    
vt = abs(v1 - v2); vt = vt[vt != 0];  
if sum(vt) == 0: vt = 0
print(f"\n La diferencia mitjana entre el potencial de les dues " \
      f"neurones és de {round(np.mean(vt), 6)}")
        
plt.figure(dpi=300)
plt.plot(v2, '-', label='Neurona 1')
plt.plot(v1, '-', label='Neurona 2')
# plt.plot(u1, '-', label='Neurona 1')
# plt.plot(u2, '--', label='Neurona 2')
# plt.legend()
plt.show()
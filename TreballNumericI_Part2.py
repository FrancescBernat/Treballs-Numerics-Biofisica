# -*- coding: utf-8 -*-
"""
Created on Fri May 27 09:19:53 2022

@author: Francesc Bernat Bieri Tauler
"""

import numpy as np
from tqdm import tqdm
from time import time
import scipy.signal as sgn
from functools import partial
from Funcions import Euler, Grafic

# Importam de la primera part les funcions i el nom de la subcarpeta
from TreballNumericI_Part1 import dv, du, NomPath
   
def main(du = du, dv = dv, NomPath = NomPath, guardarGrafics = False, 
         Euler = Euler(1)):

    t0 = time() # Per després saber el temps que ha tardat el programa
    
    # Extreim les següents funcions dels mòduls
    find_peaks = sgn.find_peaks
    
    # Interval de temps
    ti, tf = 0, 800
    h = 5e-3 # h és el pas d'integració
    t = np.arange(ti, tf, h)
    
    # Arrays on guardam els resultats
    u = np.zeros(len(t)); v = u.copy()
    
    # Les 2 variacions de corrent del problema
    I1 = np.linspace(15, 30, 150)
    I2 = np.linspace(-4, 4, 80)
    
    # Si definim la solució de les equacions com una funció aconseguim que el 
    # programa duri uns segons menys
    def ResolEqu(u, v, t):
        for i, t in enumerate(t):
            
            u[i] = Euler(u[i-1], v[i-1], u[i-1], h, du)
            v[i] = Euler(v[i-1], v[i-1], u[i-1], h, dv)
                    
            if v[i] >= 30:
                v[i] = c
                u[i] += d
                
        return u, v
    
    # Iteram per als 2 tipus de neurones
    for a, b, c, d, v0, I, tipus in zip( [0.02, 0.2], [-0.1, 0.26], [-55, -65],
                                 [6, 0], [-60, 0], [I1, I2], ["I", "II"]):
        
        # Actualitzam la equació amb els parametres actuals
        du = partial(du, a = a, b = b)
        
        # Cream array per a guardar els resultats
        Freq = np.zeros(len(I))
        
        # Definim les condicions inicials
        u[0] = 0
        v[0] = v0     
        
        
        n = 0
        for I0 in tqdm(I): # Iteram per cada una de les Intensitats
            
            dv = partial(dv, I = I0) # Actualitzam la funció amb la intensitat
            
            u, v = ResolEqu(u, v, t)
                    
            # Cercam els pics a partir de quan estan estables
            peaks, _ = find_peaks(v[2500:])
            
            # Si hi ha més de 2 pics, cercam la freqüencia
            if len(peaks) > 2:
                # Agafam la mitjana de les distancies
                Freq[n] = 1 /  np.mean(np.diff(peaks)) 
                
            n += 1
            
        # Representam la freqüencia en funció de la corrent
        Grafic(I, Freq, ylabel = r"$\nu$ (KHz)", xlabel = r"I (mV)",
               title = f"Excitable clase {tipus} \n", 
               path = NomPath, guardar = guardarGrafics)
    
    print(f"\n Han passat {round( (time()-t0)/60 , 2)} minuts")
    
if __name__ == "__main__": 
    main()
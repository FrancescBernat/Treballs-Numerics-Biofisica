# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 12:41:13 2022

@author: Francesc Bernat Bieri Tauler
"""

import numpy as np
from Funcions import RK, RK2, VarisGrafics
from functools import partial


def main(Grafics=True, guardar=False, NomPath="Resultats_Treball_2"):
    
    t0, tf = 0, 150
    h = 5e-3 # h és el pas d'integració
    
    # Interval de temps
    t = np.arange(t0, int(tf), h)
    
    # Condicions inicials
    I0 = 15
    a = 0.02
    
    # Funcions 
    def dr(i, r, v, tau=5.6):
        if v > 30:  
            return -r/tau + 0.5 
        else:   
            return -r/tau
    
    def Ic(r, v, E, gmax = 0.1): return  gmax * r * (v - E) # Intensitat conectora
    
    def dv(t, v, u, Ic, I = I0): return 0.04*v**2 + 5*v + 140 - u + I + Ic
    
    def du(t, v, u, a, b): return a*(b*v - u)
    
    
    # Assignam les array on guardarem les dades
    u1 = np.zeros(len(t));
    v1, u2, v2, r12, r21 = map(np.copy, [u1]*5)
    
    # Condicions inicials dels resultats
    v1[0] = v2[0] = - 70
    
    apartat = 0 # Per guardar el nom de l'apartat on esteim
    
    # Iteram per a les 2 variacions de dades
    for b, c, d in zip([0.2,0.25], [-65,-50], [8,2]):
    
        apartat += 1
        
        du = partial(du, a=a, b=b)
        
        # Iteram per 2 variacions de les condicions inicials
        for u10, u20 in zip([-2, -2.5], [-1, -2.5]):
            
            u1[0] = u10
            u2[0] = u20
            
            nn = -1
            for I12, I21, titol in zip([None, None, 0, -65], [0, -65, 0, -65], 
                                       ["Unidireccional excitatorio", 
                                        "Unidireccional inhibitorio",
                                        "Bidireccional excitatorio",
                                        "Bidireccional inhibitorio"]):
                nn += 1
                
                for i in (range( len(t) )):
                    
                    if i == 0: # En cas de estar a l'inici, no fer res
                        continue
                    
                    else:
                        
                        ####################################
                        ### Condició de dispar neurona 1 ### <-- Aixì m'assegur que s'ha tengut en compte a cada passa
                        ####################################
                        if v1[i-1] >= 30:
                            dr1 = partial(dr, v=v1[i-1])
                            v1[i-1] = c
                            u1[i-1] += d
                            
                        else:
                            dr1 = partial(dr, v=v1[i-1])
                            
                        ####################################    
                        ### Condició de dispar neurona 2 ###
                        ####################################
                        if v2[i-1] >= 30:
                            dr2 = partial(dr, v=v2[i-1])
                            v2[i-1] = c
                            u2[i-1] += d
                            
                        else:
                            dr2 = partial(dr, v=v2[i-1])
                        
                        r12[i] = RK(dr1, i, r12[i-1], h, r12[i-1])
                        r21[i] = RK(dr2, i, r21[i-1], h, r21[i-1])
                        
                        # Definim el tipus de conexió entre neurones 
                        if not I12:
                            dv1 = partial(dv, Ic=0)
                        else:
                            dv1 = partial(dv, Ic=Ic(r21[i-1], v1[i-1], I12))
                        
                        dv2 = partial(dv, Ic=Ic(r12[i-1], v2[i-1], I21))
                        
                        
                        v2[i], u2[i] = RK2(dv2, du, i, v2[i-1], u2[i-1], h)
                        v1[i], u1[i] = RK2(dv1, du, i, v1[i-1], u1[i-1], h)
                    
                    
                # Trobam la diferencia entre els 2 potencials
                vt = abs(v1 - v2)
                vt = vt[vt != 0] 
                
                if sum(vt) == 0: # Per no tenir problemes amb nan
                    vt = 0
                
                print(f"\n La diferencia mitjana entre el potencial de les dues " \
                      f"neurones, pel cas {titol} i condicions de l'apartat {apartat}, "\
                      f"és de {round(np.mean(vt), 6)} \n")
                        
               
                if Grafics:
                    
                    # Per saber quina neurona és quina
                    neurona1 = ["excitadora", "inhibitora", 
                                "excitadora", "inhibitora"]
                    neurona2 = [None, None, "excitadora", "inhibitora"]
                    
                    # Donam un titol per a cada gràfic
                    titul_v = f"Potencial de membrana apartat" \
                        f" {apartat} \n" + titol + f"u10={u10} u20={u20}"
                        
                    titul_u = f"Variable de recuperació apartat " \
                        f" {apartat} \n" + titol + f"u10={u10} u20={u20}"
                    
                    # Cambiam el tipus de linia, per a poder veure en tot moment
                    # els potencials de les dues neurones
                    if u10 == u20:  
                        line2 = "--"
                    else:   
                        line2 = "-"
    
                    # Representam els resultats                    
                    VarisGrafics(t, v1, v2, title=titul_v, line2=line2, 
                                 color1="midnightblue", color2="aquamarine", 
                                 path=NomPath, guardar=guardar, 
                                 label1=neurona1[nn], label2=neurona2[nn])
                    
                    VarisGrafics(t, u1, u2, title=titul_u, line2=line2, 
                                 ylabel="u (mV)", color2="aquamarine", 
                                 path=NomPath, guardar=guardar, 
                                 label1=neurona1[nn], 
                                 label2=neurona2[nn])


if __name__ == "__main__": 
    main()
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 15:37:03 2022

@author: Francesc Bernat Bieri Tauler
"""

def Euler(*arg):
    '''
    Funció que retorna el mètode d'Euler i totes les variacions del mètode que
    varem veure a Física Computacional. S'ha usat la mateixa nomenclatura que
    se va usar allà.

    Parameters
    ----------
    *arg : 1 o "EulerMillorat"  --> Per obtenir la versió que varem anomenar
                                    Euler Millorat
           2 o "EulerModificat" --> Per obtenir la versió que varem anomenar
                                    Euler Modificat
           Cualsevol altra cosa --> Obtendrem Euler Simple
           
    Returns
    -------
    FUNCTION
        Mètode d'Euler amb la variació especificada, per defecte retorna 
        l'Euler simple.

    '''
    # Per a l'Euler Millorat
    if arg == 'EulerMillorat' or 1:
        return lambda CI, x, y, h, f: CI + (h/2)*(f(h,x,y) + f(h, x + h, y + h*f(h,x,y))) 
        
    # Per a l'Euler Modificat
    elif arg == 'EulerModificat' or 2:
        return lambda CI, x, y, h, f: CI + h*f(h, x + h/2, y + (h/2)*f(h, x,y) )
    
    # Per a l'Euler Simple
    else:
        return lambda CI, x, y, h, f: CI + h*f(h, x, y)
        
    

def RK(f, x0, y0, h, CI): 
    '''    
    Runge-Kutta de 5 orden para una funcion.
    Devuelve la solución  para el instante t.

    Función creada originalmente para la pràctica 1 de 
    Física Computacional.


    Parameters
    ----------
    f : FUNCTION
        La función que queremos resolver.
    x0 : FLOAT
        Primer variable de la función.
    y0 : FLOAT
        Segunda variable de la función.
    h : FLOAT
        Paso de integración.
    CI : FLOAT
        Estado inicial de la función para este paso.

    Returns
    -------
    SolO5 : FLOAT
        Solución de 5 orden para el instante t.

    '''
    
    f0 = f(x0, y0)
    
    f1 = f(x0 + (h/4), 
           y0 + (h/4)*f0)
    
    f2 = f(x0 + ( (3*h) / 8 ), 
           y0 + ( (3*h) / 32 )*f0 + ( (9*h) / 32)*f1)
    
    f3 = f(x0 + ( (12*h) / 13), 
           y0 + ( (1932*h) / 2197 )*f0 - ( (7200*h) / 2197 )*f1 + ( (7296*h) / 2197 )*f2 )
        
    f4 = f(x0 + h, 
           y0 + ( (439*h) / 216 )*f0 - (8*h)*f1 + ( (3680*h) / 513 )*f2 - ( (845*h) / 4104)*f3)
    
    f5 = f(x0 + (h/2), 
           y0 - ( (8*h) / 27 )*f0 + (2*h)*f1 - ( (3544*h) / 2565 )*f2  \
           + ( (1859*h) / 4104)*f3 - ((11*h)/40)*f4)
        
    
    SolO5 = CI + h*( (16 / 135)*f0 + (6656 / 12825)*f2 + (28561 / 56430)*f3 - (9 / 50)*f4 + (2 / 55)*f5) # Solució de 5 ordre
        
    return SolO5 

def RK2(f, g, t0, x0, v0, h):
    '''
    Runge-Kutta de 5 orden para 2 funciones.
    Devuelve la solución solamente para el instante t.
    
    Función creada originalmente para la pràctica 1 de 
    Física Computacional.

    Parameters
    ----------
    f : FUNCTION
        Una de las funciones que queremos resolver.
    g : FUNCTION
        La otra funcion que queremos resolver.
    t0 : FLOAT
        Primer parametro de la funcion para este paso.
    x0 : FLOAT
        Segundo parametro de la funcion para este paso.
    v0 : FLOAT
        Tercer parametro de la funcion para este paso.
    h : FLOAT
        Paso de integración.

    Returns
    -------
    SolF : FLOAT
        Solucion de la función F, para el paso de tiempo t, 
        aplicando un Runge Kutta de 5 orden.
    SolG : FLOAT
        Solucion de la función G, para el paso de tiempo t, 
        aplicando un Runge Kutta de 5 orden.

    '''
    
    f0 = f(t0, x0, v0)
    g0 = g(t0, x0, v0)
    
    f1 = f(t0 + (h/4),  x0 + (h/4)*f0, v0 + (h/4)*g0)
    g1 = g(t0 + (h/4),  x0 + (h/4)*f0, v0 + (h/4)*g0)
    
    f2 = f(t0 + ( (3*h) / 8 ), 
           x0 + ( (3*h) / 32 )*f0 + ( (9*h) / 32)*f1, 
           v0 + ( (3*h) / 32 )*g0 + ( (9*h) / 32)*g1)
    g2 = g(t0 + ( (3*h) / 8 ), 
           x0 + ( (3*h) / 32 )*f0 + ( (9*h) / 32)*f1, 
           v0 + ( (3*h) / 32 )*g0 + ( (9*h) / 32)*g1)
    
    f3 = f(t0 + ( (12*h) / 13), 
           x0 + ( (1932*h) / 2197 )*f0 - ( (7200*h) / 2197 )*f1 + ( (7296*h) / 2197 )*f2,
           v0 + ( (1932*h) / 2197 )*g0 - ( (7200*h) / 2197 )*g1 + ( (7296*h) / 2197 )*g2)
    g3 = g(t0 + ( (12*h) / 13), 
           x0 + ( (1932*h) / 2197 )*f0 - ( (7200*h) / 2197 )*f1 + ( (7296*h) / 2197 )*f2,
           v0 + ( (1932*h) / 2197 )*g0 - ( (7200*h) / 2197 )*g1 + ( (7296*h) / 2197 )*g2)
        
    f4 = f(t0 + h, 
           x0 + ( (439*h) / 216 )*f0 - (8*h)*f1 + ( (3680*h) / 513 )*f2 - ( (845*h) / 4104)*f3,
           v0 + ( (439*h) / 216 )*g0 - (8*h)*g1 + ( (3680*h) / 513 )*g2 - ( (845*h) / 4104)*g3)
    g4 = g(t0 + h, 
           x0 + ( (439*h) / 216 )*f0 - (8*h)*f1 + ( (3680*h) / 513 )*f2 - ( (845*h) / 4104)*f3,
           v0 + ( (439*h) / 216 )*g0 - (8*h)*g1 + ( (3680*h) / 513 )*g2 - ( (845*h) / 4104)*g3)
    
    
    f5 = f(t0 + (h/2), 
           x0 - ( (8*h) / 27 )*f0 + (2*h)*f1 - ( (3544*h) / 2565 )*f2  \
           + ( (1859*h) / 4104)*f3 - ((11*h)/40)*f4,
           v0 - ( (8*h) / 27 )*g0 + (2*h)*g1 - ( (3544*h) / 2565 )*g2  \
           + ( (1859*h) / 4104)*g3 - ((11*h)/40)*g4)
    g5 = g(t0 + (h/2), 
           x0 - ( (8*h) / 27 )*f0 + (2*h)*f1 - ( (3544*h) / 2565 )*f2  \
           + ( (1859*h) / 4104)*f3 - ((11*h)/40)*f4,
           v0 - ( (8*h) / 27 )*g0 + (2*h)*g1 - ( (3544*h) / 2565 )*g2  \
           + ( (1859*h) / 4104)*g3 - ((11*h)/40)*g4)
        
    SolF = x0 + h*( (16 / 135)*f0 + (6656 / 12825)*f2 + (28561 / 56430)*f3 - (9 / 50)*f4 + (2 / 55)*f5)
    
    SolG = v0 + h*( (16 / 135)*g0 + (6656 / 12825)*g2 + (28561 / 56430)*g3 - (9 / 50)*g4 + (2 / 55)*g5)
    
    return SolF, SolG

def NovaSubcarpeta(nomSubcarpeta = 'Subfolder'):
    '''
    Aquesta funció generarà una subcarpeta del nom especificat
    en cas de que no existeixi.
    
    Retorna la localització d'aquesta subcarpeta.
    '''
    
    import os
    
    # Per a trobar el dirrectori complet d'aquest arxiu python
    dir_path = os.path.dirname(os.path.realpath(__file__))
    
    # Per a crear una subcarpeta en aquest nou directori
    os.makedirs(dir_path + '/' + nomSubcarpeta, exist_ok=True) 
    
    # Retornam el directori complet de la nova subcarpeta
    return os.path.realpath(nomSubcarpeta)
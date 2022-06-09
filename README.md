# Treballs-Numerics-Biofisica

Resolució dels treballs numérics de l'assignatura de Biofísica d'un grau de Física

- ***Funcions***

> Arxiu que conte funcions necessaries per fer funcionar la resta d'aquests arxius.
>> Depen dels mòduls *matplotlib*, *re* i *os*

- ***TreballNumericI_part1*** 

> Simulació d'una única neurona (potencial de membrana i variable de recuperació) usant el model de Izhikevich.
>> Depen dels mòduls *matpltolib*, *numpy*, i *functools* 

- ***TreballNumericI_part2***
 
> Representació de la freqüencia de dispar d'una sola neurona frent a un corrent injectat, usant el model de Izhikevich.
>> Depen dels mòduls *numpy*, *scipy.signal*, *functools*, *tqdm* i *time* 

- ***TreballNumericII*** 

> Estudi de l'acoplament de dues neurones usant el model de Izhikevich.
>> Depen dels mòduls *numpy* i *functools* 

- ***TreballNumericIII*** 

> Integració de les equacions de Lotka-Volterra, per uns quants casos.
>> Depen dels mòduls *matpltolib*, *numpy*, i *functools* 


A part dels mòduls mencionats abans, els arxius de la "resolució" de cada una d'aquestes pràctiques depenen de l'arxiu *Funcions.py*, 
on he inclós funcions importants com el mètode de Runge-Kutta.
L'arxiu *TreballNumericI_part2* també depen de *TreballNumericI_part1* per 3 coses, les dues equacions del model i el nom de la carpeta on és guardaran les gràfiques de així
especificar-ho.

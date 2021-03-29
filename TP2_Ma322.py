# -*- coding: utf-8 -*-
"""
Nom : EL KHALIL
Prénom : Ali
Module : Ma322
Titre : TP2_Ma322
"""

import numpy as np
import matplotlib.pyplot as plt



#------------------------------------------
#Oscillations du pendule
#------------------------------------------
#Déclaration des constantes du problème.
g=9.81
"""
L=1

w0=(g/L)**0.5
f=2*np.pi/w0
T=1/f

h=0.04
A=np.pi/2
phi=0
"""
def Pendule_linearise(w0, A, phi, h, Ta):        #Ta est le temps d'acquisition
    N=int(Ta/h)
    t=np.linspace(0,4,N)
    teta=A*np.cos(w0*t+phi)
    plt.plot(t,teta)
    plt.title("Tracé de teta(t) linéarisé")
    plt.xlabel("Temps (s)")
    plt.ylabel("Angle teta (rad)")
    plt.show()
    

    




# -*- coding: utf-8 -*-
"""
Nom : EL KHALIL
Prénom : Ali
Module : Ma322
Titre : TP2_Ma322
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as si


#------------------------------------------
#Oscillations du pendule
#------------------------------------------
#Déclaration des constantes du problème.
g=9.81
L=1

w0=(g/L)**0.5
f=2*np.pi/w0
T=1/f

h=0.04
A=np.pi/2
phi=0

Y_0=np.array([0, np.pi/2])

def Pendule_linearise(h, Ta, figure=False):
    N=int(Ta/h)
    t=np.linspace(0,Ta,N+1)
    teta=A*np.cos(w0*t+phi)
    if (figure):
        plt.plot(t,teta)
        plt.title("Tracé de teta(t) linéarisé")
        plt.xlabel("Temps (s)")
        plt.ylabel("Angle teta (rad)")
        plt.show()
    return (teta)
    

    
def pendule(Y,t):
    return np.array([-(w0**2)*np.sin(Y[1]), Y[0]])

def Euler(h, Ta, figure=False):
    N=int (Ta/h)
    Y_e=np.zeros((N+1,2))
    Y_e[0]=Y_0
    t=[0]
    for i in range(1,N+1):
        Y_e[i]=Y_e[i-1]+h*pendule(Y_e[i-1], t[i-1])
        t.append(h*i)
    if (figure):
        plt.plot(t, Y_e[:,1])
        plt.title("Tracé de teta(t) Euler")
        plt.xlabel("Temps (s)")
        plt.ylabel("Angle teta (rad)")
        plt.show()
    return (Y_e)


def Runge_Kutta_4(h,Ta,figure=False):
    N=int (Ta/h)
    Y_rk=np.zeros((N+1,2))
    Y_rk[0]=Y_0
    t=[0]
    for i in range(1, N+1):
        k1=pendule(Y_rk[i-1], t[i-1])
        k2=pendule(Y_rk[i-1]+(h/2)*k1,t[i-1]+h/2)
        k3=pendule(Y_rk[i-1]+(h/2)*k2,t[i-1]+h/2)
        k4=pendule(Y_rk[i-1]+h*k3,t[i-1]+h)
        Y_rk[i]=Y_rk[i-1]+(h/6)*(k1+2*k2+2*k3+k4)
        t.append(h*i)
    if (figure):
        plt.plot(t, Y_rk[:,1])
        plt.title("Tracé de teta(t) Runge-Kutta")
        plt.xlabel("Temps (s)")
        plt.ylabel("Angle teta (rad)")
        plt.show()
    return (Y_rk)
    
    
def figures_teta(h,Ta):
    teta_e=Euler(h,Ta)[:,1]
    teta_l=Pendule_linearise(h, Ta)
    teta_rk=Runge_Kutta_4(h, Ta)[:,1]
    N=int(Ta/h)
    t=np.linspace(0,Ta,N+1)
    Yode =si.odeint(pendule, Y_0, t)
    plt.plot(t,teta_e,label="Teta Euler")
    plt.plot(t, teta_l, label="Teta linéarisé")
    plt.xlabel("Temps (s)")
    plt.ylabel("Angle teta (rad)")
    plt.plot(t, teta_rk, label="Teta Runge-Kutta")
    plt.plot(t, Yode[:,1],"+",label="Odeint")
    plt.title("Différent tracé de teta")
    plt.legend()
    plt.show()
    
def portrait_de_phase(h, Ta):
    Y_e=Euler(h,Ta)
    Y_rk=Runge_Kutta_4(h, Ta)
    N=int(Ta/h)
    t=np.linspace(0,Ta,N+1)
    Yode =si.odeint(pendule, Y_0, t)
    plt.plot(Y_e[:,1], Y_e[:,0])
    plt.title("Portrait de phase avec la méthode d'Euler")
    plt.xlabel("Teta")
    plt.ylabel("Teta'")
    plt.axis("equal")
    plt.show()
    plt.plot(Y_rk[:,1], Y_rk[:,0])
    plt.title("Portrait de phase avec la méthode de Runge-Kutta d'ordre 4")
    plt.xlabel("Teta")
    plt.ylabel("Teta'")
    plt.axis("equal")    
    plt.show()
    plt.plot(Yode[:,1], Yode[:,0])
    plt.title("Portrait de phase avec la fonction odeint")
    plt.xlabel("Teta")
    plt.ylabel("Teta'")
    plt.axis("equal")
    plt.show()


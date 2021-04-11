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
        plt.title("Tracé de θ(t) linéarisé")
        plt.xlabel("Temps (s)")
        plt.ylabel("Angle θ (rad)")
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
        plt.title("Tracé de θ(t) Euler")
        plt.xlabel("Temps (s)")
        plt.ylabel("Angle θ (rad)")
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
        plt.title("Tracé de θ(t) Runge-Kutta 4")
        plt.xlabel("Temps (s)")
        plt.ylabel("Angle θ (rad)")
        plt.show()
    return (Y_rk)
    

def Runge_Kutta_2(h,Ta, figure=False):
    N=int (Ta/h)
    Y_rk=np.zeros((N+1,2))
    Y_rk[0]=Y_0
    t=[0]
    for i in range(1, N+1):
        k=Y_rk[i-1]+pendule(Y_rk[i-1], t[i-1])*h*0.5
        Y_rk[i]=Y_rk[i-1]+h*pendule(k,t[i-1]+h/2)
        t.append(h*i)
    if (figure):
        plt.plot(t, Y_rk[:,1])
        plt.title("Tracé de θ(t) Runge-Kutta 2")
        plt.xlabel("Temps (s)")
        plt.ylabel("Angle θ (rad)")
        plt.show()
    return (Y_rk)  

    
def figures_teta(h,Ta):
    teta_e=Euler(h,Ta)[:,1]
    teta_l=Pendule_linearise(h, Ta)
    teta_rk4=Runge_Kutta_4(h, Ta)[:,1]
    teta_rk2=Runge_Kutta_2(h, Ta)[:,1]
    N=int(Ta/h)
    t=np.linspace(0,Ta,N+1)
    Yode =si.odeint(pendule, Y_0, t)
    plt.plot(t,teta_e,label="θ Euler")
    plt.plot(t, teta_l, label="θ linéarisé")
    plt.xlabel("Temps (s)")
    plt.ylabel("Angle θ (rad)")
    plt.plot(t, teta_rk4, label="θ Runge-Kutta 4")
    plt.plot(t, teta_rk2,label="θ Runge-Kutta 2")
    plt.plot(t, Yode[:,1],"+",label="Odeint")
    plt.title("Différent tracé de θ")
    plt.legend()
    plt.show()
    
  

def portrait_de_phase(h, Ta):
    Y_e=Euler(h,Ta)
    Y_rk4=Runge_Kutta_4(h, Ta)
    Y_rk2=Runge_Kutta_2(h, Ta)
    N=int(Ta/h)
    t=np.linspace(0,Ta,N+1)
    Yode =si.odeint(pendule, Y_0, t)
    plt.plot(Y_e[:,1], Y_e[:,0])
    plt.title("Portrait de phase avec la méthode d'Euler")
    plt.xlabel("θ")
    plt.ylabel("θ'")
    plt.axis("equal")
    plt.show()
    plt.plot(Y_rk4[:,1], Y_rk4[:,0])
    plt.title("Portrait de phase avec la méthode de Runge-Kutta d'ordre 4")
    plt.xlabel("θ")
    plt.ylabel("θ'")
    plt.axis("equal")    
    plt.show()
    plt.plot(Y_rk2[:,1], Y_rk2[:,0])
    plt.title("Portrait de phase avec la méthode de Runge-Kutta d'ordre 2")
    plt.xlabel("θ")
    plt.ylabel("θ'")
    plt.axis("equal")    
    plt.show()
    plt.plot(Yode[:,1], Yode[:,0])
    plt.title("Portrait de phase avec la fonction odeint")
    plt.xlabel("θ")
    plt.ylabel("θ'")
    plt.axis("equal")
    plt.show()
 
#------------------------------------------
#Suspension de vehicule
#------------------------------------------

Y_0_v=np.array([0,0,0,0])
C2=1200
M1=15
M2=200
K1=50000
K2=5000
f=-1000

def suspension(Y,t):
    y1=Y[0]
    y2=Y[1]
    y3=Y[2]
    y4=Y[3]
    L1=(1/M1)*(-C2*y1+C2*y2-(K1+K2)*y3+K2*y4)
    L2=(1/M2)*(f+C2*y1-C2*y2+K2*y3-K2*y4)
    L3=y1
    L4=y2
    return np.array([L1, L2, L3, L4])

def resolution_suspension(h,Ta, figure=False):
    N=int(Ta/h)
    t=np.linspace(0,Ta,N+1)
    Yode =si.odeint(suspension, Y_0_v, t)
    if (figure):
        plt.plot(t,Yode[:,2], label="roue x1(t)")
        plt.plot(t,Yode[:,3], label="caisse x2(t)")
        plt.xlabel("Temps (s)")
        plt.ylabel("Affaissement (m)")
        plt.title("Affaissement de la roue et de la caisse")
        plt.legend()
        plt.show()
    
#------------------------------------------
#Programme principal
#------------------------------------------

if __name__ == "__main__":
    print("La pulsation ω0 vaut : ω0 = " +str(w0))
    print("La fréquence f du mouvement est : f= "+str(f))
    print("La période T est : T = "+str(T))
    print("La constante A de l'équation du mouvement linéarisé est : A = "+str(A))
    print("La constante ϕ de l'équation du mouvement linéarisé est : A = "+str(phi))
    h=0.04
    Ta=4
    Pendule_linearise(h, Ta, figure=True)
    Euler(h, Ta, figure=True)
    Runge_Kutta_4(h, Ta, figure=True)
    Runge_Kutta_2(h,Ta, figure=True)
    figures_teta(h, Ta)
    portrait_de_phase(h, Ta)
    resolution_suspension(h, Ta, figure=True)

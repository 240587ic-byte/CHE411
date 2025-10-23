

# notation used
# eps is epsilon (porosity)
# P is pressure
# C is conc
# T is temp
# z is distance in reactor
# U velocity
# d is diameter
# om is omega (area)
# rho (density)
# F flow rate

# import modules
from scipy.integrate import solve_ivp
from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt
import math as m



# define some functions

def eps(dt,dp):
    eps=0.38+0.073*(1-(dt/dp-2)**2/(dt/dp)**2)
    return eps

def C(P,T):
    C=P/(8.314*T)
    return C

def P(C,T):
    P=C*8.314*T
    return P

def F(C,U,om):
    F=C*U*om
    return F

def om(dt):
    om=m.pi*dt**2/4
    return om

def U(Fo,To,Po,om):
    Uo=Fo*8.314*To/(om*Po)
    return Uo

def rhog(P,Mbar,T):
    rhog=P*Mbar/(8.314*T)
    return rhog

## Initial Conditions
L=12
dt=0.1
dp=0.01
rhop=2355.2
Cpp=950
lamp=0.3489
T0=793.15
P0=25.69
## CH4,CO,CO2,H2,H2O,N2
F0=np.array(5.17,0,0.29,0.63,17.35,0.85)
F0mols=np.array(18.612,0,1.044,2.268,62.46,3.06)
## Gas Constant
R=8.314


Ftotal0=F0mols.sum()
Uz0=(Ftotal0*R*T0)/(m.pi*dt**2/4*P0)

## Eta value (User input)
eta=0.1


def aaa():
    

def reaction_rates(Tp,P_partial):
    K=np.array(0,0,0)
    Hrxn=np.array(0,0,0)
    AK=np.array(4.707*10**12,1.142*10**-2,5.375*10**10)
    Hrxn298K=np.array(206.1,-41.15,164.9)

    for i in range(0,3):
        Hrxn[i]=Hrxn298K[i]+quad(aaa)
        K[i]=AK[i]*m.exp(-Hrxn[i])

    DEN=1+







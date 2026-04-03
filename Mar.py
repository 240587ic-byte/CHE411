import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def dpdw(alpha, p, FT, FT0):
    return -(alpha/(2*p))*(FT/FT0)

def rates(p, F, P0):
    F_CH4,F_H2O,F_CO,F_H2,F_CO2=F
    FT=sum(F)
    # Mole fractions
    y=[Fi/FT for Fi in F]
    # Actual pressure
    P=p*P0
    # Partial pressures
    P_CH4=y[0]*P
    P_H2O=y[1]*P
    P_CO =y[2]*P
    P_H2=y[3]*P
    P_CO2=y[4]*P
    # Rates
    R=0.08315*10**5 # UNITS   J  K**-1 kmol**-1
    T=1000 # Kelvins
    Ak=[4.225*10**15,1.955*10**6,1.02*10**15]
    AK=[4.707*10**12,1.142*10**-2,5.375*10**10]
    AKc=[8.23*10**-5,6.12*10**-9,6.65*10**-4,1.77*10**5]
    E=[240.1*10**6,67.13*10**6,243.9*10**6] # J/kmol
    del948_Hrxn=[224*10**6,-37.3*10**6,187.5*10**6] # J/kmol
    del_Hi=[-38.28,88.68,-70.65,-82.9]

    k=[Ak[i]*np.exp(-E[i]/(R*T)) for i in [0,1,2]]
    K=[AK[i]*np.exp(-del948_Hrxn[i]/(R*T)) for i in [0,1,2]]
    Kc=[AKc[i]*np.exp(-del_Hi[i]/(R*T)) for i in [0,1,2,3]]

    DEN=1+Kc[0]*P_CH4+Kc[1]*P_H2O/P_H2+Kc[2]*P_CO+Kc[3]*P_H2

    r1=k[0]/P_H2**2.5*(P_CH4*P_H2O-P_H2**3*P_CO/K[0])*(DEN**-2)
    r2=k[1]/P_H2*(P_CO*P_H2O-P_H2*P_CO2/K[1])*(DEN**-2)
    r3=k[2]/P_H2**3.5*(P_CH4*P_H2O**2-P_H2**4*P_CO/K[2])*(DEN**-2)

    return r1,r2,r3

def f_prime(W, y):
    p=y[0]
    F=y[1:]

    alpha = 0.001
    P0=50 # bar
    FT0=4.3
    FT=sum(F)
    r1,r2,r3=rates(p,F,P0)

    dF_CH4=-(r1+r3)
    dF_H2O=-(r1+r2+2*r3)
    dF_CO=(r1-r2)
    dF_CO2=(r2+r3)
    dF_H2=(3*r1+r2+4*r3)

    dp = dpdw(alpha,p,FT,FT0)

    return [dp,dF_CH4,dF_H2O,dF_CO,dF_CO2,dF_H2]

# Initial conditions
p0=1

# Feed (kmol/s)
F_CH4_0=1
F_H2O_0=3 
F_CO_0=0.1
F_CO2_0=0.1
F_H2_0=0.1


# Solve
L=3000
Wspan=(0,L)
W_eval=np.linspace(0,L,1000)

sol=solve_ivp(f_prime,Wspan,[p0,F_CH4_0,F_H2O_0,F_CO_0,F_CO2_0,F_H2_0],t_eval=W_eval)

p_sol=sol.y[0]
F_sol=sol.y[1:]

X=(F_CH4_0-F_sol[0])/F_CH4_0

plt.plot(sol.t, p_sol)
plt.xlabel('Catalyst Weight (kg)')
plt.ylabel('p = P/P0')
plt.title('Pressure Profile')
plt.show()

plt.plot(sol.t, X)
plt.xlabel('Catalyst Weight (kg)')
plt.ylabel('Conversion of CH4')
plt.title('Conversion Profile')
plt.show()

labels = ['CH4','H2O','CO','CO2','H2']
for i in range(5):
    plt.plot(sol.t, F_sol[i],label=labels[i])
plt.xlabel('Catalyst Weight (kg)')
plt.ylabel('Molar Flow (kmol/s)')
plt.legend()
plt.title('Species Profiles')
plt.show()
## Notes
# Word doc with equations and derivations
# Darft of presentation
# Link to repository
# convert W to L before graphing




import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# with heater

# reactor data
ID=0.1 #m
Length=12 #m
vol=ID**2*np.pi/4*Length
den_cat=2355.2 # kg/m**3
W_tot=vol*den_cat


# stoic coeff
nu1=[-1,-1,1,3,0]
nu2=[0,-1,1,1,1]
nu3=[-1,-2,0,4,1]
nu=[[-1,-1,1,3,0],[0,-1,1,1,1],[-1,-2,0,4,1]]
# Cp/R data assume ideal and constant ref 298K
Cp_data=[4.217,4.038,3.507,3.468,4.476]
Cp_CH4=4.217
Cp_H2O=4.038
Cp_CO=3.507
Cp_H2=3.468
Cp_CO2=4.476

def dpdw(alpha, p, FT, FT0):
    return -(alpha/(2*p))*(FT/FT0)

def rates(p, F, P0, T):
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
    Ak=[4.225*10**15,1.955*10**6,1.02*10**15]
    AK=[4.707*10**12,1.142*10**-2,5.375*10**10]
    AKc=[8.23*10**-5,6.12*10**-9,6.65*10**-4,1.77*10**5]
    E=[240.1*10**6,67.13*10**6,243.9*10**6] # J/kmol
    del948_Hrxn=[224*10**6,-37.3*10**6,187.5*10**6] # J/kmol
    del_Hi_init=[-38.28,88.68,-70.65,-82.9]
    del_Hrxn=[]
    del_Hrxn.append(sum(Cp_data[0]*8.314*nu[0][i]*(T-298)*10**3 for i in [0,1,2,3,4])+del948_Hrxn[0])
    del_Hrxn.append(sum(Cp_data[1]*8.314*nu[1][i]*(T-298)*10**3 for i in [0,1,2,3,4])+del948_Hrxn[1])
    del_Hrxn.append(sum(Cp_data[2]*8.314*nu[2][i]*(T-298)*10**3 for i in [0,1,2,3,4])+del948_Hrxn[2])
    del_Hi=[del_Hi_init[i]+Cp_data[i]*(T-298) for i in [0,1,2,3]]
    k=[Ak[i]*np.exp(-E[i]/(R*T)) for i in [0,1,2]]
    K=[AK[i]*np.exp(-del_Hrxn[i]/(R*T)) for i in [0,1,2]]
    Kc=[AKc[i]*np.exp(-del_Hi[i]/(R*T)) for i in [0,1,2,3]]

    DEN=1+Kc[0]*P_CH4+Kc[1]*P_H2O/P_H2+Kc[2]*P_CO+Kc[3]*P_H2


    r1=k[0]/P_H2**2.5*(P_CH4*P_H2O-P_H2**3*P_CO/K[0])*(DEN**-2) # units in kmol/(h kg catalyst)
    r2=k[1]/P_H2*(P_CO*P_H2O-P_H2*P_CO2/K[1])*(DEN**-2)
    r3=k[2]/P_H2**3.5*(P_CH4*P_H2O**2-P_H2**4*P_CO/K[2])*(DEN**-2)

    return r1,r2,r3,del_Hrxn


def temp(r1,r2,r3,F,T,del_Hrxn):
    ror=[r1,r2,r3]
    RH=sum(ror[i]*del_Hrxn[i] for i in [0,1,2])
    FCp=sum(F[i]*Cp_data[i]**8.314*1**1 for i in [0,1,2,3,4])
    Ta=1500
    Ua=1
    heater=Ua*(T-Ta)
    #heater=0
    dtdw=(RH-heater)/FCp
    return dtdw


def f_prime(W, y):
    p=y[0]
    T=y[1]  # Kelvins
    F=y[2:]

    alpha = 0.00004
    P0=1 # bar
    FT0=40.3
    FT=sum(F)
    r1,r2,r3,del_Hrxn=rates(p,F,P0,T)
    #r1,r2,r3,del_Hrxn=rates(p,F,P0,T)

    dF_CH4=-(r1+r3)
    dF_H2O=-(r1+r2+2*r3)
    dF_CO=(r1-r2)
    dF_CO2=(r2+r3)
    dF_H2=(3*r1+r2+4*r3)

    dp = dpdw(alpha,p,FT,FT0)
    dT = temp(r1,r2,r3,F,T,del_Hrxn)
    #print(dT)
    return [dp,dT,dF_CH4,dF_H2O,dF_CO,dF_H2,dF_CO2]

# Initial conditions
p0=1
T0=1000
# Feed (kmol/s)
F_CH4_0=10
F_H2O_0=30 
F_CO_0=0.1
F_CO2_0=0.1
F_H2_0=0.1


# Solve
L=5000
Wspan=(0,W_tot)
W_eval=np.linspace(0,W_tot,1000)

sol=solve_ivp(f_prime,Wspan,[p0,T0,F_CH4_0,F_H2O_0,F_CO_0,F_H2_0,F_CO2_0],t_eval=W_eval)

p_sol=sol.y[0]
T_sol=sol.y[1]
F_sol=sol.y[2:]
L=[W*4/(den_cat*np.pi*ID**2) for W in sol.t]
X=(F_CH4_0-F_sol[0])/F_CH4_0
x=(F_sol[4]-F_H2_0)/(F_CH4_0-F_sol[1])

plt.plot(L, p_sol)
plt.xlabel('Catalyst Weight (kg)')
plt.ylabel('p = P/P0')
plt.title('Pressure Profile')
plt.show()

plt.plot(sol.t, T_sol)
plt.ticklabel_format(style='plain',useOffset=False)
plt.xlabel('Catalyst Weight (kg)')
plt.ylabel('Temp (K)')
plt.title('Temperature Profile')
plt.show()

plt.plot(sol.t, X)
plt.xlabel('Catalyst Weight (kg)')
plt.ylabel('Conversion of CH4')
plt.title('Conversion Profile')
plt.show()

plt.plot(sol.t, x)
plt.xlabel('Catalyst Weight (kg)')
plt.ylabel('Ratio H2 Generated / CH4 Used')
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
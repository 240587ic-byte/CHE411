import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

### steam methane reaction
### pressure and conc/conversion profile
### elem non-rev rate law, use partial pressures
### rev rate

### Ergun Equation
def dpdz(P,P0,T,T0,epsilon,X,mu,rho,gc,Dp,phi,MU):
    return (-rho*mu/(rho*gc*Dp)*((1-phi)/phi**3))((150*(1-phi)*MU/Dp)+1.75*rho*mu)*P0*T/(P*T0)*(1+epsilon*X)


## Isothermal
def dpdw(alpha,p,epsilon,X):
    return (-alpha/(2*p))*(1+epsilon*X)

# Rate of rxn
# def ra(k,conc):
#     cinit=1000
#     return k*cinit*(1-X)/()

def dxdw(p,X,cinit,k,FA0,epsilon):
    return (k*cinit)*(1-X)/(1+epsilon*X)*p/FA0

## for ethane ---> ethylene + hydrogen rxn at 1000K

def f_prime(z,data):
    alpha=0.001
    epsilon=1
    p0=50 #atm
    k=0.072
    v0=10 #liters/sec
    cinit=p0/(0.082057*1000)
    fa0=cinit*v0
    p,X=data
    return [dpdw(alpha,p,epsilon,X),dxdw(p,X,cinit,k,fa0,epsilon)]


L=300 # kg of catalyst
dx=np.linspace(0,L,100)

sol=solve_ivp(f_prime,t_span=(0,L),y0=[1,0.0],t_eval=dx)
plt.plot(sol.t,sol.y[0])
plt.show()
plt.plot(sol.t,sol.y[1])
plt.show()
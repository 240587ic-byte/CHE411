### Reactions
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

### mainly PFR as CSTR is uniform


### assume perfect radial mixing no axial mixing
# from scipy.integrate import solve_ivp
# import numpy as np
# import matplotlib.pyplot as plt

# def ror(conc):
#     return -0.05*conc**2


# def system(t,conc):
#     rate_of_rxn=ror(conc)
#     return rate_of_rxn
# init=100
# sol=solve_ivp(system,t_span=(0,1),y0=[100],t_eval=np.linspace(0,1,100))
# print(sol)
# a=sol.y[0],sol.y[0]
# print(len(a))
# print(len(a[0]))
# print(sol.t)
# plt.plot(sol.t,sol.y[0])
# plt.show()
# plt.imshow(a, extent=[0,1,0,0.05], cmap='inferno',aspect='auto')
# plt.show()


### Ergun Equation
def dpdz(P,P0,T,T0,epsilon,X,mu,rho,gc,Dp,phi,MU):
    return (-rho*mu/(rho*gc*Dp)*((1-phi)/phi**3))((150*(1-phi)*MU/Dp)+1.75*rho*mu)*P0*T/(P*T0)*(1+epsilon*X)

def dpdw(alpha,p,T,T0,epsilon,X):
    return (-alpha/(2*p))*(1+epsilon*X)*(T/T0)

# Rate of rxn
def ra(k,conc):
    return k*conc

def dxdw(FA0,k,conc):
    return -ra(k,conc)/FA0







### use ergun equation, use interms of conversion and PBR. page 187

### no radial mixing and a no slip boundary layer with laminar profile
### liquid phase so change in moles will not be needed



def ror(conc):
    return -0.005*conc**1.5

def dconcdx(x,conc,vmax,r,R):
    return ror(conc)/v(vmax,r,R)

def v(vmax,r,R):
    return vmax*(1-(r/R)**2)



L=10 
R=0.1
cinit=1000 
vmax=0.1


dr=np.linspace(R-0.000001,0,100)
dx=np.linspace(0,L,1000)
bothalf=[]

for r in dr:
    sol=solve_ivp(dconcdx,t_span=(0,L),y0=[cinit],t_eval=dx,args=(vmax,r,R))
    bothalf.append(sol.y[0])

tubesol=bothalf

for i in range(0,len(dr)-1):
    tubesol.append(bothalf[len(dr)-1-i])

plt.imshow(bothalf,extent=[0,L,-R,R],cmap="inferno",aspect="auto")
plt.colorbar(label="Concentration")
plt.xlabel("Length (m)")
plt.ylabel("Radius (m)")
plt.show()


# import numpy as np
# import matplotlib.pyplot as plt

# # Example 3D concentration field
# x = np.linspace(0, 1, 50)
# y = np.linspace(0, 1, 50)
# z = np.linspace(0, 1, 50)
# X, Y, Z = np.meshgrid(x, y, z)

# C = np.exp(-10*((X-0.5)**2 + (Y-0.5)**2 + (Z-0.5)**2))

# # Take a slice at z = 0.5
# slice_index = 25
# C_slice = C[:, :, slice_index]
# print(C_slice)
# plt.imshow(C_slice, extent=[0,1,0,1], origin='lower', cmap='inferno')
# plt.colorbar(label="Concentration")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.title("Concentration Slice at z = 0.5")
# plt.show()

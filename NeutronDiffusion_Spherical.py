import numpy as np
from scipy.integrate import odeint, solve_bvp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

'''----------------------------------------------------------------------'''

# k as a function of alfa
def k(a):
    return np.sqrt((eta+a)/mu)

# alfa as a function of k
def a(k):
    return mu*k**2-eta

#Analytical relation k(a)
def f(k,r):
    lam_t=0.0360
    #lam_t=0.0411 for Plutonium Analysis
    return k*r*(np.tan(k*r))**(-1)+3/2*r/lam_t-1

#Criticality condition alfa=0
def f_crit(r):
    k_crit=k(0)
    f_crit=f(k_crit,r)
    return f_crit

#Useful function to find the sovracritical k
def f_k(k, R):
    f_k=f(k,R)
    return f_k

#Setting for the Radial ODE
def equation(r,R):
    dR=R[1]
    ddR=-(2/r)*R[1]-k_sc**2*R[0]
    return [ dR, ddR]

#This sets the final condition as dR=-3/2*R/lam and the initial condition (at the origin) = k as the limit of the analytic solution
def boundary_conditions(Ra, Rb):
    return np.array([Ra[0]-k_sc, Rb[1]+(3/(2*lam_t))*Rb[0]], dtype=object)

#Model of the temporal differential equation
def model(y, t, alfa):
    dydt=-alfa*y
    return dydt

#Time ODE solver
def T_ode(t, alfa): 
    t_vec=np.linspace(0, t, nt)
    #initial condition
    y0=1
    T_solution=np.squeeze(odeint(model, y0, t_vec, args=alfa))
    return T_solution

'''----------------------------------------------------------------------'''

#Uranium
mu=2.3446*10**5
eta=1.8958*10**8
lam_t=0.0360

'''
For the analysis with Plutonium instead of Uranium it is sufficient to change
these three constants and to adjust properly the 
values of the boundaries using as a guide the new resulting critical lenght
mu=2.6786*10**5
eta=3.0055*10**8
lam_t=0.0411

'''

#Initialization of the problem
nr, nt=100, 100
Radius=0.09 #Radius
t_f=7*10**(-7) #Final time

#Discretization of space and time
r=np.linspace(0.0001,Radius,nr)
t=np.linspace(0,t_f,nt)

#Critical radius
r_crit=fsolve(f_crit,0.05)

#Parameters for a sovracritical radius
k_sc=fsolve(f_k, k(0), args=(Radius))
a_sc=a(k_sc)

#Initial Guess for R and dR
R_guess=np.zeros((2,nr))
R_guess[0]=1

#Solution of the differential equation for R
solution= solve_bvp(equation, boundary_conditions, r, R_guess)
R_sol=solution.sol(r)[0]

#Solution of the time differential equation
alfa=tuple(a_sc)
T_sol=T_ode(t_f, alfa)

#Combined solution prepared for the grid
n_final=np.outer(T_sol, R_sol)

#Plot of the surface
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection='3d')
X, T = np.meshgrid(r,t)  
surf = ax.plot_surface(X, T, n_final, cmap='inferno')
ax.set_xlabel('R [m]', fontsize=15)
ax.set_ylabel('$t \cdot 10^{-7} [s]$', fontsize=15)
ax.set_zlabel('n(t,R)', fontsize=15)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()







    


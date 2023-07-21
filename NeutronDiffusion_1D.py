import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, trapz

'''----------------------------------------------------------------------'''

#Spatial Ordinary Differential Equation solver. Here we solve the Eigenvalue problem
def X_ode(L, mu, N): 
    #Discretization of space
    x=np.linspace(0,L,nx)
    h=x[1]-x[0]
    #Creation of Kinetic Matrix
    T=np.zeros((nx-2)**2).reshape(nx-2,nx-2)
    for i in range(nx-2):
        for j in range(nx-2):
            if i==j:
                T[i,j]=-2
            elif np.abs(i-j)==1:
                T[i,j]=1
            else:
                T[i,j]=0
    T*=(-mu/(h**2))
    #Eigenvalues and Eigenvectors
    val, vec=np.linalg.eig(T)
    #Ordering of Eigenvalues
    z=np.argsort(val)
    #First N Eigenvalues
    z=z[0:N]
    # Addition of zeros at the borders of Eigenvectors
    vec_new = np.zeros((nx, len(z)))
    for i in range(len(z)):
       vec_new[:, i] = np.concatenate(([0], vec[:, z[i]], [0]))
    return val[z], vec_new

#Usefull rewrite of Eigenvectors and rescaling
def separation(Xi,N):
    X_q=[]
    for i in range(N):
        X_q_i=Xi[:,i]
        X_q.append(X_q_i)
        #Here we rescale the Eigenvectors as the sin function 
        scale=1/(max(abs(X_q[i])))
        X_q[i]*=scale
    return X_q

#Model of the temporal differential equation
def model(y, t, eta, alfa):
    dydt=(eta-alfa)*y
    return dydt

#Time ODE solver
def T_ode(t, eta, alfa): #t is a float Number
    t_vec=np.linspace(0, t, nt)
    #Initial condition
    y0=1
    T_solution=np.squeeze(odeint(model, y0, t_vec, args=(eta,alfa)))
    return T_solution

#Critical length finder
def find_the_boom(mu, eta):
    #accuracy
    e=10**(-4)
    #Shift in L (1mm)
    dL=0.001 
    #L temporal for loops
    L_0, L_boom=0, 0
    #Loop that search for the change in the slope of the time solution
    while dL>=e:
        L_boom=L_0+dL
        #First we solve the X ode to find the minimum eigenvalue
        alfas=X_ode(L_boom, mu, N)[0]
        #Minimum eigenvalue: the worst
        alfa_min=alfas[0]
        T_0=T_ode(t_f, eta, alfa_min)
        T_0_der=np.gradient(T_0)
        sign_check=np.average(T_0_der)
        #Evaluation of the change in the slope of the solution
        if sign_check>0:
            dL=dL/2
        else:
            L_0=L_boom
    return L_boom

#Theoretical prediction
def L_critical(mu,eta):
    return np.pi*np.sqrt(mu/eta)

#Initial condition
def f(x):
    lam=100
    l=L/2
    f=np.exp(-lam*((x-l)**2/(l**2)))
    return f

#Coefficients of the series expansion
def coefficients(f,X,N):
    ap=np.zeros(N)
    cost=2/L
    for i in range(N):
        fun= f*X[i]
        coeff= cost * trapz(fun, x)
        ap[i]= coeff
    return ap

#Combination of space an time solutions (without coefficients yet)
def Space_time_eig(t, X, A, N):
    for  i in range(N):
        X[i]*=T_ode(t, eta, A[i])
    return X

#Complete solution with coefficients
def function_exp(ap,t,X,A, N):
    X_st=Space_time_eig(t,X,A,N)
    for i in range(N):
        X_st[i]*=ap[i]
    return X_st

'''----------------------------------------------------------------------'''

#Uranium
mu=2.3446*10**5
eta=1.8958*10**8

'''
For an analysis of the Plutonium it is sufficient to change
these two values of the constants and to adjust properly the 
values of the boundaries using as a guide the new resulting critical lenght
mu=2.6786*10**5
eta=3.0055*10**8

'''

#Initialization of the problem
nx, nt=100, 100  #number of space and time points
L=0.12 #Interval
t_f=4*10**(-7) #Final time 
N=30 #Number of eigenvalues and eigenvectors

#Discretization of space and time
x=np.linspace(0,L,nx)
t=np.linspace(0,t_f,nt)

#What's the critical length?
L_crit=find_the_boom(mu, eta)
print('Theoretical critical Length:  %.3f' %(L_critical(mu, eta)*100), 'cm')
print('Numerical Solution:  %.3f' %(L_crit*100), 'cm')

#Boundary condition
f=f(x)

#Differential Spatial Equation
A_q,Xi=X_ode(L, mu, N)
X_q=separation(Xi, N)

#Coefficients
ap=coefficients(f, X_q, N)

#Solution, function expansion of the form n_p*a_p 
n_p=np.squeeze(function_exp(ap, t_f, X_q, A_q, N))

#Summation of the seires expansion
n_t=np.sum(n_p,axis=0)

#Starting smooth function that respect BCs
F=np.sin(np.pi*(x/L))
n_mat=np.outer(n_t,F)

#Plot of the evolution of the surface
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection='3d')
X, T = np.meshgrid(x,t)  
surf = ax.plot_surface(X, T, n_mat, cmap='inferno')
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel('x [m]', fontsize=15)
ax.set_ylabel('$t \cdot 10^{-7} [s]$', fontsize=15)
ax.set_zlabel('n(t,x)', fontsize=15)
plt.show()




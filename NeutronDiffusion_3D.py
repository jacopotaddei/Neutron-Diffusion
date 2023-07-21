import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, trapz

'''----------------------------------------------------------------------'''

#Spatial Ordinary Differential Equation solver. 
#Here we solve the Eigenvalue problem for the X ordinary differential equation
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
    q=np.argsort(val)
    #First N Eigenvalues
    q=q[0:N]
    #Addition of zeros
    vec_new = np.zeros((nx, len(q)))
    for i in range(len(q)):
       vec_new[:, i] = np.concatenate(([0], vec[:, q[i]], [0]))
    return val[q], vec_new

#Here we solve the Eigenvalue problem for the Y ordinary differential equation
def Y_ode(L, mu, N): 
    #Discretization of space
    y=np.linspace(0,L,ny)
    h=y[1]-y[0]
    #Creation of Kinetic Matrix
    T=np.zeros((ny-2)**2).reshape(ny-2,ny-2)
    for i in range(ny-2):
        for j in range(ny-2):
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
    q=np.argsort(val)
    #First N Eigenvalues
    q=q[0:N]
    # Addition of zeros
    vec_new = np.zeros((ny, len(q)))
    for i in range(len(q)):
       vec_new[:, i] = np.concatenate(([0], vec[:, q[i]], [0]))
    return val[q], vec_new

#Here we solve the Eigenvalue problem for the Z ordinary differential equation
def Z_ode(L, mu, N): 
    #Discretization of space
    z=np.linspace(0,L,nz)
    h=z[1]-z[0]
    #Creation of Kinetic Matrix
    T=np.zeros((nz-2)**2).reshape(nz-2,nz-2)
    for i in range(nz-2):
        for j in range(nz-2):
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
    q=np.argsort(val)
    #First N Eigenvalues
    q=q[0:N]
    # Addition of zeros
    vec_new = np.zeros((nz, len(q)))
    for i in range(len(q)):
       vec_new[:, i] = np.concatenate(([0], vec[:, q[i]], [0]))
    return val[q], vec_new

#Usefull rewrite of Eigenvectors and rescaling
def separation(Xx):
    X=[]
    for i in range(N):
        X_i=Xx[:,i]
        X.append(X_i)
        #Here we rescale the Eigenvectors as the sin function 
        scale=1/(max(abs(X[i])))
        X[i]*=scale
    return np.squeeze(X)

#Model of the temporal differential equation
def model(y, t, eta, alfa):
    dydt=(eta-alfa)*y
    return dydt

#Time ODE solver
def T_ode(t, L, eta, alfa):
    t_vec=np.linspace(0,t,nt)
    #Initial condition
    y0=1
    y1=np.squeeze(odeint(model, y0, t_vec, args=(eta,alfa)))
    return y1

#Critical length finder
def find_the_boom(mu, eta):
    #Accuracy
    e=10**(-4)
    #Shift in L (1mm)
    dL=0.001 
    #L temporal for loops
    L_0, L_boom=0, 0
    #Loop that search for the change in the slope of the time solution
    while dL>=e:
        L_boom=L_0+dL
        #First we solve the X, Y, Z ODEs to find the minimum eigenvalue
        alfa_x=X_ode(L_boom, mu, N)[0]
        alfa_y=Y_ode(L_boom, mu, N)[0]
        alfa_z=Z_ode(L_boom, mu, N)[0]
        alfas=alfa_x+alfa_y+alfa_z
        #Minimum eigenvalue: the worst
        alfa_min=alfas[0]
        T_0=T_ode(t_f,L_boom,eta, alfa_min)
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
    return np.pi*np.sqrt(3*mu/eta)

#Initial condition
def f(x,y,z):
   return (8/L**3)*x*y*z*(1-x/L)*(1-y/L)*(1-z/L)

#Initial condition restricted at z=L/2
def f2(x,y):
    z=L/2
    return f(x,y,z)

#Coefficients of the series expansion
def coef_finder(N):
    ap=np.zeros((N,N,N))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                fx=x*(1-x/L)*X[i]
                fy=y*(1-y/L)*Y[j]
                fz=z*(1-z/L)*Z[k]
                ap[i, j, k] =(8/L**3)**2 * trapz(fx,x)*trapz(fy,y)*trapz(fz,z)
                #print(i+1,j+1,k+1, '---------->', ap[i,j,k])
    return ap

#Combination of all of the spatial solutions
def combined_eigenvec(X,Y,Z,N):
    R=np.zeros((N,N,N,100))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                R[i,j,k]=X[i]*Y[j]*Z[k]
    return R

#Combination of space an time solutions (without coefficients yet)
def Space_time_eig(t,R,Ax,Ay,Az,N):
    for i in range(N):
        for j in range(N):
            for k in range(N):
                alfa=Ax[i]+Ay[j]+Az[k]
                R[i,j,k]*=T_ode(t, L, eta, alfa) 
    return R

#Complete solution with coefficients
def function_exp(ap,t,R,Ax,Ay,Az,N):
    R_st=Space_time_eig(t, R, Ax, Ay, Az, N)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                R_st[i,j,k]*=ap[i,j,k]
    return R_st

#Function that sums all of the components of n(t,x,y,z)
def summed(n_p):
    n_sum=np.sum(n_p, axis=2)
    n_sum=np.sum(n_sum, axis=1)
    n_sum=np.sum(n_sum,axis=0)
    return n_sum

'''----------------------------------------------------------------------'''

#Uranium
mu=2.3446*10**5
eta=1.8958*10**8

'''
For an analysis of the Plutonium it is sufficient to change
these two values of the constants and to adjust properly the 
values of the boundaries using as a guide the new value 
of the critical lenght
mu=2.6786*10**5
eta=3.0055*10**8

'''

#Initializzation of the problem
nx, ny, nz, nt=100,100,100,100 #number of space and time points
L=0.25 #Interval
t_f=8*10**(-7) #Final time 
N=4 #Number of eigenvalues and eigenvectors

#Discretization of space and time
x=np.linspace(0,L,nx)
y=np.linspace(0,L,ny)
z=np.linspace(0,L,nz)
t=np.linspace(0,t_f,nt)

#What's the critical length?
L_crit=find_the_boom(mu, eta)
print('Theoretical critical Length:  %.3f' %(L_critical(mu, eta)*100), 'cm')
print('Numerical Solution:  %.3f' %(L_crit*100), 'cm')

#Differential Spatial Equations
Ax,Xx=X_ode(L, mu, N)
X=separation(Xx)
Ay,Yy=Y_ode(L, mu, N)
Y=separation(Yy)
Az,Zz=Z_ode(L, mu, N)
Z=separation(Zz)

#Coefficients
ap=coef_finder(N)

#Combined eigenvectors
R=combined_eigenvec(X, Y, Z, N)

#Space-time combined eigenvectors
n_p=function_exp(ap, t_f, R, Ax, Ay, Az, N)

#Total solution
n_t=summed(n_p)


#Creation of space x-y grid
X_g, Y_g = np.meshgrid(x, y)
#Plot initial condition
#Slice at z=L/2 of the initial function and plot of the latter t=0
Z0=f2(X_g,Y_g)
plt.figure()
plt.imshow(Z0, extent=[0, L, 0, L], origin='lower', cmap='inferno')
plt.colorbar(label='f(x,y,L/2)')
plt.xlabel('x [m]', fontsize=15)
plt.ylabel('y [m]', fontsize=15)
plt.show()

#Plot n(t,x,y,z) in a slice z=L/2 
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection='3d')
'''Here we consider the evolution considering a particular value of n_t which is the max. 
This procedure is just to have a clear view of the function on the x-y plane.
If you want to find the value at a particular time you can search inside 
the vector n_t the precise value of the density relative to that time and multiply by
this particular value, indeed n_t represent the pure time evolution of the neutron density'''
Z_evolved=f2(X_g,Y_g)*max(n_t)
#Time associated time to this max
t_max= t[np.argmax(n_t)]
surf = ax.plot_surface(X_g, Y_g, Z_evolved, cmap='inferno')
ax.set_xlabel('x[m]', fontsize=15)  
ax.set_ylabel('y[m]', fontsize=15)  
ax.set_zlabel('n(x,y)', fontsize=15)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()





# Neutron-Diffusion

This project is about the diffusion of neutrons inside radioactive materials, such as $U^{235}$ or $Pu^{239}$. During the propagation of neutrons the collisions with nucleai result in a secondary realease of these particles, and if the quantity of fissile material is above a critical amount this secondary realease lead to an exponential growth in the neutron density and thus in a nuclear explosion. 
The equation that governs this proces is the following:
$$\cfrac{\partial n}{\partial t}= \mu \nabla^2 n + \eta n  $$
In this repository we analyzed the solution to the latter equation in three different systems of coordinates: 1D Cartesian, 3D Cartesian and 3D Spherical.

## Neutron diffusion 1D Cartesian coordinates - Dirichlet BC

This first system is made by a one dimensional chain of atoms of fissile material. As reported in the title here we used Dirichlet boundary conditions, i.e. the neutrons cannot escape from the system. 

Every part of the code is commented in order to have an idea while reading it of the functioning of each part, however in general the program is based on two important parameters that can be modified by hand: the space and time intervals.
The program always computes the value of the critical length after which we obtain an uncontrolled growth of the neutron diffusion, then it is possible to modify directly in the code the values L (by default setted as L=0.12) that represent the length of the system and t_f ( by default as t_f=4.0E-07) which is the final time. The solution will be computed and plotted in an interval $x\in[0,L]$ and $t\in[0,t_f]$.

For details regarding the Theoretical considerations of the problem and for a brief description of the Numerical analysis it is possible to see the pdf named Neutron_Diffusion.pdf, for a deeper description of the problem it is recomended the article used as a reference.

The code is written in Python and you need to install (if not present yet) the following libraries:
+ numpy
+ matplotlib
+ scipy

## Neutron diffusion 3D Cartesian coordinates - Dirichlet BC

This is the straight generalization fo the first case but in a three dimensional volume. We consider thus a cubic volume of fissile material and we impose also here the Dirichlet boundary conditions, i.e. no neutron can escape from the system.

The structure of the code is very similar to the one dimensional case and the equations are solved almost in the same way. Also here every part of the code is commented for the sake of clarity. The two more relevant parameters are the length L (by default L=0.25) and t_f (by default t_f=8.0E-07). In this case the length L represent the interval for all of the three coordinates, in the sense that now $x\in[0,L]$, $y\in[0,L]$, $z\in[0,L]$ , and as in the previous case $t\in[0,t_f]$. The critical length is computed too and can be converted into the information about the mass of the fissile material needed to create the explosion simply by using the Uranium density $\rho_U=18.7 g/cm^3$ (or the Plutonium density $\rho_{Pu}=15.6 g/cm^3$ if the analysis is performed on this material).

For details regarding the Theoretical considerations of the problem and for a brief description of the Numerical analysis it is possible to see the pdf named Neutron_Diffusion.pdf, for a deeper description of the problem it is recomended the article used as a reference.

The code is written in Python and you need to install (if not present yet) the following libraries:
+ numpy
+ matplotlib
+ scipy

## Neutron diffusion 3D Spherical coordinates - Neumann BC

This case is diffrent from the previous two. At first the volume of the fissile material is of Spherical shape and this will produce different equations (in form) to be solved, then here we impose Neumann bondary conditions instead of the previous simple Dirichlet. The Neumann boundary conditions are of the form $\cfrac{dR}{dr}=-\cfrac{3}{2}\cfrac{R}{\lambda_t}$ where $R$ is the function that describes the spatial evolution of the neutron density and $\lambda_t$ is a parameter specified in the refences and its value is reported directly inside the code. 

The two relevant parameters here are R that represent the value of the radius of the sphere (by default R=0.09) and as always the limit of the time interval t_f (by default t_f=7.0E-07), these can be changed directly inside the code to see how does the neutron density changes as a function of these two parameters.
The solution is visualized in an interval $t\in[0,t_f]$ and $r\in[0,R]$. This time the critical value of the radius is founded using analytical considerations.

For details regarding the Theoretical considerations of the problem and for a brief description of the Numerical analysis it is possible to see the pdf named Neutron_Diffusion.pdf, for a deeper description of the problem it is recomended the article used as a reference.

The code is written in Python and you need to install (if not present yet) the following libraries:
+ numpy
+ matplotlib
+ scipy

# References
Graham Griffiths. Neutron diffusion.02 2018. URL https://www.researchgate.net/publication/323035158_Neutron_diffusion





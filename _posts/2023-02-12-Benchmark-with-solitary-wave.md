---
title: 'Benchmark with solitary wave'
date: 2023-02-12
permalink:
tags:
  - Fluid dynamics
  - solitary wave
---



It is critical to ensure that our code is performing correctly.

We can verify it by benchmarking our code to other problems.

Simpson and Spiegelman (2018) provide **solitary wave benchmarks** in Magma Dynamics, which are excellent for benchmarking our code. The primary characteristic of a solitary wave is that it **maintains its shape and speed** as it propagates in space. 

If we put the single solitary wave as initial condition, with the compaction pressure, the wave should rise up without changing the shape or velocity. Therefore all that remains is to examine the distortion of the transported waveform.

## Governing equations for solitary wave

The governing equations for solitary wave are similar to the previously stated conservation of mass and Darcy’s law.

And If we assume no shear and no melting (reaction), the eqs are as follows,

$$
-{\partial \phi \over \partial t} + \nabla \cdot v_s-\nabla \cdot (\phi v_s)=0
$$

$$
{\partial \phi \over \partial t} + \nabla \cdot (\phi v_l)=0
$$

$$
\phi (v_l-v_s) = K ((1-\phi) \Delta\rho g\hat z -\nabla P)
$$

These equations could be dimensionless in a different way than the previous post.

## **Characteristic scales**

These equations have an intrinsic length scale, the compaction length,

$$
\delta = \sqrt{K_0 \zeta}
$$

and an intrinsic velocity scale, the Darcy’s flux (percolation flux)

$$
\phi_0 w_0 = K_0 \Delta \rho g
$$

where $$K_0$$ is the mobility at porosity $$\phi_0$$,

$$v_l=w_0 v_l'$$,

$$v_s=\phi_0 w_0 v_s’$$,

The compaction pressure $$P = \zeta C$$, where the compaction rate $$C=\nabla \cdot v_s$$,

Natural scaling,

$$(x,z) = \delta (x’,z') $$, 

$$\phi = \phi_0 \phi$$,

$$t={\delta \over w_0} t’$$,

$$K =K_0 K’$$,

$$ P = P_0 P’ = {\zeta \phi_0 w_0  \over\delta}  P'$$,



## Non-dimensionalization

### 1.The dimensionless mass conservation in the solid

eq(1), $$-{\partial \phi \over \partial t} + \nabla \cdot v_s-\nabla \cdot (\phi v_s)=0$$ becomes,

$$
-{\phi_0 w_0 \over \delta}{\partial \phi \over \partial t} + {\phi_0 w_0\over \delta}\nabla \cdot v_s-{\phi_0^2 w_0 \over \delta}\nabla \cdot (\phi v_s)=0
$$

$$\phi_0 $$ is small, we can neglect the last term and get

$$
{\partial \phi \over \partial t} = \nabla \cdot v_s
$$

###  2.The dimensionless mass conservation in the liquid

eq(2), $${\partial \phi \over \partial t} + \nabla \cdot (\phi v_l)=0$$ becomes, 

$$
{\phi_0 w_0 \over \delta} {\partial \phi \over \partial t} +{\phi_0 w_0\over \delta}\nabla \cdot (\phi v_l)=0
$$

divided by $$\phi_0 w_0$$,

$$
{\partial \phi \over \partial t} + \nabla \cdot (\phi v_l)=0
$$

### 3.The dimensionless Darcy’s law

eq(3) $$\phi (v_l-v_s) = K ((1-\phi) \Delta\rho g\hat z -\nabla P)$$ becomes,

$$
\phi_0 w_0 \phi (v_l-\phi_0 v_s) = K_0 K ((1-\phi_0\phi) \Delta\rho g\hat z -{\zeta \phi_0 w_0  \over\delta}  \nabla P)
$$

$$\phi_0 $$ is small, we can neglect the term $$ -\phi_0^2 w_0 \phi v_s$$, and $$-K_0 K \phi_0 \phi \Delta \rho g \hat z$$,

$$
\phi_0 w_0 \phi v_l= K_0 K (\Delta\rho g\hat z -{\zeta \phi_0 w_0  \over\delta^2}  \nabla P)
$$

velocity scale $$\phi_0 w_0 = K_0 \Delta \rho g$$, so eq(11) becomes,

$$
\phi v_l= K (\hat z - {\zeta K_0\over\delta^2} \nabla P)
$$

compaction length $$\delta = \sqrt{\zeta K_0}$$, 

$$
\phi v_l= K (\hat z - \nabla P)
$$


## Dimensionless Compressible Flow Equations

$$
{\partial \phi \over \partial t} = C
$$

$$
C + \nabla \cdot [K (\hat z - \nabla P)]=0
$$

That is, porosity only changes by compaction. The compaction rate is controlled by the divergence of the melt flux and the viscous resistance of the matrix to volume changes.

## Weak form of the governing equations

$$
\int_\Omega \phi R dx =\int_\Omega [\phi^n +P^n \Delta t] R dx
$$

$$
\int_\Omega P P_1 -\int_\Gamma P_1 K \nabla P\cdot n + \int_\Omega K \nabla P_1 \cdot \nabla P = \int_\Omega K {\partial P_1 \over \partial z} -\int_\Gamma P_1 K(\widehat z\cdot n)
$$

## Implementation

We use Marc Spiegelman's initSolitaryWave.py and Gideon Simpon's PySolwave libraries to test our code,

First, the fenics and libraries are imported,

```
from fenics import *
import numpy
import os, sys, subprocess
sys.path[0] += "/PySolwave"
from magmasinc import solwave_mck
from magmasinc import solwave_con
from sinc_eo import sincinterp_e
```

Then, we define solitary wave,

```
def dolfinsolwave(Q,n,m,h_on_delta,c,d,wd,z0,numcol):
    # Q -- Dolfin FunctionSpace
    # n -- coefficient in permeability relationship
    # m -- coefficient in bulk viscosity relationship
    # h_on_delta -- Domain width over compaction length
    # c -- solitary wave speed
    # d -- space dimension
    # wd -- wave dimension
    # z0 -- fractional position of wave in box
    # numcol -- number of collocation points
    
    # set up "radial function" in wd dimensions
    if wd == 1:
        str = "sqrt(pow((x[%d-1]-%g),2))" % (d,z0) 
    elif wd == 2:
        str = "sqrt(pow((x[0]-0.5),2) + pow((x[%d-1]-%g),2))" % (d,z0)
    elif wd == 3:
        str = "sqrt(pow((x[0]-0.5),2) + pow((x[1]-0.5),2) + pow((x[2]-%g),2))" % (z0)

    fP2 = Expression(str, degree=2) # define radial function on quadratic elements
    f = interpolate(fP2,Q) # just interpolate

    fx = f.vector()
    # scale the vector by h_on_delta
    fx *= h_on_delta

    n=int(n) # cast as integers
    m=int(m)
    # use the sinc solution to compute the solitary wave porosity ff at the collocation points rr
    if n==3 and m==0:
        rr, ff = solwave_mck(c, wd, numcol)
    elif n==2 and m==1:
        rr, ff = solwave_con(c, wd, numcol, verbose = False)

    #get numpy array for radial function and interpolate using sinc
    rmesh = fx[:]; # get numpy array
    phi = sincinterp_e(rr, ff - 1.0, rmesh) + 1.0

    for i in range(fx.size()):
        fx[i] = phi[i]
    return f
```

Define buoyancy-driven compaction problem

```
def compaction(V, h_n, n):
    # Test and Trial function
    P_n = TrialFunction(V)
    P_test = TestFunction(V)
    k_n = pow(h_n, n)

    # define linear and bilinear form
    aP = (P_n * P_test + k_n * inner(grad(P_n), grad(P_test))) * dx
    LP = k_n * Dx(P_test, 1) * dx - k_n * P_test * dot(z, nn) * ds

    return (aP, LP)
```

Define Time-dependent  problem (\partial \phi / \partial t = P)

```
def react(V, h_n, P_n, dt):
    h = TrialFunction(V)      # porosity \phi
    h_test = TestFunction(V)  # test function

    ah = h * h_test * dx
    Lh = (h_n + P_n * dt) * h_test * dx

    return (ah, Lh)
```

Define the main and parameters,

```
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    N = 64             # Number of grid points in x & y
    n = 3              # Coefficient in permeability relationship
    m = 0              # Coefficient in porosity relationship
    h_on_delta = 120.0  # Number of compaction lengths per box side
    courant = 1.0      # prefactor in time step criterion 
    
    mesh = UnitSquareMesh(N,N,"left/right")   # use a symmetric mesh    
    nn = FacetNormal(mesh)  # n vector is the facet normal
    V = FunctionSpace(mesh, "CG", 1)   # space for porosity
    z = Constant((0.0, 1.0))  # unit vector z

    c = 10         # solitary wave speed
    d = 2         # space dimension
    wd = 2        # wave dimension
    z0 = 0.5      # fractional position of wave in box
    numcol = 200  # number of collocation points
    shift_vel = c

    T = 50  # total simulation time
    t_num = 1000  # number of time steps
    dt = T / t_num  # time stepping

    # Single solitary wave initial condition -- uses Gideon's code
    f0 = dolfinsolwave(V,n,m,h_on_delta,c,d,wd,z0,numcol)
    h_n = project(f0, V)
```

It is important to rescale the mesh to be in number of compaction lengths

```
    # Expand mesh - Work in units scaled to the compaction length
    for x in mesh.coordinates():
        x[0] *= h_on_delta
        x[1] *= h_on_delta
```

Define the form of equations and time loop

```
    P_n = Function(V)
    h = Function(V)

    # define variation form
    aP, LP = compaction(V, h_n, n)
    ah, Lh = react(V, h_n, P_n, dt)

    plt.figure()
    Pfile = File("output/Pressure_testP.pvd")
    hfile = File("output/porosity_testP.pvd")

    h.assign(h_n)

    # main loop
    t = 0.0
    for n in range(t_num):
        # Solve variational problem
        problem1 = LinearVariationalProblem(aP, LP, P_n)
        solver1 = LinearVariationalSolver(problem1)
        solver1.solve()

        # Save file
        Pfile << P_n
        hfile << h

        problem2 = LinearVariationalProblem(ah, Lh, h)
        solver2 = LinearVariationalSolver(problem2)
        solver2.solve()

        # update current time
        t += dt

        # Update previous solution
        h_n.assign(h)
```





## Some thoughts from this benchmarks

1.The solitary wave code is difficult to find online. I don't think they published their code with their article, so I couldn't find any. However, the worst mistake I made was not talking to my supervisor about it in a timely manner. I couldn’t do my project and I‘m at a loss on what to do. So it wasted me loads of time.  **LESSON 1: CONTACT YOUR SUPERVISOR IF YOU HAVE A PROBLEM**.

2.When comparing two different codes, make sure that our parameters and the parameters they use are indicating the same physical property (eg, symbol g only represent gravity not anything else), and that they are consistent through the entire code (e.g, not a=1 at first and a=2 at the following code). Make sure that no parameter is defined more than once. **LESSON 2: Carefully check every parameters**.

3.The initial code I used was scaled by the height H, but the solitary wave code was scaled using the comapction length $$\delta$$, thus we must scale it with the same characteristic scale. Also remember to scale the mesh. **LESSON 3: Scale the equations in a same way**.

4.The solitary wave should be small enough to avoid the effect of the boundary. 

5.Since we’re working on a finite domain and a discrete mesh, the solitary wave will never be perfect. When we increase the mesh resolution and decrease the time step, the differences should to be less apparent.

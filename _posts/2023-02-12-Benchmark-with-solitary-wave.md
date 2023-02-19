---
title: 'Benchmark with solitary wave'
date: 2023-02-06
permalink:
tags:
  - Fluid dynamics
  - solitary wave
---



It is important to make sure our code is doing the right thing. A good practice is to benchmark it with the solitary wave code. 

## Governing equations for solitary wave

The governing equations is conservation of mass and Darcy’s law,

And If we assume no shear, no melting (reaction), the eqs become,

$$
-{\partial \phi \over \partial t} + \nabla \cdot v_s-\nabla \cdot (\phi v_s)=0
$$

$$
{\partial \phi \over \partial t} + \nabla \cdot (\phi v_l)=0
$$

$$
\phi (v_l-v_s) = K ((1-\phi) \Delta\rho g\hat z -\nabla P)
$$

These equations could be dimensionless in a different way compared to the previous post.

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

```
```

half done…

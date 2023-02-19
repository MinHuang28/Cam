---
title: 'FEniCS-weak form of RII problem PDE'
date: 2023-02-09
permalink:
tags:
  - FEniCS
  - PDE
---



Create the weak form of RII governing equations for FEniCS project

## Governing equations for RII problems

$$
{\partial \phi \over \partial t} = P + \chi
$$

$$
MP+\nabla \cdot q =0
$$

$$
q=\phi v_l=K(\widehat z-S \nabla P)
$$

$$
q \cdot ({\nabla \chi \over Da} - \hat z)= {1\over DaPe} {\partial \over \partial x}(\phi{\partial \chi \over \partial x})-\chi
$$



## Step_1: Solve Buoyancy-driven compaction (P)

Fluid flow is given by Darcy's law (3), and the mass conservation for fluid is given by (2). We obtain the **Strong form** of the equation by combining (2) and (3),

$$
MP+\nabla \cdot q = MP+\nabla \cdot [K(\widehat z-S \nabla P)]=0
$$

Multiply the whole equation by a **test function** P_test, and integrate the whole equation over the domain $$\Omega$$, 

$$
\int_\Omega MP P_1+\int_\Omega P_1 {\nabla \cdot [K(\widehat z-S\nabla P)]}=0
$$

the product rule of differentiation implies that,

$$
P_1 (\nabla \cdot [K(\widehat z - S\nabla P)] = \nabla \cdot(P_1 K(\widehat z - S\nabla P)) - \nabla P_1 \cdot  K(\widehat z - S\nabla P)
$$

So the Eq (6) becomes, 

$$
\int_\Omega MP P_1 + \int_\Omega [\nabla \cdot(P_1 K(\widehat z - S\nabla P)) - \nabla P_1 \cdot  K(\widehat z - S\nabla P) ]=0
$$

Apply the divergence theorem ($ \iiint_V (\nabla \cdot {P})dV = \iint_S (P \cdot \widehat n)dS $) to generate boundary integrals, surface normal $ \widehat n $, surface domain $ \Gamma = \partial \Omega $, 

$$
\int_\Omega MP P_1 + \int_\Gamma P_1 K(\widehat z-S\nabla P) \cdot n - \int_\Omega \nabla P_1 \cdot K(\widehat z -S\nabla P) =0
$$

The weak form is 

$$
\int_\Omega MP P_1 -\int_\Gamma P_1 K S \nabla P\cdot n + \int_\Omega KS \nabla P_1 \cdot \nabla P = \int_\Omega K {\partial P_1 \over \partial z} -\int_\Gamma P_1 K(\widehat z\cdot n)
$$

If we apply the boundary condition $${\partial P \over \partial n} =0$$, the second term would vanish.


## Step_2: Solve reaction-infiltration($$\chi$$)

With Darcy flux $$q = K(\widehat z-S \nabla P)$$, the equation becomes,

$$
q \cdot ({\nabla \chi \over Da} - \hat z)= {1\over DaPe} {\partial \over \partial x}(\phi{\partial \chi \over \partial x})-\chi
$$

Multiply the whole equation by a **test function** C_test, and integrate the whole equation over the domain $$\Omega$$. 

$$
\int_\Omega q \cdot ({\nabla \chi \over Da} - \hat z) C_1 - \int_\Omega[{1\over DaPe} {\partial \over \partial x}(\phi{\partial \chi \over \partial x})-\chi] C_1 =0
$$

$$
\int_\Omega q\cdot {\nabla \chi \over Da} C_1 - \int_\Omega q \cdot \hat z C_1 - \int_\Omega{1 \over DaPe} {\partial \over \partial x} (\phi{\partial \chi \over \partial x}) C_1+\int_\Omega \chi C_1 =0
$$

Apply the product rule,

$$
{\partial \over \partial x} (\phi{\partial \chi \over \partial x}) C_1 = \phi{\partial \over \partial x}({\partial \chi \over \partial x})C_1+{\partial \phi \over \partial x}{\partial \chi \over \partial x} C_1
$$

Integrating by parts,

$$
\int_\Omega {\partial \over \partial x}({\partial \chi \over \partial x})C_1 = [{\partial \chi \over \partial x} C_1]_{\partial \Omega} -\int_\Omega {\partial \chi \over \partial x} {\partial C_1 \over \partial x}
$$

The eq(13) becomes,

$$
\int_\Omega q\cdot {\nabla \chi \over Da} C_1 - \int_\Omega q \cdot \hat z C_1 - \int_\Omega{1 \over DaPe} [( \phi{\partial \chi \over \partial x} C_1)_{\partial \Omega} - \phi {\partial \chi \over \partial x} {\partial C_1 \over \partial x}+{\partial \phi \over \partial x}{\partial \chi \over \partial x} C_1]+\int_\Omega \chi C_1 =0
$$

So the weak form is, 

$$
\int_\Omega q \cdot {\nabla \chi \over Da} C_1 + \int_\Omega{1 \over DaPe} [\phi {\partial \chi \over \partial x} {\partial C_1 \over \partial x} -{\partial \phi \over \partial x}{\partial \chi \over \partial x} C_1]+\int_\Omega \chi C_1 -{1 \over DaPe} [(\phi {\partial \chi \over \partial x} C_1)_{\partial \Omega}  = \int_\Omega q \cdot \hat z C_1
$$

## Step_3: Upgrade porosity with time

Governing equation is 

$$
{\partial \phi \over \partial t} = P + \chi
$$

Discretize the time derivative by a finite difference approximation

$$
{\partial\over \partial t} \phi ^n = (P+\chi) ^n
$$

The time-derivative can be approximated by a difference quotient, here we use a backward difference:

$$
{\partial\over \partial t} \phi ^n \approx  {\phi ^{n+1}- \phi^n \over \Delta t}
$$

Insert (20) in (19) yields

$$
{\phi ^{n+1}- \phi^n \over \Delta t} =  (P+ \chi)^n
$$

So we got a discretization.

$$
\phi ^0 =\phi_0,\\
\phi ^{n+1}- \phi^n =  (P+ \chi)^n \Delta t
$$

Multiply the whole equation by a **test function** , and integrate the whole equation over the domain $$\Omega$$, Introducing the symbol $$\phi$$ for $$\phi^{n+1}$$, the weak form can be written,

$$
a(\phi,R)=\int_\Omega \phi R dx \\
L_n(R) = \int_\Omega [\phi^n +(P+\chi)^n\Delta t] R dx
$$

Initial condition $$a_0(\phi, R) = L_0(R)$$

$$
a_0(\phi, R) = \int_\Omega \phi R dx,\\
L_0(R) =\int_\Omega \phi_0 R dx
$$

## Domain and boundary condition

Domain is a rectangle area, with 1.0 in height and 2.0* $$\pi$$/10 in length

Boundary condition is 

$$
\phi =1, \chi =1, {\partial P\over \partial z} =0,  (z =0)\\
{\partial P\over \partial z} =0,  (z =1)
$$





<script type="text/javascript" src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

---
title: 'FEniCS-Variation formulation of RII problem'
date: 2023-02-08
permalink:
tags:
  - FEniCS
  - PDE
---

Create the variation form of RII governing equations for FEniCS project

Governing equations

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

Fluid flow is given by Darcy's law, Dimensionless equation: (Jones, 2018)

$$
q=\phi v_l=K(\widehat z-S \nabla P)
$$

$$
MP+\nabla \cdot q =0
$$

Introduced $$M= \beta H/\alpha$$,  $$S= M\delta^2/ H^2 =\beta \zeta K_0/\alpha H$$.

The **Strong form** of the equation is,

$$
MP+\nabla \cdot q = MP+\nabla \cdot [K(\widehat z-S \nabla P)]=0
$$

Weak formulation:

Multiply the whole equation by a **test function** P_test, and integrate the whole equation over the domain $$\Omega$$. 

$$
\int_\Omega MP P_1+\int_\Omega P_1 {\nabla \cdot [K(\widehat z-S\nabla P)]}=0
$$

Apply the boundary condition $$\frac{\partial P}{\partial n} =0$$, on $$\Gamma$$, we get the weak form,

$$
\int_\Omega MP P_1+\int_\Gamma P_1 K(\widehat z-S\nabla P) \cdot n - \int_\Omega \nabla P_test \cdot K(\widehat z -S\nabla P) =0
$$

$$
\int_\Omega MP P_1 + \int_\Omega KS \nabla P_1 \cdot \nabla P = \int_\Omega K {\partial P_1 \over \partial z} -\int_\Gamma P_1 K(\widehat z\cdot n)
$$

## Step_2: Solve reaction-infiltration(\chi)

With Darcy flux $$q = K(\widehat z-S \nabla P)$$, the equations becomes,

$$
q \cdot ({\nabla \chi \over Da} - \hat z)= {1\over DaPe} {\partial \over \partial x}(\phi{\partial \chi \over \partial x})-\chi
$$

Multiply the whole equation by a **test function** C_test, and integrate the whole equation over the domain $$\Omega$$. 

$$
\int_\Omega q \cdot ({\nabla \chi \over Da} - \hat z) C_1 - \int_\Omega[{1\over DaPe} {\partial \over \partial x}(\phi{\partial \chi \over \partial x})-\chi] C_1 =0
$$

$$
\int_\Omega q\cdot {\nabla \chi \over Da} C_1 - \int_\Omega q \cdot \hat z C_1 - \int_\Omega{\phi \over DaPe} {\partial \over \partial x} ({\partial \chi \over \partial x}) C_1+\int_\Omega \chi C_1 =0
$$

Integrating by parts,

$$
\int_\Omega {\partial \over \partial x}({\partial \chi \over \partial x})C_1 = [{\partial \chi \over \partial x} C_1]_{\partial \Omega} -\int_\Omega {\partial \chi \over \partial x} {\partial C_1 \over \partial x}
$$

So the weak form is, 

$$
\int_\Omega q \cdot {\nabla \chi \over Da} C_1 - \int_\Omega q \cdot \hat z C_1 + \int_\Omega{\phi \over DaPe} {\partial C_1 \over \partial x} {\partial \chi \over \partial x}+\int_\Omega \chi C_1 =0
$$

$$
\int_\Omega q \cdot {\nabla \chi \over Da} C_1+ \int_\Omega{\phi \over DaPe} {\partial C_1 \over \partial x} {\partial \chi \over \partial x}+\int_\Omega \chi C_1 = \int_\Omega q \cdot \hat z C_1
$$

## Step3: Upgrade porosity with time

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

Insert (14) in (13 ) yields

$$
{\phi ^{n+1}- \phi^n \over \Delta t} =  (P+ \chi)^n
$$

So we got a *forward Euler* discretization.

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

Domain is a rectangle area

Boundary condition is

$$
\phi =1, \chi =1, {\partial P\over \partial z} =0,  (z =0)\\
{\partial P\over \partial z} =0,  (z =1)
$$





<script type="text/javascript" src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

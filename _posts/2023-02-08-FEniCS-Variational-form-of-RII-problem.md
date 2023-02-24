---
title: 'FEniCS-Weak form of RII problem'
date: 2023-02-09
permalink:
tags:
  - FEniCS
  - PDE
---



Create the weak form of RII governing equations for FEniCS project

As we previously discussed, the governing equations for RII problems are as follows,

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

As we can see, the strong form PDEs (1)-(4) are of the second order. 

This requires a high level of smoothness in the solution. 

That means the second derivative of the displacement must exist and be continuous! 

This also implies needs for parameters that are not influenceable such as geometry (sharp edges) and material parameters.

To develop the finite element formulation, the PDEs must be restated in an integral form known as the **weak form**. The weak form and the strong form are **equivalent**.

The Steps of creating the weak form would be sth like this,

1. Multiply the governing equation and the traction boundary condition by an **arbitrary** function (Called **test function**) and then integrate over the domains

   Note: Arbitrariness of the test function is crucial for the weak form. Otherwise the strong form is **NOT** equivalent to the weak form.

2. Transform eqs into a form containing only first derivatives. (Usually applying the integration by parts)

Let’s start!


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

Apply the divergence theorem ($$ \iiint_V (\nabla \cdot {P})dV = \iint_S (P \cdot \widehat n)dS $$) to generate boundary integrals, surface normal $$ \widehat n $$, surface domain $$ \Gamma = \partial \Omega $$, 

$$
\int_\Omega MP P_1 + \int_\Gamma P_1 K(\widehat z-S\nabla P) \cdot n - \int_\Omega \nabla P_1 \cdot K(\widehat z -S\nabla P) =0
$$

The weak form is 

$$
\int_\Omega MP P_1 -\int_\Gamma P_1 K S \nabla P\cdot n + \int_\Omega KS \nabla P_1 \cdot \nabla P = \int_\Omega K {\partial P_1 \over \partial z} -\int_\Gamma P_1 K(\widehat z\cdot n)
$$

If we apply the boundary condition $${\partial P \over \partial n} =0$$ on $$\Gamma$$, the second term would vanish.


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

Integrating by parts,

$$
\int_\Omega {\partial \over \partial x} (\phi{\partial \chi \over \partial x}) C_1 = [\phi{\partial \chi \over \partial x} C_1]_{\partial \Omega}- \int_\Omega  \phi {\partial \chi \over \partial x}{\partial C_1 \over \partial x}
$$

The eq(13) becomes,

$$
\int_\Omega q\cdot {\nabla \chi \over Da} C_1 - \int_\Omega q \cdot \hat z C_1 - \int_\Omega{1 \over DaPe} [( \phi{\partial \chi \over \partial x} C_1)_{\partial \Omega} - \phi {\partial \chi \over \partial x} {\partial C_1 \over \partial x}]+\int_\Omega \chi C_1 =0
$$

So the weak form is, 

$$
\int_\Omega q \cdot {\nabla \chi \over Da} C_1 + \int_\Omega{1 \over DaPe} \phi {\partial \chi \over \partial x} {\partial C_1 \over \partial x}+\int_\Omega \chi C_1 -{1 \over DaPe} [(\phi {\partial \chi \over \partial x} C_1)_{\partial \Omega}  = \int_\Omega q \cdot \hat z C_1
$$

If we apply the boundary condition $${\partial \chi \over \partial n} =0$$ on $$\Gamma$$$, the boundary term would vanish?

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
P = - 1, \chi =1, {\partial P\over \partial z} =0,  (z =0)\\
{\partial P\over \partial z} =0,  (z =1)
$$





<script type="text/javascript" src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

---
title: 'Governing equations of RII problem'
date: 2023-02-06
permalink:
tags:
  - Fluid dynamics
  - PDE
---


### The full system of PDEs that govern the fluid mechanics

#### Governing equations

1. **Conservation of mass** in both phases

$$
{\partial (1-\phi) \over \partial t} + \nabla \cdot ((1-\phi)v_s)=-\Gamma
$$

$$
{\partial \phi \over \partial t} + \nabla \cdot (\phi v_l)=\Gamma
$$

where $$\phi$$ is the porosity, $$(1-\phi)$$ is the fraction of solid phase, t is time, $$v_l$$ is fluid velocity, and $$v_s$$ is solid velocity, $$\Gamma$$ is the melting rate. We assume that $$\Gamma$$ is proportional to the undersaturation of soluble component in the melt, so $$\Gamma = -R(c_l-c_l^{eq})$$, where R is a kinetic coefficient with units 1/times.

2. **Darcy’s law** (conservation of momentum for liquid)

$$
\phi (v_l-v_s) = K ((1-\phi) \Delta\rho g\hat z -\nabla P)
$$

A Darcy flux$$\phi(v_l-v_s)$$ is driven by gravity $$g\hat z$$ associated with the density difference $$\Delta \rho$$, and by compaction pressure gradient $$\nabla P$$. The mobility K, representing the permeability k divided by liquid viscosity $$\mu$$, depends on the porosity. $$K =K_0 (\phi/\phi_0)^n$$, where $$K_0$$ is a reference mobility at a reference porosity $$\phi_0$$, and n is a constant. 

3. **Chemical component conservation** equation, consisting of three contributions: diffusion, advection and chemical source term.

$$
{\partial \over \partial t}(\phi c_l) +\nabla \cdot (\phi v_l c_l) = \nabla \cdot (\phi D \nabla c_l) + \Gamma c_\Gamma
$$


where D is diffusivity in the liquid phase; diffusion through the solid phase is negligible, $$C_\Gamma$$, is the concentration of reactively produced melts. 

Expand the partial derivatives and simplify (4) with (2)b to get
$$
\phi{\partial c_l\over \partial t} +\phi v_l \cdot \nabla  c_l = \nabla \cdot (\phi D \nabla c_l) + (c_\Gamma - c_l)\Gamma
$$

#### Non-dimensional Scales

**Solubility** is assumed to be a linear function of height, $${\partial c_{eq} \over \partial z}=\beta$$, $$\beta$$ is the solubility gradient($$m^{-1}$$)

Assuming zero solubility at the base area z=0, $$c_{eq}= \beta z$$, 

the **concentration** $$c_\Gamma$$ is offset from the equilibrium concentration by $$\alpha$$, $$c_\Gamma=\beta z +\alpha$$,

**Characteristic scales**:

​				$$(x,z) = H (x’,z') $$, 

​				$$ \phi = \phi_0 \phi’$$,

​				$$v_l=w_0 v_l'=K_0 \Delta \rho g/\phi_0 v_l' $$,

​				$$v_s=\phi_0 w_0 v_s’$$,

​				$$t=\alpha/(w_0\beta)t’$$,

​				$$ P = P_0 P’ =\zeta \phi_0 w_0 \beta/\alpha  P'$$,

​				$$c_l=\beta H c_l’$$,

​				$$\Gamma = \phi_0 w_0 \beta / \alpha \Gamma’$$,



#### Non-dimensionalization step by step

##### 1. The dimensionless mass conservation in the solid

Starting with the first equation (1), $${\partial (1-\phi) \over \partial t} + \nabla \cdot ((1-\phi)v_s)=-\Gamma$$, 

$$
-{\partial \phi \over \partial t} + \nabla \cdot v_s-\nabla \cdot (\phi v_s)=-\Gamma
$$

The compaction rate $$\nabla \cdot v_s$$ is related to the compaction pressure P according to the linear constitutive law.

$$
\nabla \cdot v_s = P/\zeta
$$

so the eq (7) becomes,

$$
-{\partial \phi \over \partial t} + P/\zeta - \nabla \cdot (\phi v_s) = -\Gamma
$$

put the dimensionless parameter for the dimensional one,

$$
-{\partial (\phi_0 \phi) \over \partial (\alpha/(w_0\beta)t)} + \zeta \phi_0 w_0 \beta/\alpha P/\zeta - {1\over H} \nabla \cdot (\phi_0 \phi (\phi_0 w_0 v_s)) = -\phi_0 w_0 \beta / \alpha \Gamma
$$

$$
-{w_0\beta\phi_0 \over \alpha}{\partial \phi \over \partial t} + {\phi_0 w_0 \beta \over \alpha} P- {\phi_0^2  w_0\over H} \nabla \cdot ( \phi v_s) = -{\phi_0 w_0 \beta \over \alpha} \Gamma
$$

Multiple by $${ \alpha\over\phi_0 w_0 \beta }$$, we can get,

$$
-{\partial \phi \over \partial t} + P- {\phi_0 \alpha \over \beta H} \nabla \cdot ( \phi v_s) = -\Gamma
$$

The dimensionless reactive melting rate $$\Gamma$$ is equal to the scaled undersaturation $$\chi$$, 

$$
{\partial \phi \over \partial t} - P + {\phi_0 \alpha \over \beta H} \nabla \cdot ( \phi v_s) = \chi
$$

We simplify the equations by taking the limit of small porosity $$\phi_0 <<M << 1$$, and $$M = {\beta H \over \alpha}$$, so $${\phi_0 \alpha \over \beta  H} << 1$$, we can neglect the third term and get,

$$
{\partial \phi \over \partial t}= P + \chi
$$

##### 2. The dimensionless mass conservation in the liquid

The scaled eq(2) $${\partial \phi \over \partial t} + \nabla \cdot (\phi v_l)=\Gamma$$ becomes,

$$
{\partial \phi_0\phi \over \partial (\alpha/(w_0\beta)t)} + {1 \over H}\nabla \cdot (\phi_0 \phi w_0 v_l)= \phi_0 w_0 \beta / \alpha \Gamma
$$

$$
{\phi_0 w_0\beta \over \alpha}{\partial \phi \over \partial t} + {\phi_0 w_0 \over H}\nabla \cdot ( \phi v_l)= {\phi_0 w_0\beta \over \alpha} \Gamma
$$

Multiple by $${ \alpha\over\phi_0 w_0 \beta }$$ and get,

$$
{\partial \phi \over \partial t} + {\alpha \over \beta H}\nabla \cdot ( \phi v_l)= \chi
$$

$$M= {\beta H \over \alpha}$$,  so $${\alpha\over \beta H} = 1/M$$, 

$$
M {\partial \phi \over \partial t} + \nabla \cdot ( \phi v_l)= M\chi
$$

##### 3. The dimensionless Darcy’s law

eq (3) Darcy’s law, $$\phi (v_l-v_s) = K ((1-\phi) \Delta\rho g\hat z -\nabla P)$$, 

$$
\phi_0\phi (w_0v_l-\phi_0 w_0 v_s) = K_0 K ((1-\phi_0\phi) \Delta\rho g\hat z - {1 \over H} \nabla ((\zeta \phi_0 w_0 \beta/\alpha) P))
$$

Taking the limit of small porosity $$\phi_0 <<M << 1$$, we can neglect the $$\phi_0 w_0 v_s$$ and $$\phi_0\phi \Delta\rho g\hat z $$, 

$$
\phi_0 w_0 \phi v_l = K_0 K  (\Delta\rho g\hat z - {1 \over H} \nabla ({\zeta \phi_0 w_0 \beta \over\alpha} P)
$$

Divided by $$\phi_0 w_0$$, 

$$
\phi v_l = K_0 K  ({\Delta\rho g \over \phi_0 w_0 }\hat z - {\zeta \beta \over \alpha H} \nabla  P)
$$

We have $$w_0=K_0 \Delta \rho g/\phi_0 $$, so $$K_0 = \phi_0 w_0 / \Delta \rho g$$, 

$$
\phi v_l =K  (\hat z - {K_0 \zeta \beta \over \alpha H}  \nabla P)
$$

Stiffness $$S = M  {\delta^2 \over H^2} = {\beta H \over \alpha } {K_0 \zeta \over H^2} = {\beta K_0 \zeta \over \alpha H} $$, 

$$
\phi v_l =K  (\hat z - S \nabla P)
$$

##### 4. The dimensionless Chemical component conservation

eq (4) $$\phi{\partial c_l\over \partial t} +\phi v_l \cdot \nabla c_l = \nabla \cdot (\phi D \nabla c_l) + (c_\Gamma - c_l)\Gamma$$, 

Scaled eq(4),

$$
\phi_0\phi{\partial (\beta H c_l)\over \partial (\alpha/(w_0\beta) t)} +\phi_0w_0 \phi v_l  {1 \over H}\nabla \cdot (\beta H c_l) =  {1 \over H} \nabla \cdot (\phi_0 \phi D  {1 \over H} \nabla (\beta H c_l)) + (\beta z +\alpha - \beta z) \phi_0 w_0 \beta / \alpha \Gamma
$$

$$
{\phi_0 w_0\beta\beta H  \over \alpha}\phi{\partial c_l\over \partial t} +\phi_0w_0 
\beta \phi v_l\nabla \cdot c_l =  {\phi_0 \beta D  \over H} \nabla \cdot ( \phi \nabla c_l) + \phi_0 w_0 \beta\Gamma
$$

Given $$M= {\beta H \over \alpha} <<1$$, we can neglect the first term $$\phi_0 w_0\beta{\beta H  \over \alpha}\phi{\partial c_l\over \partial t} $$,  and divided the whole equation by $$\phi_0 w_0 \beta$$, 

$$
\phi v_l\cdot  \nabla c_l =  {D  \over w_0 H } \nabla \cdot ( \phi \nabla c_l) + \chi
$$

Undersaturation $$\chi = Da(z-c_l)$$, $$c_l = z - {\chi \over Da}$$, $$ \nabla c_l = \nabla(z- {\chi \over Da})=z-{\nabla\chi \over Da}$$

And $$\nabla \cdot (\phi \nabla c_l)={\partial \over \partial x} (\phi {\partial c_l\over \partial x}) + {\partial \over \partial z} (\phi {\partial c_l\over \partial z})$$, $$c_l$$ is meaningful in x axis, so $$\nabla \cdot (\phi \nabla c_l)={\partial \over \partial x} (\phi {\partial c_l\over \partial x}) =-{1\over Da} {\partial \over \partial x} (\phi {\partial \chi\over \partial x})$$

So the equation (29) becomes,

$$
\phi v_l\cdot  [z-{\nabla\chi \over Da}] =  -{1  \over Da Pe } {\partial \over \partial x} (\phi {\partial \chi\over \partial x}) + \chi
$$

$$
\phi v_l\cdot  [{\nabla\chi \over Da}-z] =  {1  \over Da Pe } {\partial \over \partial x} (\phi {\partial \chi\over \partial x}) - \chi
$$

#### Dimensionless equations

So the four dimensionless governing equations (1)-(4) become,
$$
{\partial \phi \over \partial t} = P + \chi
$$

$$
MP+\nabla \cdot q =M\chi
$$

$$
q=\phi v_l=K(\widehat z-S \nabla P)
$$

$$
q \cdot ({\nabla \chi \over Da} - \hat z)= {1\over DaPe} {\partial \over \partial x}(\phi{\partial \chi \over \partial x})-\chi
$$


<script type="text/javascript" src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

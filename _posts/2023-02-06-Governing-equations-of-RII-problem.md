---
title: 'Governing equations of RII problem'
date: 2023-02-06
permalink:
tags:
  - Fluid dynamics
  - PDE
---



The full system of PDEs that govern the fluid mechanics



## Governing equations

1.**Conservation of mass** in both phases

$$
{\partial (1-\phi) \over \partial t} + \nabla \cdot ((1-\phi) \mathbf{v_s})=-\Gamma
$$

$$
{\partial \phi \over \partial t} + \nabla \cdot (\phi \mathbf{v_l})=\Gamma
$$

where $$\phi$$ is the porosity, $$(1-\phi)$$ is the fraction of solid phase, t is time, $$v_l$$ is fluid velocity, and $$v_s$$ is solid velocity, $$\Gamma$$ is the melting rate. We assume that $$\Gamma$$ is proportional to the undersaturation of soluble component in the melt, so $$\Gamma = R A(\phi,c_s)(c_l^{eq}-c_l)$$, where R is a kinetic coefficient with units 1/times.

2.**Darcy’s law** (conservation of momentum for liquid)

$$
\phi (\mathbf{v_l}-\mathbf{v_s}) = {k \over \mu} ((1-\phi) \Delta\rho g \mathbf{\hat z} -\nabla \mathcal{P})
$$

A Darcy flux $$\phi (\mathbf{v_l}-\mathbf{v_s})$$ is driven by gravity $$g\hat \mathbf {z}$$ associated with the density difference $$\Delta \rho$$, and by compaction pressure gradient $$\nabla \mathcal{P}$$. The permeability k is divided by liquid viscosity $$\mu$$, depending on the porosity. $$k =k_0 (\phi/\phi_0)^n$$, where $$k_0$$ is a reference permeability at a reference porosity $$\phi_0$$, and n is a constant. 

3.**Chemical component conservation** equation, consisting of three contributions: diffusion, advection and chemical source term.

$$
{\partial \over \partial t} (\phi c_l)+  \nabla \cdot (\phi \mathbf{v_l} c_l) = \nabla \cdot (\phi D \nabla c_l) + c_\Gamma \Gamma
$$


where D is diffusivity in the liquid phase; diffusion through the solid phase is negligible, $$C_\Gamma$$, is the concentration of reactively produced melts. 

The solid conservation,

$$
{\partial \over \partial t} ((1-\phi) c_s) + \nabla \cdot ((1-\phi) \mathbf{v_s} c_s ) =- c_\Gamma \Gamma
$$



## Non-dimensional Scales

**Solubility** is assumed to be a linear function of height, $${\partial c_{eq} \over \partial z}=\beta$$, $$\beta$$ is the solubility gradient($$m^{-1}$$)

Assuming zero solubility at the base area z=0, $$c_{eq}= \beta z$$, 

The **concentration** $$c_\Gamma$$ is offset from the equilibrium concentration by $$\alpha$$, $$c_\Gamma=\beta z +\alpha$$,

Characteristic velocity,

$$
w_0 = {k_0 \Delta \rho g \over  \mu \phi_0}
$$

Other characteristic scales:

$$(x,z) = H (x’,z') $$,  $$\nabla$$= $${1\over H } \nabla$$,

$$ \phi = \phi_0 \phi’$$,

$k=k_0 (d')^2 (\phi') ^n$

$$\mathbf{v_l}=w_0 \mathbf{v_l}'=(k_0 \Delta \rho g/\mu \phi_0) \mathbf{v_l}' $$,

$$\mathbf{v_s}=\phi_0 w_0 \mathbf{v_s}’$$,

$$t=\alpha/(w_0\beta)t’$$,                   

$$ \mathcal{P} = \mathcal{P_0} \mathcal{P}’ =(\zeta_0 \phi_0 w_0 \beta/\alpha) \mathcal{P}'$$,

$$c_l=\beta H c_l’$$,

$$\Gamma = (\phi_0 w_0 \beta / \alpha) \Gamma’$$,

$c_s = c_{s0} c_s’$,

$d = d_0 d'=bc_{s0} d'$, 

$$\zeta=\zeta_0 \zeta’$$,

**Scaled reaction rate**

Reaction rate $$\Gamma = R A(\phi,c_s)(c_l^{eq}-c_l)$$,

$$
{\phi_0 w_0 \beta \over \alpha}\Gamma = R A'(\beta H z-\beta H c_l)\\
$$

So the scaled Reaction rate $$\Gamma’$$ is,

$$
\Gamma={R \alpha \over \phi_0 w_0 \beta} A' (\beta H z-\beta H c_l)\\
\Gamma={\alpha R H \over \phi_0 w_0} A'(z-c_l) \\
\Gamma=Da A'(z-c_l), A’= {c_s (1-\phi_0 \phi) \over c_{s0}(1-\phi_0)}\\
$$



## Non-dimensionalization step by step

### 1. The dimensionless mass conservation in the solid

Starting with the first equation (1), $${\partial (1-\phi) \over \partial t} + \nabla \cdot ((1-\phi) \mathbf{v_s)}=-\Gamma$$, 

$$
-{\partial \phi \over \partial t} + (1- \phi ) \nabla \cdot \mathbf{v_s} - \mathbf{v_s}\nabla \phi =-\Gamma
$$

The compaction rate $$\nabla \cdot \mathbf{v_s}$$ is related to the compaction pressure $\mathcal{P}$ according to the linear constitutive law.

$$
\nabla \cdot \mathbf{v_s} = \mathcal{P} / \zeta = \mathcal{C}
$$

Eq (1) becomes,

$$
-{\partial \phi \over \partial t} + (1-\phi) P/\zeta - \mathbf{v_s} \nabla \cdot \phi = -\Gamma
$$

Put the dimensionless parameter for the dimensional one,

$$
-{\partial (\phi_0 \phi) \over \partial (\alpha/(w_0\beta)t)} + (1-\phi_0 \phi){\zeta_0 \phi_0 w_0 \beta/\alpha P \over \zeta_0 \zeta} - {\phi_0 w_0 \mathbf{v_s} \over H} \cdot\nabla  (\phi_0 \phi ) = -{\phi_0 w_0 \beta \over \alpha} \Gamma
$$

$$
-{w_0\beta\phi_0 \over \alpha}{\partial \phi \over \partial t} + (1-\phi_0 \phi){\phi_0 w_0 \beta \over \alpha} {P \over \zeta}- {\phi_0^2  w_0  \mathbf{v_s}\over H}\cdot \nabla  \phi = -{\phi_0 w_0 \beta \over \alpha} \Gamma
$$

Multiple by $${ \alpha\over\phi_0 w_0 \beta }$$, and get,

$$
-{\partial \phi \over \partial t} + (1-\phi_0 \phi) C - {\phi_0 \over M} \mathbf{v_s} \cdot \nabla \phi= -\Gamma
$$

$$
{\partial \phi \over \partial t} +{\phi_0 \over M} \mathbf{v_s} \cdot \nabla \phi = (1- \phi_0 \phi) C +\Gamma
$$

The dimensionless reactive melting rate $$\Gamma$$ is equal to the scaled undersaturation $$\chi$$, 

We simplify the equations by taking the limit of small porosity $$\phi_0 <<M << 1$$, and $$M = {\beta H \over \alpha}$$, so $${\phi_0 \over M} << 1$$. We also assume $\zeta$ is a constant, the equation becomes,

$$
{\partial \phi \over \partial t}= \mathcal{P} + \chi
$$



### 2. The dimensionless mass conservation in the liquid

The scaled eq(2) $${\partial \phi \over \partial t} + \nabla \cdot (\phi \mathbf{v_l})=\Gamma$$ becomes,

$$
{\partial (\phi_0\phi) \over \partial (\alpha/(w_0\beta)t)} + {1 \over H}\nabla \cdot (\phi_0 \phi w_0 \mathbf{v_l})= \phi_0 w_0 \beta / \alpha \Gamma
$$

$$
{\phi_0 w_0\beta \over \alpha}{\partial \phi \over \partial t} + {\phi_0 w_0 \over H}\nabla \cdot ( \phi \mathbf{v_l})= {\phi_0 w_0\beta \over \alpha} \Gamma
$$

Multiple by $${ \alpha\over\phi_0 w_0 \beta }$$ and get,

$$
{\partial \phi \over \partial t} + {\alpha \over \beta H}\nabla \cdot ( \phi \mathbf{v_l})= \Gamma
$$

$$M= {\beta H \over \alpha}$$,  so $${\alpha\over \beta H} = 1/M$$, 

$$
M {\partial \phi \over \partial t} + \nabla \cdot ( \phi \mathbf {v_l})= M\chi
$$



### 3. The dimensionless Darcy’s law

eq (3) Darcy’s law, $$\phi (\mathbf{v_l}-\mathbf{v_s}) = {k \over \mu} ((1-\phi) \Delta\rho g\hat z -\nabla \mathcal{P})$$, 

$$
\phi_0\phi (w_0v_l-\phi_0 w_0 v_s) = K_0 K ((1-\phi_0\phi) \Delta\rho g\hat z - {1 \over H} \nabla ((\zeta \phi_0 w_0 \beta/\alpha) \mathcal{P}))
$$

Taking the limit of small porosity $$\phi_0 <<M << 1$$, we can neglect the $$\phi_0 w_0 v_s$$ and $$\phi_0\phi \Delta\rho g\hat z $$, 

$$
\phi_0 w_0 \phi v_l = K_0 K  (\Delta\rho g\hat z - {1 \over H} \nabla ({\zeta \phi_0 w_0 \beta \over\alpha} P)
$$

Divided by $$\phi_0 w_0$$, 

$$
\phi(\mathbf{v_l} -\phi_0 \mathbf{v_s}) ={k_0 k \over \mu} ({1\over \phi_0 w_0 }(1- \phi_0\phi) \Delta \rho g \mathbf{\hat{z}} - {1 \over H} \nabla ({\zeta_0 \beta \over\alpha} \mathcal{P}))
$$

Darcy’s flux $$\phi_0 w_0= {k_0 \Delta \rho g \over \mu}$$, so $$k_0 = {\phi_0 w_0\mu \over \Delta \rho g}$$, 

$$
\phi(\mathbf{v_l} -\phi_0 \mathbf{v_s}) = k((1-\phi_0 \phi) \mathbf{\hat z}- {k_0 \zeta_0 \beta  \over\alpha H} \nabla \mathcal{P})
$$

Stiffness $$S = M  {\delta^2 \over H^2} = {\beta H \over \alpha } {K_0 \zeta \over H^2} = {\beta k_0 \zeta \over \alpha H} $$, 

$$
\phi(\mathbf{v_l} -\phi_0 \mathbf{v_s}) = k((1-\phi_0 \phi)\mathbf{\hat z}- S\nabla \mathcal{P})
$$

Taking the limit of small porosity $$\phi_0 <<M << 1$$, we can neglect the $$\phi \phi_0$$ term and get,

$$
\phi \mathbf {v_l} = k (\hat z - S \nabla \mathcal{P})
$$



### 4. The dimensionless Chemical component conservation

Eq (5)  $${\partial \over \partial t} (\phi c_l)+  \nabla \cdot (\phi \mathbf{v_l} c_l) = \nabla \cdot (\phi D \nabla c_l) + c_\Gamma \Gamma$$

$$
\phi{\partial c_l\over \partial t} +c_l{\partial \phi\over \partial t}+  \phi \mathbf{v_l} \nabla \cdot c_l +  c_l \nabla \cdot (\phi \mathbf{v_l})  = \nabla \cdot (\phi D \nabla c_l) + c_\Gamma \Gamma
$$

Use eq(2) $${\partial \phi \over \partial t} + \nabla \cdot (\phi \mathbf{v_l})=\Gamma$$ to simply eq(28)

$$
\phi{\partial c_l\over \partial t} +\phi \mathbf{v_l} \cdot \nabla c_l = \nabla \cdot (\phi D \nabla c_l) + (c_\Gamma - c_l)\Gamma
$$

Scaled eq(4),

$$
\phi_0\phi{\partial (\beta H c_l)\over \partial (\alpha/(w_0\beta) t)} +\phi_0w_0 \phi \mathbf{v_l} \cdot {1 \over H}\nabla  (\beta H c_l) =  {1 \over H} \nabla \cdot (\phi_0 \phi D  {1 \over H} \nabla (\beta H c_l)) + \alpha \phi_0 w_0 \beta / \alpha \Gamma
$$

$$
{\phi_0 w_0\beta\beta H  \over \alpha}\phi{\partial c_l\over \partial t} +\phi_0w_0 
\beta \phi \mathbf{v_l} \cdot \nabla c_l =  {\phi_0 \beta D  \over H} \nabla \cdot ( \phi \nabla c_l) + \phi_0 w_0 \beta \Gamma
$$

Divide the whole equation by $$\phi_0 w_0 \beta$$, and $$M= {\beta H \over \alpha}$$,  

$$
M\phi{\partial c_l\over \partial t} +\phi \mathbf{v_l} \cdot \nabla  c_l =  {D \over w_0 H} \nabla \cdot ( \phi \nabla c_l) +\Gamma
$$

Since $$Pe = {w_o H \over D}$$, eq(28) becomes,

$$
\phi M{\partial c_l\over \partial t} +\phi \mathbf{v_l} \cdot \nabla  c_l =  {1 \over Pe} \nabla \cdot ( \phi \nabla c_l) +\Gamma
$$

Given $$M= {\beta H \over \alpha} <<1$$, we can neglect the first term. Given undersaturation $$\chi = Da(z-c_l)$$, $$c_l = z - {\chi \over Da}$$,

$$
\nabla c_l = \nabla(z- {\chi \over Da})=\hat {z}-{\nabla\chi \over Da}
$$

$$\nabla \cdot (\phi \nabla c_l)={\partial \over \partial x} (\phi {\partial c_l\over \partial x}) + {\partial \over \partial z} (\phi {\partial c_l\over \partial z})$$,but $$c_l$$ is meaningful in x axis,

$$
\nabla \cdot (\phi \nabla c_l)={\partial \over \partial x} (\phi {\partial c_l\over \partial x}) =-{1\over Da} {\partial \over \partial x} (\phi {\partial \chi\over \partial x})
$$

So the equation becomes,

$$
\phi \mathbf{v_l} \cdot  [\hat z-{\nabla\chi \over Da}] =  -{1  \over Da Pe } {\partial \over \partial x} (\phi {\partial \chi\over \partial x}) + \chi
$$

$$
\phi \mathbf {v_l} \cdot  [{\nabla\chi \over Da}- \hat z] =  {1  \over Da Pe } {\partial \over \partial x} (\phi {\partial \chi\over \partial x}) - \chi
$$



### 5. The dimensionless Chemical component conservation

eq (5) $${\partial \over \partial t} ((1-\phi) c_s) + \nabla \cdot ((1-\phi) \mathbf{v_s} c_s ) =- c_\Gamma \Gamma$$, 

$$
c_s{\partial (1-\phi) \over \partial t} +(1-\phi) {\partial c_s \over \partial t} + (1-\phi) \mathbf{v_s} \cdot \nabla  c_s + c_s \nabla \cdot ((1-\phi) \mathbf{v_s}) =- c_\Gamma \Gamma
$$

Use eq(1) $${\partial (1-\phi) \over \partial t} + \nabla \cdot ((1-\phi) \mathbf{v_s})=-\Gamma$$ to simply eq(35),

$$
(1-\phi) {\partial c_s \over \partial t} + (1-\phi) \mathbf{v_s} \cdot \nabla c_s =(c_s- c_\Gamma) \Gamma
$$

Scale eq(36),

$$
(1-\phi_0\phi) {\partial c_{s0} c_s \over \partial (\alpha /w_0 \beta)t} + (1-\phi_0 \phi) \phi_0 w_0 {c_{s0} \over H} \mathbf{v_s} \cdot \nabla c_s =(c_{s0} c_s- c_\Gamma') (\phi_0 w_0 \beta /\alpha) \Gamma
$$

$$
{w_0 \beta c_{s0} \over \alpha}(1-\phi_0\phi) {\partial c_s \over \partial t} + {\phi_0 w_0 c_{s0} (1-\phi_0 \phi) \over H} \mathbf{v_s} \cdot \nabla c_s ={\phi_0 w_0 \beta \over \alpha}(c_{s0} c_s- c_\Gamma') \Gamma
$$

Divided by $$w_0 \beta /\alpha $$, 

$$
c_{s0}(1-\phi_0\phi) {\partial c_s \over \partial t} + {\phi_0 c_{s0}\alpha \over\beta H} (1-\phi_0 \phi) \mathbf{v_s} \cdot \nabla c_s ={\phi_0}(c_{s0} c_s- c_\Gamma') \Gamma
$$

$$
{\partial c_s \over \partial t} +{\phi_0 \over M} \mathbf{v_s}\cdot \nabla c_s = {\phi_0 \over (1-\phi_0 \phi) c_{s0}} (c_{s0} c_s-c_\Gamma') \Gamma
$$

If we take the limit $$\phi_0 <<M<<1$$, 

$$
{\partial c_s \over \partial t} = {\phi_0 \over (1-\phi_0 \phi) c_{s0}} (c_{s0} c_s-\beta H (z-{\chi \over Da})-\alpha) \chi
$$

## Dimensionless equations
So the five governing equations are,

$$
{\partial \phi \over \partial t} + {\phi_0 \over M} \mathbf{v_s} \cdot \nabla \phi= (1-\phi_0 \phi) C+\chi
$$

$$
M (1- \phi_0 \phi) C +  \phi \cdot \nabla v_l+ (v_l - \phi_0 v_s) \cdot \nabla \phi = 0
$$

$$
\phi(\mathbf{v_l} -\phi_0 \mathbf{v_s}) = k((1-\phi_0 \phi)\mathbf{\hat z}- S\nabla \mathcal{P})
$$

$$
{M \over Da}\phi{\partial \chi\over \partial t} +\phi v_l \cdot [{\nabla\chi \over Da}-z] =  {1 \over Pe} \nabla \cdot (\phi({\nabla \chi \over Da}-z)) - \chi
$$

$$
{\partial c_s \over \partial t} +{\phi_0 \over M} \mathbf{v_s}\cdot \nabla c_s = {\phi_0 \over (1-\phi_0 \phi) c_{s0}} (c_{s0} c_s-(\beta H + \alpha) \Gamma
$$


Taking the limit of small porosity, the five dimensionless governing equations (1)-(4) become,

$$
{\partial \phi \over \partial t} = \mathcal{P} + \chi
$$

$$
M{\partial \phi \over \partial t}+\nabla \cdot (\phi \mathbf {v_l}) =M\chi
$$

$$
\phi \mathbf{v_l}=k (\hat z-S \nabla \mathcal P)
$$

$$
\phi \mathbf {v_l} \cdot ({\nabla \chi \over Da} - \hat z)= {1\over DaPe} {\partial \over \partial x}(\phi{\partial \chi \over \partial x})-\chi
$$

$$
{\partial c_s \over \partial t} = {\phi_0 \over (1-\phi_0 \phi) c_{s0}} (c_{s0} c_s-\beta H (z-{\chi \over Da})-\alpha) \chi
$$




<script type="text/javascript" src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

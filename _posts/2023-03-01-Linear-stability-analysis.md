---
title: 'Linear stability analysis'
date: 2023-03-01
permalink:
tags:
  - Fluid dynamics
  - PDE
---

We can also benchmark our result with linear stability analysis.

Linear stability analysis is a mathematical method used to study the stability of a system by examining small perturbations from an equilibrium or steady state. 

We expand the governing equations with a z-dependent term and an (x, z, t) dependent perturbation,

$$
\phi= \phi_0(z) +\phi_1(x,z,t)\\
P= P_0(z) +P_1(x,z,t)\\
\chi= \chi_0(z) +\chi_1(x,z,t)\\
v= w_0(z) +v_1(x,z,t)\\
$$

The base state of governing equations are

$$
0=P_0+\chi_0
$$

$$
{d \over dz}(\phi_0 w_0) = M \chi_0
$$

$$
\phi_0 w_0 = K_0 [1-S {dP_0 \over dz}]
$$

$$
\phi_0 w_0 [{1\over Da} {d \chi_0 \over dz}-1]=-\chi_0
$$

Since M <<1, we can get,

$$
-P_0 =\chi_0=\phi_0 =w_0 =1
$$



## Governing equations for perturbations

The governing equations for perturbations can be written as,

$$
{\partial \phi_1 \over \partial t} = P_1 + \chi_1
$$

$$
M {\partial \phi_1 \over \partial t} + \phi_0 \nabla \cdot v_l + w_0 {\partial \phi_1 \over \partial z} =M\chi_1
$$

$$
\phi_0 v_1 = -SK_0 \nabla P_1 + (n-1)w_0 \phi_1 z
$$

$$
(\phi_0 w_1 +\phi_1 w_0)  ({1 \over Da} {d\chi_0 \over dz} - 1) + {\phi_0 w_0 \over Da} {\partial \chi_1 \over dz} = {\phi_0 \over DaPe} {\partial^2 \chi_1\over \partial x}-\chi_1
$$

We eliminate $$\chi_1$$ using eq(7) and $$v_1$$ using eq(9). 

Eq (8) becomes,

$$
M (P_1+\chi_1) + \phi_0 \nabla \cdot ({-SK_0 \over \phi_0} \nabla P_1 + (n-1)w_0 {\phi_1 \over \phi_0} z) + w_0 {\partial \phi_1 \over \partial z} =M\chi_1
$$

$$
-SK_0\nabla P_1^2  + (n-1)w_0 {\partial \phi_1 \over \partial z} + w_0 {\partial \phi_1 \over \partial z} = -M P_1
$$

$$
-SK_0\nabla P_1^2  + nw_0 {\partial \phi_1 \over \partial z} = -M P_1
$$

Then we simplify eq(10) and get,

$$
(\phi_0 w_1 +\phi_1 w_0)  ({1 \over Da} {d\chi_0 \over dz} - 1) =(-{\phi_0 w_0 \over Da} {\partial \over  \partial z}+{\phi_0 \over DaPe} {\partial^2 \over \partial x}-1)({\partial \phi_1 \over \partial t} - P_1 )
$$

From the base state eq(5), we get $${1\over Da} {d \chi_0 \over dz}-1={-\chi_0 \over \phi_0 w_0 }$$. 

$$
(\phi_0 w_1 +\phi_1 w_0)  ({-\chi_0 \over \phi_0 w_0 }) =(-{\phi_0 w_0 \over Da} {\partial \over  \partial z}+{\phi_0 \over DaPe} {\partial^2 \over \partial x}-1)({\partial \phi_1 \over \partial t} - P_1 )
$$

$$w_1$$ is the vertical component of $$v_1$$, so $$\phi_0 w_1$$ =$$\phi_0 v_1 \hat{z}$$ = $$-SK_0  {\partial P_1 \over \partial z} + (n-1)w_0 \phi_1 $$ 

$$
\phi_0 w_1 + \phi_1 w_0 = -SK_0  {\partial P_1 \over \partial z} + (n-1)w_0 \phi_1+\phi_1 w_0 = SK_0  {\partial P_1 \over \partial z} + nw_0 \phi_1
$$

eq(15) becomes,

$$
(-SK_0  {\partial P_1 \over \partial z} + nw_0 \phi_1)  ({-\chi_0 \over \phi_0 w_0 }) =(-{\phi_0 w_0 \over Da} {\partial \over  \partial z}+ {\phi_0 \over DaPe} {\partial^2 \over \partial x}-1)({\partial \phi_1 \over \partial t} - P_1 )
$$



### Perturbation Pressure

we substitute the constant base state, neglect the $$O(M)$$ term, eliminate $$\chi_1$$, and obtain,

$$
[{1\over Da} \partial _{tz} - {1 \over Da Pe} \partial_{txx} + \partial_t -n] \nabla^2 P_1 = {n \over S}[({1\over Da} -S)\partial_z -{1 \over Da Pe} \partial_{xx} +1] \partial_z P_1
$$

We seek normal-mode solution $$P_1$$ is proportional to $$exp(\sigma t +ikx +m_jz)$$, where the $$\sigma$$ is the growth rate and the k is a horizontal wavenumber. 

$$
[{1\over Da} \sigma m  - {1 \over Da Pe} \sigma (ik)^2+ \sigma -n] (m^2-k^2)e^{\sigma t +ikx +m_jz} = {n \over S}[({1\over Da} -S)m -{1 \over Da Pe} (ik)^2 +1] m e^{\sigma t +ikx +m_jz}
$$

Eq(19) becomes,

$$
[{\sigma m\over Da}   + {\sigma k^2 \over Da Pe} + \sigma -n] (m^2-k^2) = {n \over S}[{m\over Da} -Sm +{k^2 \over Da Pe} +1] m
$$

We set $$\mathcal{K}=1+{k^2 \over Da Pe}$$, 

$$
[{\sigma m\over Da}   + \sigma \mathcal{K} - n] (m^2-k^2) = {n \over S}[{m\over Da} -Sm +\mathcal{K}] m
$$

Then we can get the characteristic polynomial,

$$
\begin{eqnarray}
\frac{\unicode[STIX]{x1D70E}}{Da}m^{3}+\left(\unicode[STIX]{x1D70E}{\mathcal{K}}-\frac{n}{Da{\mathcal{S}}}\right)m^{2}-\left(\frac{n{\mathcal{K}}}{{\mathcal{S}}}+\frac{\unicode[STIX]{x1D70E}}{Da}k^{2}\right)m+\left(n-\unicode[STIX]{x1D70E}{\mathcal{K}}\right)k^{2}=0,
\end{eqnarray}
$$

Eq(22) has three roots $$m_j$$ (*j*=1,2,3), and hence the compaction pressure perturbation will be given by

$$
P_1=\mathop{\sum }_{j=1}^{3} A_j\exp (\unicode[STIX]{x1D70E}t+\text{i}kx+m_{j}z)
$$

The coefficient  prefactors $$A_j$$ could be determined by the boundary conditions.

Proper boundary conditions are,

$$
P_1 =0, {\partial P_1 \over \partial z} =0,  (z =0)
$$

$$
{\partial P_1 \over \partial z} =0,  (z =1)
$$

With $$P_1 =0$$ at $$z=0$$, we get,

$$
A_1 e^{\sigma t +ikx}+A_2 e^{\sigma t +ikx}+A_3 e^{\sigma t +ikx}=0
$$

With $${\partial P_1 \over \partial z} =0$$, at $$z =0$$, 

$$
m_1A_1 e^{\sigma t +ikx}+m_2 A_2 e^{\sigma t +ikx}+m_3 A_3 e^{\sigma t +ikx}=0
$$

With  $${\partial P_1 \over \partial z} =0$$, at $$z =1$$, 

$$
m_1 A_1 e^{\sigma t +ikx+m_1}+m_2 A_2 e^{\sigma t +ikx+m_2}+m_2 A_2 e^{\sigma t +ikx+m_2}=0
$$

These can be explained in matrix as,

$$
\begin{eqnarray}
\left(\begin{array}{@{}ccc@{}}1 & 1 & 1\\ m_{1} & m_{2} & m_{3}\\ m_{1}\text{e}^{m_{1}} & m_{2}\text{e}^{m_{2}} & m_{3}\text{e}^{m_{3}}\end{array}\right)\left(\begin{array}{@{}c@{}}A_{1}\\ A_{2}\\ A_{3}\end{array}\right)=\left(\begin{array}{@{}c@{}}0\\ 0\\ 0\end{array}\right).
\end{eqnarray}
$$

A necessary condition for a non-trivial solution $$A_j$$ to exist is that the boundary-condition matrix M has zero determinant. (i.e, det M=0)



### Perturbation Porosity

We can get the porosity field using eq(13), 

$$
-SK_0 \sum_{j=1}^3[(m_j^2-k^2) A_j e^{\sigma t +ikx +m_jz}] +n w_0 {\partial \phi_1 \over \partial z} = -M e^{\sigma t +ikx +m_jz}
$$

Since $$M <<1$$, we neglect the $$-MP_1$$ term on the right. Eq(30) becomes,

$$
{\partial \phi_1 \over \partial z} = {S K_0 \over n w_0 } (m_j^2-k^2) \sum_{j=1}^3 A_j e^{\sigma t +ikx +m_jz}
$$

Since $$K_0$$ and $$w_0$$ is 0,

$$
\phi_1 = {S \over n } {m_j^2-k^2 \over m_j} \sum_{j=1}^3 A_j e^{\sigma t +ikx +m_jz} + f(t, x)
$$

The non z-dependent term could be determined by boundary condition $$\phi_1 =0$$, at $$z =0$$,

$$
f(t, x)= {S \over n } {m_j^2-k^2 \over m_j} \sum_{j=1}^3 A_j e^{\sigma t +ikx }
$$

So, insert eq(33) to eq(32)

$$
\phi_1 = {S \over n } {m_j^2-k^2 \over m_j} \sum_{j=1}^3 A_j e^{\sigma t +ikx } (e^{m_jz} -1)
$$



### Perturbation undersaturation

The perturbation undersaturation field is obtained by eq(7), $${\partial \phi_1 \over \partial t} = P_1 + \chi_1$$,

$$
\chi_1 = {\partial \phi_1 \over \partial t} - P_1 = {S \sigma\over n } {m_j^2-k^2 \over m_j} \sum_{j=1}^3 A_j e^{\sigma t +ikx } (e^{m_jz} -1)-\sum_{j=1}^3 A_j e^{\sigma t +ikx +m_jz}
$$

$$
\chi_1 = {S \sigma\over n } {m_j^2-k^2 \over m_j} \sum_{j=1}^3 A_j e^{\sigma t +ikx } (e^{m_jz} -1)-\sum_{j=1}^3 A_j e^{\sigma t +ikx +m_jz}
$$

The perturbation undersaturation field is,

$$
\chi_1 =\sum_{j=1}^3 [ ({S \sigma\over n } {m_j^2-k^2 \over m_j} e^{m_jz} -1)-e^{m_jz}] A_j e^{\sigma t +ikx }
$$



## Solving with python

According to eq(23), we need to know $$A_j$$, $$\sigma$$, $$k$$,  $$m_j$$ to obtain $$P_1$$, while $$k$$ is given to be 10. 

We can use characteristic polynomial eq(22) to solve $$m_j$$, and use eq (29) to solve $$A_j$$. But  there are two unknown parameters in eq(22), $$\sigma$$ and $$m$$, so the whole process is,

1. First, we define a function $$f(\sigma)$$=0,   
2. Then put the randomly-chosen $$\sigma$$ in eq(22) to solve $$m_1$$,  $$m_2$$,  $$m_3$$.
3. Calculate the determinant of M and return det M.
4. The “fslove” tool will keep finding the right $$\sigma$$ until det M=0, and then get the right $$m_1$$,  $$m_2$$,  $$m_3$$.
5. Put $$m_j$$ to  matrix eq(29) to get $$A_1$$,  $$A_2$$,  $$A_3$$, and have the expression of $$P_1$$ eq(23). 



First, the package is imported

```
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math
```

Define the constant value

```
Da = 100
Pe = 100
S = 1
n = 3
k = 10
K = 1 + pow(k,2) / (Da * Pe)
```

Define the polynomial eq(23) with “np.roots” and return m1, m2, m3

```
def root_m(sigma, Da, K, n, S, k):

    coeff = [sigma / Da,
             sigma * K - n / (Da * S),
             - n * K / S - sigma / Da * pow(k, 2),
             (n - sigma * K) * pow(k, 2)]  # set up the companion matrix for eq 3.7

    m = np.roots(coeff)  # Solve all the roots m
    m1, m2, m3 = np.array_split(m, 3)  # define m1 m2 m3
    m1 = complex(m1)  # define data type of m for M_array
    m2 = complex(m2)
    m3 = complex(m3)

    return m1, m2, m3
```

Define the matrix M with “np.array”.

```
def matrix_M(m1, m2, m3):
    M_array = np.array([[1, 1, 1],
                        [m1, m2, m3],
                        [m1 * np.exp(m1), m2 * np.exp(m2), m3 * np.exp(m3)]], dtype='complex')
    return M_array
```

Define the determinant of matrix M with “np.linalg.det”. It should be note that det M is purely imaginary. We need to return the imaginary part of the result.

```
def determinant(sigma):

    m1,m2,m3 = root_m(sigma, Da, K, n, S, k)
    M_array = matrix_M(m1, m2, m3)
    det = np.linalg.det(M_array)  # Calculate det(M)

    return det.imag
```

Use “opt.root_scalar” to find the sigma that makes determinant to be zero. We set a sigma scale of 2.5~2.9.

```
sigma = opt.root_scalar(determinant, bracket=[2.5, 2.9])  # find sigma
sigma = (sigma.root)
```

With the right $$\sigma$$, $$m_j$$, we can solve $$A_j$$ by eq(29) $$\text{M} \times  A_j =0$$. 

Here we use **Singular value decomposition**, if M is a $$m \times n$$ matrix, we can find decomposition $$M = U \sum V^*$$, 

$$
M = 

\begin{pmatrix}
u1, u3, u3
\end{pmatrix}

\begin{pmatrix}
\sigma_1 & 0 & 0 \\
0 & \sigma_2 & 0 \\
0 & 0& \sigma_3 
\end{pmatrix}

\begin{pmatrix}
v_1^* \\
v_2^* \\ 
v_3^* \\
\end{pmatrix}
$$

Where U is an $$ m\times m $$ complex unitary matrix, $$U U^T =I$$. 

$$ \sum$$ is an $$ m\times n $$ rectangular diagonal matrix with non-negative real numbers on the diagonal, $$\sigma_i$$ is the singluar value of M.

and $$V^*$$ is  a conjugate transpose of an $$n\times n$$ complex unitary matrix, $$V V^T= I$$. 

We can use “u, s, vh = np.linalg.svd” to get the decomposition form.  Then use “np.conj” to return the complex conjugate $$v_3$$.  (e.g., conjugate of complex number $$2+5j$$ is $$2-5j$$ )

```
m1,m2,m3 = root_m(sigma, Da, K, n, S, k)
M_array = matrix_M(m1, m2, m3)

def prefactorA(M_array):

    u, s, vh = np.linalg.svd(M_array) 
    A = np.conj(vh[2, :])   # M * v_3 =0

    return A

Aj = prefactorA(M_array)
mj = root_m(sigma, Da, K, n, S, k)
Bj = S * (np.power(mj, 2) - pow(k, 2)) / mj / n
```

Create the mesh and calculate the Pressure, porosity and undersaturation.

```
xv = np.linspace(0, 2.0* np.pi/k, 64)
yv = np.linspace(0, 1, 64)
x, y = np.meshgrid(xv, yv)

P = np.exp(1.0j * k * x) * (Aj[0] * np.exp(mj[0] * y)+ Aj[1] * np.exp(mj[1] * y) + Aj[2] * np.exp(mj[2] * y))
P = P.real

h = np.exp(1.0j * k * x) * (Aj[0] * Bj[0] * (np.exp(mj[0] * y)-1) + Aj[1] * Bj[1] * (np.exp(mj[1] * y)-1) + Aj[2] * Bj[2] * (np.exp(mj[2] * y)-1))
h = h.real

Ch = np.exp(1.0j * k * x) * (Aj[0] * (Bj[0] * sigma * (np.exp(mj[0] * y)-1) - np.exp(mj[0] * y)) + Aj[1] * (Bj[1] * sigma * (np.exp(mj[1] * y)-1) - np.exp(mj[1] * y)) + Aj[2] * (Bj[2] * sigma * (np.exp(mj[2] * y)-1) - np.exp(mj[2] * y)))
Ch =Ch.real
```




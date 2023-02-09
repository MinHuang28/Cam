---
---
title: 'FEniCS-Solve RII problem'
date: 2023-02-09
permalink: 
tags:
  - FEniCS
  - RII
---

---

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

## Step porosity as initial condition

A region of partially molten rock of height H in the z direction with a flux of liquid phase from beneath. 

The porosity is given by the piecewise constant function,

$$
\phi (z)=
\begin{cases} 
c_1 \phi_0 & \text {if } z> 1/2 \\
c_2 \phi_0 &\text{if } z \leq 1/2
\end{cases}
$$

Dimensionless mobility  $$K= \phi(z)^n$$.

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
\int_\Omega MP P_test+\int_\Omega P_test {\nabla \cdot [K(\widehat z-S\nabla P)]}=0
$$

Apply the boundary condition $$\frac{\partial P}{\partial n} =0$$, on $$\Gamma$$, we get the weak form,

$$
\int_\Omega MP P_test+\int_\Gamma P_test K(\widehat z-S\nabla P) \cdot n - \int_\Omega \nabla P_test \cdot K(\widehat z -S\nabla P) =0
$$

$$
\int_\Omega MP P_test + \int_\Omega KS \nabla P_test \cdot \nabla P = \int_\Omega K {\partial P_test \over \partial z} -\int_\Gamma P_test K(\widehat z\cdot n)
$$

## Step_2: Solve reaction-infiltration(\chi)

With Darcy flux $$q = K(\widehat z-S \nabla P)$$, the equations becomes,

$$
q \cdot ({\nabla \chi \over Da} - \hat z)= {1\over DaPe} {\partial \over \partial x}(\phi{\partial \chi \over \partial x})-\chi
$$

Multiply the whole equation by a **test function** C_test, and integrate the whole equation over the domain $$\Omega$$. 

$$
\int_\Omega q \cdot ({\nabla \chi \over Da} - \hat z) C_test - \int_\Omega[{1\over DaPe} {\partial \over \partial x}(\phi{\partial \chi \over \partial x})-\chi] C_test =0
$$

$$
\int_\Omega q\cdot {\nabla \chi \over Da} C_test - \int_\Omega q \cdot \hat z C_test - \int_\Omega{\phi \over DaPe} {\partial \over \partial x} ({\partial \chi \over \partial x}) C_test+\int_\Omega \chi C_test =0
$$

Integrating by parts,

$$
\int_\Omega {\partial \over \partial x}({\partial \chi \over \partial x})C_test = [{\partial \chi \over \partial x} C_test]_{\partial \Omega} -\int_\Omega {\partial \chi \over \partial x} {\partial C_test \over \partial x}
$$

So the weak form is, 

$$
\int_\Omega q \cdot {\nabla \chi \over Da} C_test - \int_\Omega q \cdot \hat z C_test + \int_\Omega{\phi \over DaPe} {\partial C_test \over \partial x} {\partial \chi \over \partial x}+\int_\Omega \chi C_test =0
$$

$$
\int_\Omega q \cdot {\nabla \chi \over Da} C_test + \int_\Omega{\phi \over DaPe} {\partial C_test \over \partial x} {\partial \chi \over \partial x}+\int_\Omega \chi C_test = \int_\Omega q \cdot \hat z C_test
$$

## Step3: Upgrade porosity with time

Governing equation is 

$$
{\partial \phi \over \partial t} = P + \chi
$$

Discretize the time derivative by a finite difference approximation

$$
{\partial\over \partial t} \phi ^n = P ^n
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

$$\frac{\partial P}{\partial n} =0$$, on $$\Gamma$$

$$\Omega =[0,1] \times [0,1] $$

---

## Code

```
from fenics import DirichletBC, Function, UserExpression, RectangleMesh, FunctionSpace, Point,\
    Constant, project, DOLFIN_EPS, File, LinearVariationalProblem, LinearVariationalSolver, FacetNormal, \
    TrialFunction, TestFunction, dx, grad, dot, Dx, inner, ds, plot
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

# define initial porosity by expression
class Phi_Initial(UserExpression):
    def __init__(self, amplitude, A, m, w, **kwargs):
        super().__init__(**kwargs)
        self.amplitude = amplitude
        self.A = A
        self.m = m
        self.w = w

    def eval(self, values, xx):
        x = xx[0]
        y = xx[1]
        P = np.exp(1.0j * self.w * x) * sum(
            [AA * np.exp(mm * y) for AA, mm in zip(self.A, self.m)]
        )
        values[0] = 1.0 - self.amplitude * P.real

    def value_shape(self):
        return ()

# define parameter
Da = 100                    # Damkohler number Da=100
Pe = 100                    # Peclet number Pe=100
b = 1e-6                    # Solubility gradient \beta
H = 1e4                     # Melting region depth H
M = b * H                   # Solubility M=0.01
d = 1e5                     # Compaction length \delta=10^5
S = M * pow(d,2) / pow(H,2) # Stiffness S=1
n = 3                       # n =2/3
w = 10                      # horizontal wavenumber w=10
K = 1 + pow(w,2) / (Da * Pe)

# define linear analysis
def root_m(sigma, Da, K, n, S, w):

    coeff = [sigma / Da,
             sigma * K - n / (Da * S),
             - n * K / S - sigma / Da * pow(w, 2),
             (n - sigma * K) * pow(w, 2)]  # set up the companion matrix for eq 3.7

    m = np.roots(coeff)  # Solve all the roots m
    m1, m2, m3 = np.array_split(m, 3)  # define m1 m2 m3
    m1 = complex(m1)  # define data type of m for M_array
    m2 = complex(m2)
    m3 = complex(m3)

    return m1, m2, m3

def matrix_M(m1, m2, m3):
    M_array = np.array([[1, 1, 1],
                        [m1, m2, m3],
                        [m1 * np.exp(m1), m2 * np.exp(m2), m3 * np.exp(m3)]], dtype='complex')
    return M_array

def determinant(sigma):

    m1,m2,m3 = root_m(sigma, Da, K, n, S, w)
    M_array = matrix_M(m1, m2, m3)
    det = np.linalg.det(M_array)  # Calculate det(M)

    return det.imag

sigma = opt.root_scalar(determinant, bracket=[2.5, 2.9])  # find sigma
sigma = sigma.root

m1, m2, m3 = root_m(sigma, Da, K, n, S, w)
M_array = matrix_M(m1, m2, m3)

def prefactor(M_array):

    u, s, vh = np.linalg.svd(M_array)
    A = np.conj(vh[2, :])
    A1 = A[0]
    A2 = A[1]
    A3 = A[2]

    return A1, A2, A3

A = prefactor(M_array)
m = root_m(sigma, Da, K, n, S, w)

# define mesh
mesh = RectangleMesh(Point(0, 0), Point(2.0* np.pi/10, 1), 64, 64)
V = FunctionSpace(mesh,"CG", 1)
nn = FacetNormal(mesh)       # n vector is the facet normal
z = Constant((0,1))     # unit vector z

# define initial value of \phi
phi_i = Phi_Initial(0.01, A, m, w, degree=2)
h_n = project(phi_i, V)
#plot(h_n, mesh=mesh)
#plt.show()

# define boundary condition
def on_bottom(x, on_boundary):
    return x[1] <= DOLFIN_EPS and on_boundary

C_bottom = 1.0   # \chi =1
bc_C = DirichletBC(V, C_bottom, on_bottom)

h_bottom = 1.0   # \phi =1
bc_h = DirichletBC(V, h_bottom, on_bottom)

# define simulation time
t = 0.0
T = 0.1         # total simulation time
t_num = 200      # number of time steps
dt = T /t_num         # time stepping

# define buoyancy-driven compaction problem
def compaction(V, h_n, M, S):
    # Test and Trial function
    P_n = TrialFunction(V)
    P_test = TestFunction(V)
    k_n = pow(h_n, 3)

    # define linear and bilinear form
    aP = (M * P_n * P_test + k_n * S * inner(grad(P_n), grad(P_test))) * dx
    LP = k_n * Dx(P_test, 1) * dx - k_n * P_test * dot(z, nn) * ds

    return (aP, LP)

P_n = Function(V)
aP, LP = compaction(V, h_n, M, S)

def reaction(V, Da, Pe, P_n, h_n, S):

    C_n = TrialFunction(V)  # Trial function undersaturation \chi
    C_test = TestFunction(V)  # test function
    q_n = pow(h_n,3) * (z - S * grad(P_n))

    # Define bilinear and linear forms
    aC = ((1 / Da) * dot(q_n, grad(C_n)) * C_test + (h_n / (Da * Pe)) * Dx(C_test, 0) * Dx(C_n, 0) + C_n * C_test) * dx
    LC = dot(q_n, z) * C_test * dx

    return aC, LC

# define Time-dependent infiltration problem (\partial \phi / \partial t = \chi)
def time(V, h_n, P_n, C_n, dt):

    h = TrialFunction(V)  # porosity \phi
    h_test = TestFunction(V)  # test function

    ah = h * h_test * dx
    Lh = (h_n + (P_n + C_n) * dt) * h_test * dx

    return ah, Lh

C_n = Function(V)
h = Function(V)

# define variation form
aC, LC = reaction(V, Da, Pe, P_n, h_n, S)
ah, Lh = time(V, h_n, P_n, C_n, dt)

# Output file
plt.figure()

Pfile = File("result/Pressure.pvd")
Cfile = File("result/undersaturation.pvd")
hfile = File("result/porosity.pvd")

h = h_n

# main loop, Time-stepping
for n in range(t_num):

    # Solve
    problem1 = LinearVariationalProblem(aP, LP, P_n, bc_h)
    solver1 = LinearVariationalSolver(problem1)
    solver1.solve()

    problem2 = LinearVariationalProblem(aC, LC, C_n, bc_C)
    solver2 = LinearVariationalSolver(problem2)
    solver2.solve()

    problem3 = LinearVariationalProblem(ah, Lh, h)
    solver3 = LinearVariationalSolver(problem3)
    solver3.solve()

    # save file
    Pfile << P_n
    Cfile << C_n
    hfile << h

    # update current time
    t += dt

    h_n.assign(h)
```
<script type="text/javascript" src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

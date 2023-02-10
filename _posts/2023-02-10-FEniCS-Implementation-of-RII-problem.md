## FEniCS Implementation for RII problem

FEniCS is an open-source software platform for solving partial differential equations (PDEs) using the finite element method (FEM). The primary programming language for using FEniCS is python, which is quite simple to use. 

Therefore, based on the mathematical equation that describes the RII problem in the previous blog, we could code up the PDEs with FEniCS and python.

First, the FEniCS is imported.

Note: Some older tutorials might import everything the from Fenics package (e.g., from fenics import *), which is not a good practice. It is preferable to list any commands you might use and to be fully aware of their usage in your code.

```
from fenics import DirichletBC, Function, UserExpression, RectangleMesh, FunctionSpace, Point, Constant, project, DOLFIN_EPS, File, FacetNormal,TrialFunction, TestFunction, dx, grad, dot, Dx, inner, ds, plot, solve
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
```

The plotting tool and the scipy and numpy packages are also imported.

Then, create a mesh covering the rectangle square. We will let the mesh consist of 64 $$\times$$ 64 squares and define function space.

```
mesh = RectangleMesh(Point(0, 0), Point(2.0* np.pi/10, 1), 64, 64)
V = FunctionSpace(mesh,"CG", 1)
nn = FacetNormal(mesh)       # n vector is the facet normal
z = Constant((0,1))     # unit vector z
```

Initial porosity is defined by linear stability analysis. (Shown in previous blog)

```
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
```

Define values of parameters that we need to use

```
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

# define simulation time
t = 0.0
T = 0.1         # total simulation time
t_num = 200      # number of time steps
dt = T /t_num         # time stepping
```

Solve A, m for initial porosity condition

```
# Define root finding
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
    
# Define matrix M
def matrix_M(m1, m2, m3):
    M_array = np.array([[1, 1, 1],
                        [m1, m2, m3],
                        [m1 * np.exp(m1), m2 * np.exp(m2), m3 * np.exp(m3)]], dtype='complex')
    return M_array
   
# Solve the detM=0
def determinant(sigma):

    m1,m2,m3 = root_m(sigma, Da, K, n, S, w)
    M_array = matrix_M(m1, m2, m3)
    det = np.linalg.det(M_array)  # Calculate det(M)

    return det.imag

sigma = opt.root_scalar(determinant, bracket=[2.5, 2.9])  # find sigma
sigma = sigma.root

m1, m2, m3 = root_m(sigma, Da, K, n, S, w)
M_array = matrix_M(m1, m2, m3)

# With sigma and mj, solve Aj
def prefactor(M_array):

    u, s, vh = np.linalg.svd(M_array)
    A = np.conj(vh[2, :])
    A1 = A[0]
    A2 = A[1]
    A3 = A[2]

    return A1, A2, A3

A = prefactor(M_array)
m = root_m(sigma, Da, K, n, S, w)

# define initial value of \phi
phi_i = Phi_Initial(0.01, A, m, w, degree=2)
h_n = project(phi_i, V)
#plot(h_n, mesh=mesh)
#plt.show()
```

 Define boundary conditions with
$$
\phi =1, \chi =1, {\partial P\over \partial z} =0,  (z =0)\\
{\partial P\over \partial z} =0,  (z =1)
$$

```
# define boundary condition
def on_bottom(x, on_boundary):
    return x[1] <= DOLFIN_EPS and on_boundary

C_bottom = 1.0   # \chi =1
bc_C = DirichletBC(V, C_bottom, on_bottom)

h_bottom = 1.0   # \phi =1
bc_h = DirichletBC(V, h_bottom, on_bottom)
```

Define buoyancy-driven compacting problem
$$
\int_\Omega MP P_1 + \int_\Omega KS \nabla P_1 \cdot \nabla P = \int_\Omega K {\partial P_1 \over \partial z} -\int_\Gamma P_1 K(\widehat z\cdot n)
$$
In FEniCS, the variational formulation is defined by trial and test functions. The trial function is a function that approximates the solution to the PDE, and the test function is a function that is used to integrate the weak form of the PDE over the domain.

```
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
```

Define reaction-infiltration problem
$$
\int_\Omega q \cdot {\nabla \chi \over Da} C_1+ \int_\Omega{\phi \over DaPe} {\partial C_1 \over \partial x} {\partial \chi \over \partial x}+\int_\Omega \chi C_1 = \int_\Omega q \cdot \hat z C_1
$$

```
def reaction(V, Da, Pe, P_n, h_n, S):

    C_n = TrialFunction(V)  # Trial function undersaturation \chi
    C_test = TestFunction(V)  # test function
    q_n = pow(h_n,3) * (z - S * grad(P_n))

    # Define bilinear and linear forms
    aC = ((1 / Da) * dot(q_n, grad(C_n)) * C_test + (h_n / (Da * Pe)) * Dx(C_test, 0) * Dx(C_n, 0) + C_n * C_test) * dx
    LC = dot(q_n, z) * C_test * dx

    return aC, LC
```

Define Time-dependent infiltration problem
$$
a(\phi,R)=\int_\Omega \phi R dx \\
L_n(R) = \int_\Omega [\phi^n +(P+\chi)^n\Delta t] R dx
$$

```
def time(V, h_n, P_n, C_n, dt):

    h = TrialFunction(V)  # porosity \phi
    h_test = TestFunction(V)  # test function

    ah = h * h_test * dx
    Lh = (h_n + (P_n + C_n) * dt) * h_test * dx

    return ah, Lh
```

```
C_n = Function(V)
h = Function(V)

# define variation form
aC, LC = reaction(V, Da, Pe, P_n, h_n, S)
ah, Lh = time(V, h_n, P_n, C_n, dt)
```

Define output file

```
# Output file
plt.figure()

Pfile = File("result/Pressure.pvd")
Cfile = File("result/undersaturation.pvd")
hfile = File("result/porosity.pvd")
```

Define time-stepping main loop

```
h = h_n

# main loop, Time-stepping
for n in range(t_num):

    # Solve
	solve(aP == LP, P_n, bc_h)
	solve(aC == LC, C_n, bc_C)
	solve(ah == Lh, h)

    # save file
    Pfile << P_n
    Cfile << C_n
    hfile << h

    # update current time
    t += dt# 

    h_n.assign(h)
```



<script type="text/javascript" src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

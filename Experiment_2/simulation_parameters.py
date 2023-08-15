import sys

import yaml

# setting path
sys.path.append('../')

import numpy as np
from dolfin import (File, Function, FunctionSpace, Mesh, MeshValueCollection,
                    XDMFFile, cpp)
from PML_base.PoroelasticParameters import Set3
from PML_base.Source import neumann_ricker_pulse
from PML_base.TimeStepping import TimeStepping3

# Import geometric parameters
with open('PARAMETERS.yaml') as file:
    params = yaml.load(file, Loader=yaml.FullLoader)
Lx        = params[0]['Lx']
Ly        = params[0]['Ly']  
LPML      = params[0]['LPML']

# Central frequency
fr      = params[3]['fr']
omega_p = 2*np.pi*fr
omega_c = 1.059095*omega_p    # Central frequency
omega_b = 0.577472*omega_p    # half-bandwidth

# Import mesh
mesh = Mesh()
with XDMFFile("mesh/mesh_hybrid.xdmf") as infile:
    infile.read(mesh)

# Markers for layers
mvc = MeshValueCollection("size_t", mesh, mesh.geometry().dim())
with XDMFFile("mesh/mf_hybrid_1.xdmf") as infile:
    infile.read(mvc, "markers")
mf1 = cpp.mesh.MeshFunctionSizet(mesh, mvc)

# Get cell indices for each layer
layer_1 = np.where(mf1.array()[:] == 3)[0]
layer_2 = np.where(mf1.array()[:] == 4)[0]
layer_3 = np.where(mf1.array()[:] == 5)[0]
layer_4 = np.where(mf1.array()[:] == 6)[0]
layers = [layer_1, layer_2, layer_3, layer_4]

# Function space for physical parameters
P = FunctionSpace(mesh, 'DG', 0)

# Porous medium parameters
params = Set3(V=P, layers=layers)
rho_s = params.rho_s  # solid densitiy
rho_f = params.rho_s  # fluid density
K_s   = params.K_s    # bulk moduli of fluid
K_f   = params.K_f    # bulk moduli of the solid
K_b   = params.K_b    # bulk moduli of the dry porous frame
mu_b  = params.mu_b   # shear moduli of the dry porous frame
T     = params.T      # tortuosity
kappa = params.kappa  # permeability of the fluid
phi   = params.phi    # porosity
alpha = params.alpha  # Biot-Willis coefficient
rho   = params.rho    # effective density
rho_w = params.rho_w 
eta   = params.eta    # dynamic viscosity of the fluid
M     = params.M      # fluid-solid coupling modulus

# Wave speeds
c1p = params.c1p.vector().max()
c2p = params.c2p.vector().min()
cs  = params.cs.vector().min()
# c1p = params.V1.vector().min()
# cs = params.V3.vector().min()
File('output/debug/c1p.pvd') << params.V1
File('output/debug/c2p.pvd') << params.V2
File('output/debug/cs.pvd') << params.V3
File('output/debug/fc.pvd') << params.fc

# Wave lengths
omega_ = omega_c + omega_b              # Frequency to evaluate
L1p = 2*np.pi*c1p/omega_p               # P-wavelength
L2p = 2*np.pi*c2p/omega_p               # P-wavelength
Ls = 2*np.pi*min([cs,c2p,c1p])/omega_p  # S-wavelength

# Stability evaluation
order  = 1                            # Order of the Lagrange polynomials
n      = 12                           # Number of elements wavelength
cell_d = order*(Ls/n)                 # Maxmimum cell diameter
CFL    = 0.75                         # CFL number
dt     = CFL/(c1p/cell_d)             # CFL estability condition

# Mesh element size for visualization
mesh_size = Function(P)
mesh_size.vector()[:] = 2*np.pi*order*(params.cs.vector()/omega_p)/n
File('output/debug/mesh_size.pvd') << mesh_size

# Ricker pulse
g = neumann_ricker_pulse(A=1e+04, t=0.0, omega=omega_p, degree=2)

# Dimensions for the unbounded experiment
#   - d: diagonal of the rectangle regular domain
#   - T: time that takes the slowest wave to reach the fartest boundary
#       of the regular domain
#   - R: radius of the extended domain
d = np.sqrt((Lx/2)**2 + Ly**2)
# T = d/min([cs,c1p,c2p]) + g.cutoff
T = 2.0
R = 0.5*(c1p*T + 2*d) #0.5*(c1p*T + 3*d)

# print('\n[*] Wave velocities')
# print('   Fast primary wave speed (c1p):  {:.2f} [m/s]'.format(params.c1p))
# print('   Slow primary wave speed (c2p):  {:.2f} [m/s]'.format(params.c2p))
# print('   Secondary wave speed (cs):      {:.2f} [m/s]'.format(params.cs))

print('\n[*] Extended domain')
print('   Radius of the extended domain: {:.2f} [m]'.format(R))
print('   Simulation time:               {:.2f} [s]'.format(T))

print('\n[*] Stability (CFL) conditions')
print('   Maximun element size:   {:.6f} [m]'.format(cell_d))
print('   Maximun time step:      {:.6f} [s]'.format(dt))

print('L1p = {:.1f}'.format(L1p))
print('L2p = {:.1f}'.format(L2p))
print('Ls = {:.1f}'.format(Ls))

P = FunctionSpace(mesh, 'DG', 0)
lmbda_b = Function(P)
lmbda_b.vector()[:] = K_b.vector() - 2/3*mu_b.vector() + alpha.vector()*alpha.vector()*M.vector()

nu = Function(P)
nu.vector()[:] = lmbda_b.vector()
nu.vector()[:] /= (2*(lmbda_b.vector()+mu_b.vector()))
File('output/debug/nu_b.pvd') << nu
File('output/debug/lmbda_b.pvd') << lmbda_b

import sys

# setting path
sys.path.append('../')

import numpy as np
from PML_base.PoroelasticParameters import Set1

# Poisson coefficient
nu = 0.35

# Shear wave velocities
# cs = np.array([300, 400, 500, 700]) # m/s
cs = np.array([1000, 1330, 1660, 2330]) # m/s

# Porous medium parameters
params = Set1()
rho_s = params.rho_s.values()[0]  # solid densitiy
rho_f = params.rho_s.values()[0]  # fluid density
K_s   = params.K_s.values()[0]    # bulk moduli of fluid
K_f   = params.K_f.values()[0]    # bulk moduli of the solid
K_b   = params.K_b.values()[0]    # bulk moduli of the dry porous frame
mu_b  = params.mu_b.values()[0]   # shear moduli of the dry porous frame
T     = params.T.values()[0]      # tortuosity
kappa = params.kappa.values()[0]  # permeability of the fluid
phi   = params.phi.values()[0]    # porosity
alpha = params.alpha.values()[0]  # Biot-Willis coefficient
rho   = params.rho.values()[0]    # effective density
rho_w = params.rho_w.values()[0] 
eta   = params.eta.values()[0]    # dynamic viscosity of the fluid
M     = params.M.values()[0]      # fluid-solid coupling modulus

# Approximation of Lame parameters
mu_b = params.Z1*np.power(cs,2)/rho_w
K_b  = 2*mu_b*nu/(1.0-2.0*nu)
print(mu_b/mu_b[0])
print(K_b/K_b[0])

# Print wave velocities
print('\n[*] Wave velocities')
print('   Fast primary wave speed (c1p):  {:.2f} [m/s]'.format(params.c1p))
print('   Slow primary wave speed (c2p):  {:.2f} [m/s]'.format(params.c2p))
print('   Secondary wave speed (cs):      {:.2f} [m/s]'.format(params.cs))

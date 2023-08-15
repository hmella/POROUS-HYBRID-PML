# Copyright (C) 2019 Hernan Mella

import numpy as np
from dolfin import (DOLFIN_EPS, Constant, SpatialCoordinate, UserExpression,
                    conditional, interpolate)


#############################
# Uniform PML functions
#############################
# Get PML functions
def PMLFunctions(domain,PhysParams,form=1,R=1e-8,m=2,alpha_max=None,beta_max=None,PMLType='PML',Pxy=0.15,Pyx=0.15):

  # Retrieve geometry parameters
  Lpml = domain.LPML

  # Get characteristic length
  b = domain.b

  # Get characteristic velocity
  cp = PhysParams.V1.vector().max()
  print('[PML] Velocity (cp) considered for PML functions: {:.0f}'.format(cp))

  # PML functions parameters
  # Form 1 represent the usual polynomial form of the PML functions
  # (see for instance Kucuckovan et al - 2011), whereas form 2 denotes
  # PML functions that have been previously tested on poroelastodynamic
  # problems (see He et al - 2019)
  if form==1:
    # S. Kucukcoban and L. F. Kallivokas, “Mixed perfectly matched layers for direct transient analysis in 2D elastic heterogeneous media,” Comput. Meth. Appl. Mech. Eng., vol. 200, no. 1, pp. 57–76, 2011.
    alpha = alpha_max if alpha_max != None else (m + 1)*b/(2*Lpml)*np.log(1/R)
    beta  = beta_max if beta_max != None else (m + 1)*cp/(2*Lpml)*np.log(1/R)
  elif form==2:
    # Y. He, T. Chen, and J. Gao, “Unsplit perfectly matched layer absorbing boundary conditions for second-order poroelastic wave equations,” Wave Motion, vol. 89, pp. 116–130, 2019, doi: 10.1016/j.wavemoti.2019.01.004.
    alpha = alpha_max-1 if alpha_max != None else alpha_max - 1.0
    beta  = beta_max if beta_max != None else (2.0*cp/Lpml)*np.log(1/R)
  # Change beta coefficient for Optimal PML
  if PMLType is 'OPML':
    beta  = cp
  print('[PML] alpha considered for PML functions: {:.2f}'.format(alpha),flush=True)
  print('[PML] beta considered for PML functions: {:.2f}'.format(beta),flush=True)


  # Mesh coordinates
  x = SpatialCoordinate(PhysParams.V1.function_space().mesh())

  # Attenuation and stretching functions and tensors
  if PhysParams.V1.function_space().mesh().geometry().dim()==2:
    ax, ay, bx, by = PML(domain,x,alpha,beta,m,form,PMLType,Pxy,Pyx)
    return ax, ay, bx, by
  elif PhysParams.V1.function_space().mesh().geometry().dim()==3:
    ax, ay, az, bx, by, bz = PML3(domain,x,alpha,beta,m,form,PMLType,Pxy,Pyx)
    return ax, ay, az, bx, by, bz


# Stretching functions
def PML(domain,x,alpha,beta,m,form,PMLType,Pxy,Pyx):
  # Domain dimensions
  Lx = domain.Lx
  Ly = domain.Ly
  Lp = domain.LPML

  # Domain coordinates
  y = x[1]
  x = x[0]

  # Impose quadratic an cubic profiles for alpha and 
  # beta functions respectively
  if form==2:
    (ma, mb) = (2, 3)
  else:
    (ma, mb) = (m, m)

  # PML functions
  # \gamma(s) = 1 + \alpha_s(s) + \frac{\beta_s(s)}{i\omega}
  ax = Constant(1.0) \
      + conditional(x < -0.5*Lx, alpha*pow((abs(x) - 0.5*Lx)/Lp, ma), 0) \
      + conditional(x > 0.5*Lx, alpha*pow((abs(x) - 0.5*Lx)/Lp, ma), 0)
  ay = Constant(1.0) \
      + conditional(y < -Ly, alpha*pow((abs(y) - Ly)/Lp, ma), 0)
  bx = conditional(x < -0.5*Lx, beta*pow((abs(x) - 0.5*Lx)/Lp, mb), 0) \
      + conditional(x > 0.5*Lx, beta*pow((abs(x) - 0.5*Lx)/Lp, mb), 0)
  by = conditional(y < -Ly, beta*pow((abs(y) - Ly)/Lp, mb), 0)

  # Add M-PML functions
  if PMLType is 'MPML':
    ax += Constant(Pyx)*(ay - Constant(1.0))
    ay += Constant(Pxy)*(ax - Constant(1.0))
    bx += Constant(Pyx)*by
    by += Constant(Pxy)*bx
    print('[PML] Factors considered for multiaxial PML functions (Pxy, Pyx): ({:.4f},{:.4f})'.format(Pxy,Pyx),flush=True)

  return ax, ay, bx, by
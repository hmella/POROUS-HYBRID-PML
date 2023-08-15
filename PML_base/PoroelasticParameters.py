import numpy as np
from dolfin import Function


# Base class for physical parameters
class PhysicalParameters:
  def __init__(self, rho_s=[], rho_f=[], K_s=[], K_f=[], K_b=[], mu_b=[],
               T=[], kappa=[], phi=[], eta=[], V=[], layers=[]):
    self.V = V            # Function space
    self.layers = layers
    self.rho_s_arr = rho_s # solid densitiy [kg/m^3]
    self.rho_f_arr = rho_f # fluid density [kg/m^3]
    self.K_s_arr = K_s     # bulk moduli of the solid [N/m^2]
    self.K_f_arr = K_f     # bulk moduli of fluid [N/m^2]
    self.K_b_arr = K_b     # bulk moduli of the dry porous frame [N/m^2]
    self.mu_b_arr = mu_b   # shear moduli of the dry porous frame [N/m^2]
    self.T_arr = T         # tortuosity [a.u.]
    self.kappa_arr = kappa # permeability of the fluid [m^2]
    self.phi_arr = phi     # porosity [a.u.]
    self.eta_arr = eta     # dynamic viscosity of the fluid [N s]

  # Pass porous medium parameters across meshes
  def build_parameters(self):

    # Base parameters
    self.rho_s = LayeredParameters(self.rho_s_arr, self.layers, self.V)
    self.rho_f = LayeredParameters(self.rho_f_arr, self.layers, self.V)
    self.K_s   = LayeredParameters(self.K_s_arr, self.layers, self.V)
    self.K_f   = LayeredParameters(self.K_f_arr, self.layers, self.V)
    self.K_b   = LayeredParameters(self.K_b_arr, self.layers, self.V)
    self.mu_b  = LayeredParameters(self.mu_b_arr, self.layers, self.V) 
    self.T     = LayeredParameters(self.T_arr, self.layers, self.V) 
    self.kappa = LayeredParameters(self.kappa_arr, self.layers, self.V) 
    self.phi  = LayeredParameters(self.phi_arr, self.layers, self.V)
    self.eta = LayeredParameters(self.eta_arr, self.layers, self.V)

    # Calculate derived parameters
    self.alpha = Function(self.V)   # Biot-Willis coefficient [a.u.]
    self.rho = Function(self.V)     # effective density [kg/m^3]
    self.rho_w = Function(self.V)   # [kg/m^3]     
    self.M = Function(self.V)       # fluid-solid coupling modulus [a.u.]
    self.fc = Function(self.V)      # Characteristic frequency [Hz]

    K_b = self.K_b.vector().get_local()
    K_s = self.K_s.vector().get_local()
    K_f = self.K_f.vector().get_local()
    rho_s = self.rho_s.vector().get_local()
    rho_f = self.rho_f.vector().get_local()
    phi = self.phi.vector().get_local()
    T   = self.T.vector().get_local()
    alpha = 1.0 - K_b/K_s

    self.alpha.vector()[:] = alpha
    self.rho.vector()[:]   = rho_s*(1.0 - phi) + rho_f*phi
    self.rho_w.vector()[:] = T*rho_f/phi
    self.M.vector()[:]     = 1.0/(phi/K_f + (alpha - phi)/K_s)

    # Needed parameters for wave velocities
    rho   = self.rho.vector().get_local()
    rho_w = self.rho_w.vector().get_local()
    mu_b  = self.mu_b.vector().get_local()
    eta   = self.eta.vector().get_local()
    M     = self.M.vector().get_local()
    lmbda = K_b + alpha**2*M - 2/3*mu_b

    # Wave velocities
    self.Z1 = Function(self.V)
    self.Z1.vector()[:] = rho_w*rho - rho_f**2
    Z2 = -2*rho_f*alpha*M + rho*M + rho_w*lmbda + 2*rho_w*mu_b
    Z3 = rho*(4*alpha**2*rho_w - 4*alpha*rho_f + rho)*M**2 - 2*(2*alpha*rho_w*rho_f + rho_w*rho - 2*rho_f**2)*M*(2*mu_b + lmbda) + rho_w**2*(2*mu_b + lmbda)**2
    self.c1p = Function(self.V)
    self.c2p = Function(self.V)
    self.cs = Function(self.V)
    self.c1p.vector()[:] = np.sqrt((Z2 + np.sqrt(Z3))/(2*self.Z1.vector().get_local()))
    self.c2p.vector()[:] = np.sqrt((Z2 - np.sqrt(Z3))/(2*self.Z1.vector().get_local()))
    self.cs.vector()[:]  = np.sqrt(rho_w*mu_b/self.Z1.vector().get_local())

    # Velocities for paraxial BC
    self.V1 = self.c1p
    self.V2 = self.c2p
    self.V3 = self.cs

    # Characteristic frequency
    kappa = self.kappa.vector().get_local()
    self.fc.vector()[:] = eta*phi/(2*np.pi*T*rho_f*kappa)


class Set1(PhysicalParameters):
  def __init__(self, **kwargs):
    super(Set1, self).__init__(**kwargs)

    # Physical parameters taken from the literature
    param = 4
    if param == 1:
      # Bottom layer in Table 1 of:
      # [1] Y. He, T. Chen, and J. Gao, “Unsplit perfectly matched layer absorbing boundary conditions for second-order poroelastic wave equations,” Wave Motion, vol. 89, pp. 116–130, 2019, doi: 10.1016/j.wavemoti.2019.01.004.
      self.rho_s_arr = np.array([2700.0])
      self.rho_f_arr = np.array([1040.0])
      self.K_s_arr   = np.array([2.6e+10])
      self.K_f_arr   = np.array([2.5e+9])
      self.K_b_arr   = np.array([2.5e+10])
      self.mu_b_arr  = np.array([1.7e+10])
      self.T_arr     = np.array([2.0])
      self.kappa_arr = np.array([1e-14])
      self.phi_arr   = np.array([0.05])
      self.eta_arr   = np.array([1e-3])
    elif param == 2:
      # Top layer in Table 1 of:
      # [1] Y. He, T. Chen, and J. Gao, “Unsplit perfectly matched layer absorbing boundary conditions for second-order poroelastic wave equations,” Wave Motion, vol. 89, pp. 116–130, 2019, doi: 10.1016/j.wavemoti.2019.01.004.
      self.rho_s_arr = np.array([2670.0])
      self.rho_f_arr = np.array([1040.0])
      self.K_s_arr   = np.array([9.4e+9])
      self.K_f_arr   = np.array([2.5e+9])
      self.K_b_arr   = np.array([8.9e+9])
      self.mu_b_arr  = np.array([5.5e+9])
      self.T_arr     = np.array([2.0])
      self.kappa_arr = np.array([1e-13])
      self.phi_arr   = np.array([0.02])
      self.eta_arr   = np.array([1e-3])
    elif param==3:
      # Table 2 (high Poisson ratio)
      # [1] Y. He, T. Chen, and J. Gao, “Perfectly matched absorbing layer for modelling transient wave propagation in heterogeneous poroelastic media,” Journal of Geophysics and Engineering, Oct. 2019, doi: 10.1093/jge/gxz080.
      self.rho_s_arr = np.array([2650.0])
      self.rho_f_arr = np.array([1040.0])
      self.K_s_arr   = np.array([35.0e+9])
      self.K_f_arr   = np.array([2.4e+9])
      self.K_b_arr   = np.array([4.17e+9])
      self.mu_b_arr  = np.array([1.855e+9])
      self.T_arr     = np.array([2.0])
      self.kappa_arr = np.array([1e-12])
      self.phi_arr   = np.array([0.3])
      self.eta_arr   = np.array([1e-3])
    elif param==4:
      # Table 1
      # [1] N. F. Dudley Ward, T. Lähivaara, and S. Eveson, “A discontinuous Galerkin method for poroelastic wave propagation: The two-dimensional case,” Journal of Computational Physics, vol. 350, pp. 690–727, 2017, doi: 10.1016/j.jcp.2017.08.070.
      self.rho_s_arr = np.array([2650.0])
      self.rho_f_arr = np.array([900.0])
      self.K_s_arr   = np.array([12.0e+9])
      self.K_f_arr   = np.array([2.0e+9])
      self.K_b_arr   = np.array([10.0e+9])
      self.mu_b_arr  = np.array([5.0e+9])
      self.T_arr     = np.array([1.2])
      self.kappa_arr = np.array([1e-12])
      self.phi_arr   = np.array([0.3])
      self.eta_arr   = np.array([1e-3])  
    elif param==5:
      # Table 6
      # [1] N. F. Dudley Ward, T. Lähivaara, and S. Eveson, “A discontinuous Galerkin method for poroelastic wave propagation: The two-dimensional case,” Journal of Computational Physics, vol. 350, pp. 690–727, 2017, doi: 10.1016/j.jcp.2017.08.070.
      self.rho_s_arr = np.array([2200.0])
      self.rho_f_arr = np.array([900.0])
      self.K_s_arr   = np.array([5.0e+9])
      self.K_f_arr   = np.array([1.0e+9])
      self.K_b_arr   = np.array([3.0e+9])
      self.mu_b_arr  = np.array([1.0e+9])
      self.T_arr     = np.array([1.2])
      self.kappa_arr = np.array([1e-12])
      self.phi_arr   = np.array([0.4])
      self.eta_arr   = np.array([1e-3])  
    elif param==6:
      # Table 8
      # [1] N. F. Dudley Ward, T. Lähivaara, and S. Eveson, “A discontinuous Galerkin method for poroelastic wave propagation: The two-dimensional case,” Journal of Computational Physics, vol. 350, pp. 690–727, 2017, doi: 10.1016/j.jcp.2017.08.070.
      self.rho_s_arr = np.array([2200.0])
      self.rho_f_arr = np.array([950.0])
      self.K_s_arr   = np.array([7.0e+9])
      self.K_f_arr   = np.array([2.0e+9])
      self.K_b_arr   = np.array([6.5e+9])
      self.mu_b_arr  = np.array([3.0e+9])
      self.T_arr     = np.array([2.0])
      self.kappa_arr = np.array([1e-10])
      self.phi_arr   = np.array([0.2])
      self.eta_arr   = np.array([1e-3])  

    # Porous medium parameters
    self.build_parameters()


class Set3(PhysicalParameters):
  def __init__(self, **kwargs):
    super(Set3, self).__init__(**kwargs)

    # Physical parameters tuned for the experiment
    self.rho_s_arr = np.array([2600, 2600, 2600, 2600])
    self.rho_f_arr = np.array([1.29, 1000, 1000, 1.29])
    self.K_s_arr   = np.array([2.3e+8, 2.3e+8, 2.5e+8, 4.6e+8])
    self.K_f_arr   = np.array([1.4e+5, 2e+09, 2e+09, 1.4e+5])
    self.K_b_arr   = np.array([1.5e+8, 1.5e+8, 1.7e+8, 4.2e+8])
    self.mu_b_arr  = np.array([1.33e+8, 1.33e+8, 2.84e+08, 4.44e+08])
    self.T_arr     = np.array([1, 1, 1, 1])
    self.kappa_arr = np.array([8e-9, 8e-9, 1e-10, 1e-11])
    self.phi_arr   = np.array([0.2, 0.2, 0.2, 0.2])
    self.eta_arr   = np.array([0.02e-3, 1e-3, 1e-3, 0.02e-3])

    # Porous medium parameters
    self.build_parameters()


class Set4(PhysicalParameters):
  def __init__(self, **kwargs):
    super(Set4, self).__init__(**kwargs)

    # Physical parameters tuned for the experiment
    self.rho_s_arr = np.array([1850.0])
    self.rho_f_arr = np.array([1000.0])
    self.K_s_arr   = np.array([1.8e+8])
    self.K_f_arr   = np.array([2.0e+9])
    self.K_b_arr   = np.array([1.0e+8])
    self.mu_b_arr  = np.array([1.33e+8])
    self.T_arr     = np.array([1.0])
    self.kappa_arr = np.array([8e-9])
    self.phi_arr   = np.array([0.3])
    self.eta_arr   = np.array([1e-5]) 

    # Porous medium parameters
    self.build_parameters()


# Function to build the heterogeneous parameters
def LayeredParameters(param, layers, V):

  # Retrieve mesh
  mesh = V.mesh()

  # Physical parameter
  param_f = Function(V)
  for i, layer in enumerate(layers):
    try:
      param_f.vector()[layer] = param.values()[0]
    except:
      param_f.vector()[layer] = param[i]      

  return param_f


# Function to transfer the heterogeneous parameters to the hybrid meshes
def Transfer(param, V, V_RD, V_PML):

  # Retrieve meshes
  mesh = V.mesh()
  mesh_RD = V_RD.mesh()
  mesh_PML = V_PML.mesh()

  # Maps between meshes
  map_RD  = mesh_RD.topology().mapping()[mesh.id()].cell_map()
  map_PML  = mesh_PML.topology().mapping()[mesh.id()].cell_map()    

  # Physical parameters on individual meshes
  param_RD  = Function(V_RD)
  param_PML = Function(V_PML)
  param_RD.vector()[:] = param.vector()[map_RD]
  param_PML.vector()[:] = param.vector()[map_PML]

  return (param_RD, param_PML)
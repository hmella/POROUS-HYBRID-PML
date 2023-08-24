# Copyright (C) 2019 Hernan Mella

from fenics import *
from multiphenics import *
from PML_base.ConstitutiveRelations import (PoroelasticCompliance,
                                            PoroelasticStress)
from PML_base.Newmark import Newmark, NewmarkHybrid, update


class BaseProblem:
  def __init__(self, mesh=None, parameters=None, 
              g=Constant(0.0),
              TimeStepping=None, XDMF=True, HDF5=True,
              PVD=False, save_every_n=20,
              domains=None,
              boundaries=None,
              uD=None,
              wD=None,
              pD=None,
              element='P2-P1'):
    self.mesh = mesh
    self.parameters = parameters
    self.g = g
    self.TimeStepping = TimeStepping
    self.XDMF = XDMF
    self.HDF5 = HDF5
    self.PVD = PVD
    self.save_every_n = save_every_n
    self.domains = domains
    self.boundaries = boundaries
    self.uD = uD
    self.wD = wD
    self.pD = pD
    self.element = element


# Class for poroelastic problem
class Problem(BaseProblem):
  def __init__(self, **kwargs):
    super(Problem, self).__init__(**kwargs)
    self.W = self.function_space()
    self.set_problem()

  def function_space(self):
    # Create function spaces
    if self.element=="P2-P1":
      VE = VectorElement("CG", self.mesh.ufl_cell(), 2, dim=2)
      FE = FiniteElement("CG", self.mesh.ufl_cell(), 1)
    elif self.element=="P1b-P1":
      FE = FiniteElement("CG", triangle, 1)
      B = FiniteElement("B", triangle, 3)
      VE = MixedElement(2*[FE+B])
    elif self.element=="P1-P1":
      VE = VectorElement("CG", self.mesh.ufl_cell(), 1, dim=2)
      FE = FiniteElement("CG", self.mesh.ufl_cell(), 1)
    W = FunctionSpace(self.mesh, MixedElement([VE,VE,FE]))
    print("[Function space] Number of dofs: {:d}".format(W.dim()))
    return W

  def set_problem(self):
    # Time stepping parameters
    time_stepping = self.TimeStepping
    dt      = time_stepping.dt
    gamma   = time_stepping.gamma
    beta    = time_stepping.beta

    # Porous medium parameters
    rho_s = self.parameters.rho_s  # solid densitiy
    rho_f = self.parameters.rho_f  # fluid density
    K_s   = self.parameters.K_s    # bulk moduli of fluid
    K_f   = self.parameters.K_f    # bulk moduli of the solid
    K_b   = self.parameters.K_b    # bulk moduli of the dry porous frame
    mu_b  = self.parameters.mu_b   # tortuosity
    kappa = self.parameters.kappa  # permeability of the fluid
    phi   = self.parameters.phi    # porosity
    alpha = self.parameters.alpha  # Biot-Willis coefficient
    rho   = self.parameters.rho    # effective density
    rho_w = self.parameters.rho_w 
    eta   = self.parameters.eta    # shear moduli of the dry porous frame
    T     = self.parameters.T      # dynamic viscosity of the fluid
    M     = self.parameters.M      # fluid-solid coupling modulus

    # Stress and tensor
    stress = PoroelasticStress(mu_b=mu_b,K_b=K_b,alpha=alpha,M=M)

    # Set up initial values
    self.U0 = Function(self.W)  
    self.V0 = Function(self.W)
    self.A0 = Function(self.W)
    (u0, wu0, pu0) = split(self.U0)
    (v0, wv0, pv0) = split(self.V0)
    (a0, wa0, pa0) = split(self.A0)

    # Newmark accelerations and velocities
    Nu = Newmark(u0=u0, v0=v0, a0=a0, beta=beta, gamma=gamma, dt=dt)
    Nw = Newmark(u0=wu0, v0=wv0, a0=wa0, beta=beta, gamma=gamma, dt=dt)
    Np = Newmark(u0=pu0, v0=pv0, a0=pa0, beta=beta, gamma=gamma, dt=dt)

    # Create measure for the source term
    dx = Measure("dx", domain=self.mesh, subdomain_data=self.domains)
    ds = Measure("ds", domain=self.mesh, subdomain_data=self.boundaries)

    # Project boundary conditions if bubble elements are being used
    if self.element=="P1b-P1":
      P = FiniteElement("CG", self.mesh.ufl_cell(), 1)
      B = FiniteElement("B", self.mesh.ufl_cell(), 3)
      V = FunctionSpace(self.mesh, MixedElement(2*[P+B]))
      if self.uD is not None:
        for i, ud in enumerate(self.uD):
          self.uD[i] = [project(ud[0], V, solver_type="gmres"), ud[1]]
      if self.wD is not None:
        for i, wd in enumerate(self.wD):
          self.wD[i] = [project(wd[0], V, solver_type="gmres"), wd[1]]

    # Boundary conditions
    bcs = []
    if self.uD is not None:
      for ud in self.uD:
        bcs.append(DirichletBC(self.W.sub(0), ud[0], self.boundaries, ud[1]))
    if self.wD is not None:
      for wd in self.wD:
        bcs.append(DirichletBC(self.W.sub(1), wd[0], self.boundaries, wd[1]))
    if self.pD is not None:
      for pd in self.pD:
        bcs.append(DirichletBC(self.W.sub(2), pd[0], self.boundaries, pd[1]))
    self.bcs = bcs

    # Test and trial functions
    (u, w, p) = TrialFunctions(self.W)
    (v, x, q) = TestFunctions(self.W)

    # Normal vector
    n = FacetNormal(self.mesh)

    # Define variational problem
    F1 = rho*dot(Nu.ddot(u), v)*dx \
      + rho_f*dot(Nw.ddot(w), v)*dx \
      + inner(stress.tensor(u) - alpha*p*Identity(2), grad(v))*dx \
      - dot(self.g, v)*ds(3)
    F2 = -p*div(x)*dx + inner(p*n, x)*ds(3) + dot(rho_w*Nw.ddot(w), x)*dx \
      + dot(rho_f*Nu.ddot(u), x)*dx + dot(eta/kappa*Nw.dot(w), x)*dx
    F3 = div(M*(Nw.dot(w)+ alpha*Nu.dot(u)))*q*dx + Np.dot(p)*q*dx
    F = F1 + F2 + F3
    self.a, self.L = lhs(F), rhs(F)

    # Assemble rhs (once)
    self.A = PETScMatrix()
    self.b = PETScVector()
    self.A = assemble(self.a, tensor=self.A)

    # Create GMRES Krylov solver
    use_preconditioner = True
    self.solver = PETScKrylovSolver()
    PETScOptions.set("ksp_type", "gmres")
    PETScOptions.set("ksp_monitor")
    PETScOptions.set("ksp_max_it", 1500)
    if use_preconditioner:
      PETScOptions.set("pc_type", "hypre")
      PETScOptions.set('pc_hypre_type', 'euclid')
      PETScOptions.set("ksp_reuse_preconditioner","true")
      PETScOptions.set('ksp_pc_side', 'right')
    # PETScOptions.set("ksp_rtol", "1.0e-10")
    # PETScOptions.set("ksp_atol", "1.0e-12")
    PETScOptions.set("ksp_rtol", "1.0e-7")
    PETScOptions.set("ksp_atol", "1.0e-9")
    PETScOptions.set("ksp_initial_guess_nonzero", True)
    PETScOptions.set('ksp_view')
    self.solver.set_operator(self.A)
    self.solver.set_from_options()

  def solve(self):
    # MPI comm and rank
    comm = self.mesh.mpi_comm()
    rank = MPI.rank(comm)

    # Output files
    if self.XDMF:
      XDMFu = XDMFFile(comm, "output/xdmf/u_extended.xdmf")
      XDMFw = XDMFFile(comm, "output/xdmf/w_extended.xdmf")
      XDMFp = XDMFFile(comm, "output/xdmf/p_extended.xdmf")
      XDMFu.parameters["rewrite_function_mesh"] = False
      XDMFw.parameters["rewrite_function_mesh"] = False
      XDMFp.parameters["rewrite_function_mesh"] = False
    if self.HDF5:
      HDF5u = HDF5File(comm, "output/h5/u_extended.h5", "w")
      HDF5w = HDF5File(comm, "output/h5/w_extended.h5", "w")
      HDF5p = HDF5File(comm, "output/h5/p_extended.h5", "w")
    if self.PVD:
      PVDu = File("output/pvd/u_extended.pvd","compressed")
      PVDw = File("output/pvd/w_extended.pvd","compressed")
      PVDp = File("output/pvd/p_extended.pvd","compressed")

    # Discrete solution
    Uh = Function(self.W)
    Uh.rename('U','U')

    # Time stepping parameters
    t     = self.TimeStepping.t
    t_end = self.TimeStepping.t_end
    dt    = self.TimeStepping.dt
    gamma = self.TimeStepping.gamma
    beta  = self.TimeStepping.beta

    # Solving loop
    step = 0
    while t < t_end - 0.5*dt:

      t += float(dt)
      step += 1

      if rank == 0:
        print('\n\rtime: {:.3f} (Progress: {:.2f}%)'.format(t, 100*t/t_end), flush=True)

      # Update source term
      self.g.t = t

      # Assemble rhs and apply boundary condition
      self.b = assemble(self.L, tensor=self.b)
      [bc.apply(self.A, self.b) for bc in self.bcs]

      # Compute solution
      self.solver.solve(Uh.vector(), self.b)
      (u, w, p) = Uh.split(True)

      # Rename functions
      u.rename('u','u')
      w.rename('w','w')
      p.rename('p','p')

      # Update previous time step
      update(Uh, self.U0, self.V0, self.A0, beta, gamma, dt)

      # Text file for debugging
      if rank==0:
        with open('timesteps_extended.txt', 'a') as debug:
          debug.write('Saved timestep {:d} (time {:.4f} s)\n'.format(step,t))

      # Save solutions
      if self.XDMF:
        if (step % self.save_every_n == 0) or (step == 1):
          XDMFu.write(u, t)
          XDMFw.write(w, t)
          XDMFp.write(p, t)
      if self.HDF5:
        HDF5u.write(u, "u", t)
        HDF5w.write(w, "w", t)
        HDF5p.write(p, "p", t)
      if self.PVD:
        if (step % self.save_every_n == 0) or (step == 1):
          PVDu << (u, t)
          PVDw << (w, t)
          PVDp << (p, t)

    # Close files
    if self.XDMF:
      XDMFu.close()
      XDMFw.close()
      XDMFp.close()
    if self.HDF5:
      HDF5u.close()
      HDF5w.close()
      HDF5p.close()




# Poroelastic problem with paraxial conditions
class ParaxialProblem(BaseProblem):
  def __init__(self, solver_type="direct",
                     use_preconditioner=False,
                     set_nullspace=False,
                     **kwargs):
    super(ParaxialProblem, self).__init__(**kwargs)
    self.W = self.function_space()
    self.solver_type = solver_type
    self.use_preconditioner = use_preconditioner
    self.set_nullspace = set_nullspace
    self.set_problem()

  def function_space(self):
    # Create function spaces
    if self.element=="P2-P1":
      VE = VectorElement("CG", self.mesh.ufl_cell(), 2, dim=2)
      FE = FiniteElement("CG", self.mesh.ufl_cell(), 1)
    elif self.element=="P1b-P1":
      FE = FiniteElement("CG", triangle, 1)
      B = FiniteElement("B", triangle, 3)
      VE = MixedElement(2*[FE+B])
    elif self.element=="P1-P1":
      VE = VectorElement("CG", self.mesh.ufl_cell(), 1, dim=2)
      FE = FiniteElement("CG", self.mesh.ufl_cell(), 1)
    W = FunctionSpace(self.mesh, MixedElement([VE,VE,FE]))
    print("[Function space] Number of dofs: {:d}".format(W.dim()))
    return W

  def set_problem(self):
    # Time stepping parameters
    time_stepping = self.TimeStepping
    dt      = time_stepping.dt
    gamma   = time_stepping.gamma
    beta    = time_stepping.beta

    # Porous medium parameters
    rho_s = self.parameters.rho_s  # solid densitiy
    rho_f = self.parameters.rho_f  # fluid density
    K_s   = self.parameters.K_s    # bulk moduli of fluid
    K_f   = self.parameters.K_f    # bulk moduli of the solid
    K_b   = self.parameters.K_b    # bulk moduli of the dry porous frame
    mu_b  = self.parameters.mu_b   # tortuosity
    kappa = self.parameters.kappa  # permeability of the fluid
    phi   = self.parameters.phi    # porosity
    alpha = self.parameters.alpha  # Biot-Willis coefficient
    rho   = self.parameters.rho    # effective density
    rho_w = self.parameters.rho_w 
    eta   = self.parameters.eta    # shear moduli of the dry porous frame
    T     = self.parameters.T      # dynamic viscosity of the fluid
    M     = self.parameters.M      # fluid-solid coupling modulus

    # Stress and strain tensor
    stress = PoroelasticStress(mu_b=mu_b,K_b=K_b,alpha=alpha,M=M)
    strain = PoroelasticCompliance(mu_b=mu_b,K_b=K_b,alpha=alpha,M=M)

    # Set up initial values
    self.U0 = Function(self.W)  
    self.V0 = Function(self.W)
    self.A0 = Function(self.W)
    (u0, wu0, pu0) = split(self.U0)
    (v0, wv0, pv0) = split(self.V0)
    (a0, wa0, pa0) = split(self.A0)

    # Newmark accelerations and velocities
    Nu = Newmark(u0=u0, v0=v0, a0=a0, beta=beta, gamma=gamma, dt=dt)
    Nw = Newmark(u0=wu0, v0=wv0, a0=wa0, beta=beta, gamma=gamma, dt=dt)
    Np = Newmark(u0=pu0, v0=pv0, a0=pa0, beta=beta, gamma=gamma, dt=dt)

    # Create measure for the source term
    dx = Measure("dx", domain=self.mesh, subdomain_data=self.domains)
    ds = Measure("ds", domain=self.mesh, subdomain_data=self.boundaries)

    # Boundary conditions
    bcs = []
    self.bcs = bcs

    # Test and trial functions
    (u, w, p) = TrialFunctions(self.W)
    (v, x, q) = TestFunctions(self.W)

    # Paraxial boundary condition
    # T. Akiyoshi et al., “Absorbing boundary conditions for dynamic analysis of fluid-saturated porous media,” Soil Dynamics and Earthquake Engineering, vol. 13, no. 6, pp. 387–397, Jan. 1994, doi: 10.1016/0267-7261(94)90009-4.
    n  = FacetNormal(self.mesh)   # Normal vector
    nt = as_vector([n[1],-n[0]])  # Tangential vector
    Id = Identity(2)
    V1 = self.parameters.V1
    V3 = self.parameters.V3
    lmbda = K_b - 2/3*mu_b + alpha**2*M
    ts = -(lmbda + 2*mu_b + alpha**2/M)/V1*dot(Nu.dot(u), n)*n
    ts += -alpha/M/V1*dot(Nw.dot(w), n)*n
    ts += -mu_b/V3*dot(Nu.dot(u), nt)*nt
    tf = 1/(M*V1)*dot(alpha*Nu.dot(u) + Nw.dot(w), n)*n 

    # # # # C. Morency and J. Tromp, “Spectral-element simulations of wave propagation in porous media,” Geophysical Journal International, vol. 175, no. 1, pp. 301–345, 2008, doi: 10.1111/j.1365-246X.2008.03907.x.
    # # # n = FacetNormal(self.mesh)
    # # # Id = Identity(2)
    # # # nn = outer(n, n)
    # # # cp1 = self.parameters.V1
    # # # cp2 = self.parameters.V2
    # # # cs  = self.parameters.V3
    # # # rho = self.parameters.rho
    # # # ts = -rho*cp1*dot(Nu.dot(u), n)*n
    # # # ts += -rho_f*cp2*dot(Nw.dot(w), n)*n
    # # # ts += -(rho - rho_f*phi/T)*cs*(Id - nn)*Nu.dot(u)
    # # # tf = -rho_f*T/phi*cp2*dot(Nw.dot(w), n)*n
    # # # tf += -rho_f*cp1*dot(Nu.dot(u), n)*n

    # Define variational problem
    F1 = rho*inner(Nu.ddot(u), v)*dx \
      + rho_f*inner(Nw.ddot(w), v)*dx \
      + inner(stress.tensor(u) - alpha*p*Id, grad(v))*dx \
      - inner(self.g, v)*ds(3) \
      - inner(ts, v)*ds(1)
    F2 = -inner(p, div(x))*dx \
      + inner(tf, x)*ds(1) \
      + inner(p*n, x)*ds(3) \
      + inner(rho_f*Nu.ddot(u) + rho_w*Nw.ddot(w) + eta/kappa*Nw.dot(w), x)*dx  
    F3 = (div(M*alpha*Nu.dot(u) + M*Nw.dot(w)) + Np.dot(p))*q*dx
    F = F1 + F2 + F3
    self.a, self.L = lhs(F), rhs(F)

    # Assemble rhs (once)
    self.A = PETScMatrix()
    self.b = PETScVector()
    self.A = assemble(self.a, tensor=self.A)

    # Create GMRES Krylov solver
    if self.solver_type == "iterative":
      self.solver = PETScKrylovSolver()
      PETScOptions.set("ksp_type", "fgmres")
      PETScOptions.set("ksp_max_it", 1500)
      PETScOptions.set("ksp_monitor_true_residual")
      PETScOptions.set("ksp_rtol", "1.0e-7")
      PETScOptions.set("ksp_atol", "1.0e-9")
      # PETScOptions.set("ksp_rtol", "1.0e-15")
      # PETScOptions.set("ksp_atol", "1.0e-17")
      PETScOptions.set("ksp_initial_guess_nonzero", True)
      self.solver.set_from_options()
      self.solver.set_operator(self.A)
    else:
      self.solver = PETScLUSolver('mumps')
      PETScOptions.set("reuse_factorization",True)
      self.solver.set_from_options()
      self.solver.set_operator(self.A)

    # Remove rigid body motions (overwrites previously defined
    # solver if needed)
    if self.set_nullspace:
      # Create near null space basis (required for smoothed aggregation
      # AMG). The solution vector is passed so that it can be copied to
      # generate compatible vectors for the nullspace.
      null_space = self.build_nullspace(self.W, self.U0.vector())

      # Attach near nullspace to matrix
      as_backend_type(self.A).set_near_nullspace(null_space)

      # Create PETSC smoothed aggregation AMG preconditioner and attach near
      # null space
      pc = PETScPreconditioner("petsc_amg")

      # Use Chebyshev smoothing for multigrid
      PETScOptions.set("mg_levels_ksp_type", "chebyshev")
      PETScOptions.set("mg_levels_pc_type", "jacobi")

      # Improve estimate of eigenvalues for Chebyshev smoothing
      PETScOptions.set("mg_levels_esteig_ksp_type", "gmres")
      PETScOptions.set("mg_levels_ksp_chebyshev_esteig_steps", 50)

      # Create CG Krylov solver and turn convergence monitoring on
      self.solver = PETScKrylovSolver("gmres", pc)
      self.solver.set_operator(self.A)

  def build_nullspace(self, V, x):
      """Function to build null space for 3D elasticity"""

      # Create list of vectors for null space
      nullspace_basis = [x.copy() for i in range(3)]

      # Build translational null space basis
      V.sub(0).sub(0).dofmap().set(nullspace_basis[0], 1.0)
      V.sub(0).sub(1).dofmap().set(nullspace_basis[1], 1.0)

      # Build rotational null space basis
      V.sub(0).sub(0).set_x(nullspace_basis[2], -1.0, 1)
      V.sub(0).sub(1).set_x(nullspace_basis[2],  1.0, 0)

      for x in nullspace_basis:
          x.apply("insert")

      # Create vector space basis and orthogonalize
      basis = VectorSpaceBasis(nullspace_basis)
      basis.orthonormalize()

      return basis

  def solve(self):
    # MPI comm and rank
    comm = self.mesh.mpi_comm()
    rank = MPI.rank(comm)

    # Output files
    if self.XDMF:
      XDMFu = XDMFFile(comm, "output/xdmf/u_paraxial.xdmf")
      XDMFw = XDMFFile(comm, "output/xdmf/w_paraxial.xdmf")
      XDMFp = XDMFFile(comm, "output/xdmf/p_paraxial.xdmf")
      XDMFu.parameters["rewrite_function_mesh"] = False
      XDMFw.parameters["rewrite_function_mesh"] = False
      XDMFp.parameters["rewrite_function_mesh"] = False
    if self.HDF5:
      HDF5u = HDF5File(comm, "output/h5/u_paraxial.h5", "w")
      HDF5w = HDF5File(comm, "output/h5/w_paraxial.h5", "w")
      HDF5p = HDF5File(comm, "output/h5/p_paraxial.h5", "w")
    if self.PVD:
      PVDu = File("output/pvd/u_paraxial.pvd","compressed")
      PVDw = File("output/pvd/w_paraxial.pvd","compressed")
      PVDp = File("output/pvd/p_paraxial.pvd","compressed")

    # Discrete solution
    Uh = Function(self.W)
    Uh.rename('U','U')

    # Time stepping parameters
    t     = self.TimeStepping.t
    t_end = self.TimeStepping.t_end
    dt    = self.TimeStepping.dt
    gamma = self.TimeStepping.gamma
    beta  = self.TimeStepping.beta

    # Solving loop
    step = 0
    while t < t_end - 0.5*dt:

      t += float(dt)
      step += 1

      if rank == 0:
        print('\n\rtime: {:.3f} (Progress: {:.2f}%)'.format(t, 100*t/t_end), flush=True)

      # Update source term
      self.g.t = t

      # Assemble rhs and apply boundary condition
      self.b = assemble(self.L, tensor=self.b)
      [bc.apply(self.A, self.b) for bc in self.bcs]

      # Compute solution
      self.solver.solve(Uh.vector(), self.b)
      (u, w, p) = Uh.split(True)

      # Rename functions
      u.rename('u','u')
      w.rename('w','w')
      p.rename('p','p')

      # Update previous time step
      update(Uh, self.U0, self.V0, self.A0, beta, gamma, dt)

      # Text file for debugging
      if rank==0:
        with open('timesteps_paraxial.txt', 'a') as debug:
          debug.write('Saved timestep {:d} (time {:.4f} s)\n'.format(step,t))

      # Save solutions
      if self.XDMF:
        if (step % self.save_every_n == 0) or (step == 1):
          XDMFu.write(u, t)
          XDMFw.write(w, t)
          XDMFp.write(p, t)
      if self.HDF5:
        HDF5u.write(u, "u", t)
        HDF5w.write(w, "w", t)
        HDF5p.write(p, "p", t)
      if self.PVD:
        if (step % self.save_every_n == 0) or (step == 1):
          PVDu << (u, t)
          PVDw << (w, t)
          PVDp << (p, t)

    # Close files
    if self.XDMF:
      XDMFu.close()
      XDMFw.close()
      XDMFp.close()
    if self.HDF5:
      HDF5u.close()
      HDF5w.close()
      HDF5p.close()




# Poroelastic problem with PML (mixed form)
class MixedPMLProblem(BaseProblem):
  def __init__(self, PMLFunctions=None,
              solver_type='direct',
              use_preconditioner=False,**kwargs):
    super(MixedPMLProblem, self).__init__(**kwargs)
    self.PMLFunctions = PMLFunctions
    self.W = self.function_space()
    self.solver_type = solver_type
    self.use_preconditioner = use_preconditioner
    self.set_problem()

  def function_space(self):
    # Cell
    triangle = self.mesh.ufl_cell()

    # Create function spaces
    if self.element=='P2-P1':
        VE = VectorElement("CG", triangle, 2, dim=2)
        FE = FiniteElement("CG", triangle, 1)
        TE = TensorElement("DG", triangle, 1, shape=(2,2), symmetry=True)
    elif self.element=='P1b-P1':
        FE = FiniteElement("CG", triangle, 1)
        B = FiniteElement("B", triangle, 3)
        VE = MixedElement(2*[FE+B])
        TE = TensorElement("DG", triangle, 1, shape=(2,2), symmetry=True)
    elif self.element=='P1-P1':
        FE = FiniteElement("CG", triangle, 1)
        VE = VectorElement("CG", triangle, 1, dim=2)
        TE = TensorElement("DG", triangle, 0, shape=(2,2), symmetry=True)
    W = FunctionSpace(self.mesh, MixedElement([VE, VE, FE, TE]))
    print("[Function space] Number of dofs: {:d}".format(W.dim()))
    return W

  def set_problem(self):
    # Time stepping parameters
    time_stepping = self.TimeStepping
    dt      = time_stepping.dt
    gamma   = time_stepping.gamma
    beta    = time_stepping.beta

    # Porous medium parameters
    rho_s = self.parameters.rho_s  # solid densitiy
    rho_f = self.parameters.rho_f  # fluid density
    K_s   = self.parameters.K_s    # bulk moduli of fluid
    K_f   = self.parameters.K_f    # bulk moduli of the solid
    K_b   = self.parameters.K_b    # bulk moduli of the dry porous frame
    mu_b  = self.parameters.mu_b   # tortuosity
    kappa = self.parameters.kappa  # permeability of the fluid
    phi   = self.parameters.phi    # porosity
    alpha = self.parameters.alpha  # Biot-Willis coefficient
    rho   = self.parameters.rho    # effective density
    rho_w = self.parameters.rho_w 
    eta   = self.parameters.eta    # shear moduli of the dry porous frame
    T     = self.parameters.T      # dynamic viscosity of the fluid
    M     = self.parameters.M      # fluid-solid coupling modulus

    # Stress and tensor
    stress = PoroelasticStress(mu_b=mu_b,K_b=K_b,alpha=alpha,M=M)
    strain = PoroelasticCompliance(mu_b=mu_b,K_b=K_b,alpha=alpha,M=M)

    # Attenuation and stretching functions and tensors
    alpha_1 = self.PMLFunctions[0]
    alpha_2 = self.PMLFunctions[1]
    beta_1 = self.PMLFunctions[2]
    beta_2 = self.PMLFunctions[3]
    a_ = alpha_1*alpha_2
    b_ = alpha_1*beta_2 + alpha_2*beta_1
    c_ = beta_1*beta_2

    Le = as_tensor([[alpha_2, 0],[0, alpha_1]])
    Lp = as_tensor([[beta_2, 0],[0, beta_1]])
    Lee = as_tensor([[alpha_1, 0],[0, alpha_2]])
    Lpp = as_tensor([[beta_1, 0],[0, beta_2]])

    # Create measure for the source term
    dx = Measure("dx", domain=self.mesh, subdomain_data=self.domains)
    ds = Measure("ds", domain=self.mesh, subdomain_data=self.boundaries)

    # Project boundary conditions if bubble elements are being used
    if self.element=="P1b-P1":
      P = FiniteElement("CG", self.mesh.ufl_cell(), 1)
      B = FiniteElement("B", self.mesh.ufl_cell(), 3)
      V = FunctionSpace(self.mesh, MixedElement(2*[P+B]))
      if self.uD is not None:
        for i, ud in enumerate(self.uD):
          self.uD[i] = [project(ud[0], V, solver_type="gmres"), ud[1]]
      if self.wD is not None:
        for i, wd in enumerate(self.wD):
          self.wD[i] = [project(wd[0], V, solver_type="gmres"), wd[1]]

    # Boundary conditions
    bcs = []
    if self.uD is not None:
      for ud in self.uD:
        bcs.append(DirichletBC(self.W.sub(0), ud[0], self.boundaries, ud[1]))
    if self.wD is not None:
      for wd in self.wD:
        bcs.append(DirichletBC(self.W.sub(1), wd[0], self.boundaries, wd[1]))
    if self.pD is not None:
      for pd in self.pD:
        bcs.append(DirichletBC(self.W.sub(2), pd[0], self.boundaries, pd[1]))
    self.bcs = bcs

    # Set up initial values
    self.U0m = Function(self.W)  
    self.V0m = Function(self.W)
    self.A0m = Function(self.W)
    (u0, wu0, piu0, U0) = split(self.U0m)
    (v0, wv0, piv0, V0) = split(self.V0m)
    (a0, wa0, pia0, A0) = split(self.A0m)

    # Newmark accelerations and velocities
    Nu = NewmarkHybrid(u0=u0, v0=v0, a0=a0, beta=beta, gamma=gamma, dt=dt)
    Nw = NewmarkHybrid(u0=wu0, v0=wv0, a0=wa0, beta=beta, gamma=gamma, dt=dt)
    Npi = NewmarkHybrid(u0=piu0, v0=piv0, a0=pia0, beta=beta, gamma=gamma, dt=dt)  
    NS = NewmarkHybrid(u0=U0, v0=V0, a0=A0, beta=beta, gamma=gamma, dt=dt)

    # Create measure for the source term
    dx = Measure("dx", domain=self.mesh, subdomain_data=self.domains)
    ds = Measure("ds", domain=self.mesh, subdomain_data=self.boundaries)

    # Function and test functions
    (u, w, pi, S) = TrialFunctions(self.W)
    (du, dw, dpi, dS) = TestFunctions(self.W)

    # Normal vector and identity matrix
    n  = FacetNormal(self.mesh)
    Id = Identity(2)

    # Define variational problem
    a = rho*inner(a_*Nu.ddotl(u) + b_*Nu.dotl(u) + c_*u, du)*dx \
      + rho_f*inner(a_*Nw.ddotl(w) + b_*Nw.dotl(w) + c_*w, du)*dx \
      + inner((NS.dotl(S) - alpha*Npi.dotl(pi)*Id)*Le + (S - alpha*pi*Id)*Lp, grad(du))*dx
    L = rho*inner(a_*Nu.ddotr(u) + b_*Nu.dotr(u), du)*dx \
      + rho_f*inner(a_*Nw.ddotr(w) + b_*Nw.dotr(w), du)*dx \
      + inner((NS.dotr(S) - alpha*Npi.dotr(pi)*Id)*Le, grad(du))*dx \
      + dot(self.g, du)*ds(3)

    a += -Npi.dotl(pi)*div(dw)*dx + inner(Npi.dotl(pi)*n, dw)*ds(3) \
      + inner(Lee*(rho_f*Nu.ddotl(u) + rho_w*Nw.ddotl(w) + eta/kappa*Nw.dotl(w)), dw)*dx \
      + inner(Lpp*(rho_f*Nu.dotl(u) + rho_w*Nw.dotl(w) + eta/kappa*w), dw)*dx
    L += -Npi.dotr(pi)*div(dw)*dx + inner(Npi.dotr(pi)*n, dw)*ds(3) \
      + inner(Lee*(rho_f*Nu.ddotr(u) + rho_w*Nw.ddotr(w) + eta/kappa*Nw.dotr(w)), dw)*dx \
      + inner(Lpp*(rho_f*Nu.dotr(u) + rho_w*Nw.dotr(w)), dw)*dx  

    a += (a_*Npi.ddotl(pi) + b_*Npi.dotl(pi) + c_*pi)*dpi*dx \
      + div(Le*M*(alpha*Nu.dotl(u) + Nw.dotl(w)))*dpi*dx \
      + div(Lp*M*(alpha*u + w))*dpi*dx
    L += (a_*Npi.ddotr(pi) + b_*Npi.dotr(pi))*dpi*dx \
      + div(Le*M*(alpha*Nu.dotr(u) + Nw.dotr(w)))*dpi*dx

    a += inner(strain.tensor(a_*NS.ddotl(S) + b_*NS.dotl(S) + c_*S), dS)*dx \
      - 0.5*inner(grad(u)*Lp + Lp*grad(u).T, dS)*dx \
      - 0.5*inner(grad(Nu.dotl(u))*Le + Le*grad(Nu.dotl(u)).T, dS)*dx
    L += inner(strain.tensor(a_*NS.ddotr(S) + b_*NS.dotr(S)), dS)*dx \
      - 0.5*inner(grad(Nu.dotr(u))*Le + Le*grad(Nu.dotr(u)).T, dS)*dx
    self.a, self.L = a, L

    # Assemble rhs (once)
    self.A = PETScMatrix()
    self.b = PETScVector()
    self.A = assemble(self.a, tensor=self.A)

    # Create GMRES Krylov solver
    if self.solver_type is 'iterative':
      self.solver = PETScKrylovSolver()
      PETScOptions.set("ksp_type", "gmres")
      PETScOptions.set("ksp_monitor")
      PETScOptions.set("ksp_max_it", 250)
      if self.use_preconditioner:
        PETScOptions.set("pc_type", "hypre")
        PETScOptions.set('pc_hypre_type', 'euclid')
        PETScOptions.set("ksp_reuse_preconditioner","true")
        PETScOptions.set('ksp_pc_side', 'right')
      PETScOptions.set("ksp_rtol", "1.0e-10")
      PETScOptions.set("ksp_atol", "1.0e-12")
      PETScOptions.set("ksp_initial_guess_nonzero",True)
      self.solver.set_from_options()
      self.solver.set_operator(self.A)
    elif self.solver_type is 'direct':
      self.solver = PETScLUSolver('mumps')
      PETScOptions.set("reuse_factorization",True)
      self.solver.set_from_options()
      self.solver.set_operator(self.A)

  def solve(self, output_suffix=''):
    # MPI comm and rank
    comm = self.mesh.mpi_comm()
    rank = MPI.rank(comm)

    # Output files
    if self.XDMF:
      tmp = "output/xdmf/{:s}_{:s}{:s}.xdmf"
      XDMFu = XDMFFile(comm, tmp.format('u','MIXED',output_suffix))
      XDMFw = XDMFFile(comm, tmp.format('w','MIXED',output_suffix))
      XDMFp = XDMFFile(comm, tmp.format('p','MIXED',output_suffix))
      XDMFu.parameters["rewrite_function_mesh"] = False
      XDMFw.parameters["rewrite_function_mesh"] = False
      XDMFp.parameters["rewrite_function_mesh"] = False
    if self.HDF5:
      tmp = "output/h5/{:s}_{:s}{:s}.h5"
      HDF5u = HDF5File(comm, tmp.format('u','MIXED',output_suffix), "w")
      HDF5w = HDF5File(comm, tmp.format('w','MIXED',output_suffix), "w")
      HDF5p = HDF5File(comm, tmp.format('p','MIXED',output_suffix), "w")
    if self.PVD:
      tmp = "output/pvd/{:s}_{:s}{:s}.pvd"
      PVDu = File(tmp.format('u','MIXED',output_suffix),"compressed")
      PVDw = File(tmp.format('w','MIXED',output_suffix),"compressed")
      PVDp = File(tmp.format('p','MIXED',output_suffix),"compressed")

    # Discrete solution
    Uh = Function(self.W)
    Uh.rename('U','U')

    # Time stepping parameters
    t     = self.TimeStepping.t
    t_end = self.TimeStepping.t_end
    dt    = self.TimeStepping.dt
    gamma = self.TimeStepping.gamma
    beta  = self.TimeStepping.beta

    # Solving loop
    step = 0
    while t < t_end - 0.5*dt:

      t += float(dt)
      step += 1

      if rank == 0:
        print('\n\rtime: {:.3f} (Progress: {:.2f}%)'.format(t, 100*t/t_end), flush=True)

      # Update source term
      self.g.t = t

      # Assemble rhs and apply boundary condition
      self.b = assemble(self.L, tensor=self.b)
      [bc.apply(self.A, self.b) for bc in self.bcs]

      # Compute solution
      self.solver.solve(Uh.vector(), self.b)
      (u, w, pi, S) = Uh.split(True)

      # Rename functions
      u.rename('u','u')
      w.rename('w','w')

      # Update previous time step
      update(Uh, self.U0m, self.V0m, self.A0m, beta, gamma, dt)

      # Get current pressure on the PML layer
      (_, _, p, _) = self.V0m.split(True)
      p.rename('p','p')

      # Text file for debugging
      if rank==0:
        with open('timesteps_mixed.txt', 'a') as debug:
          debug.write('Saved timestep {:d} (time {:.4f} s)\n'.format(step,t))

      # Save solutions
      if self.XDMF:
        if (step % self.save_every_n == 0) or (step == 1):
          XDMFu.write(u, t)
          XDMFw.write(w, t)
          XDMFp.write(p, t)
      if self.HDF5:
        HDF5u.write(u, "u", t)
        HDF5w.write(w, "w", t)
        HDF5p.write(p, "p", t)
      if self.PVD:
        if (step % self.save_every_n == 0) or (step == 1):
          PVDu << (u, t)
          PVDw << (w, t)
          PVDp << (p, t)

    # Close files
    if self.XDMF:
      XDMFu.close()
      XDMFw.close()
      XDMFp.close()
    if self.HDF5:
      HDF5u.close()
      HDF5w.close()
      HDF5p.close()




class HybridPMLProblem(BaseProblem):
  def __init__(self, PMLFunctions=None,
                     solver_type="direct",
                     use_preconditioner=True,
                     **kwargs):
    super(HybridPMLProblem, self).__init__(**kwargs)
    self.PMLFunctions = PMLFunctions
    self.solver_type = solver_type
    self.use_preconditioner = use_preconditioner
    self.W = self.function_space()
    self.set_problem()

  def function_space(self):
    # Restrictions
    RD  = MeshRestriction(self.mesh, "mesh/mesh_RD_restriction.rtc.xdmf")
    PML = MeshRestriction(self.mesh, "mesh/mesh_PML_restriction.rtc.xdmf")
    I = MeshRestriction(self.mesh, "mesh/mesh_interface_restriction.rtc.xdmf")

    # Mesh cell
    triangle = self.mesh.ufl_cell()

    # Create function spaces
    if self.element=='P2-P1':
      VE = VectorElement("CG", triangle, 2, dim=2)
      FE = FiniteElement("CG", triangle, 1)
      TE = TensorElement("DG", triangle, 1, shape=(2,2), symmetry=True)
    elif self.element=='P1b-P1':
      FE = FiniteElement("CG", triangle, 1)
      B = FiniteElement("B", triangle, 3)
      VE = MixedElement(2*[FE+B])
      TE = TensorElement("DG", triangle, 1, shape=(2,2), symmetry=True)
    elif self.element=='P1-P1':
      VE = VectorElement("CG", triangle, 1, dim=2)
      FE = FiniteElement("CG", triangle, 1)
      TE = TensorElement("DG", triangle, 0, shape=(2,2), symmetry=True)
    P = FunctionSpace(self.mesh, FE)
    V = FunctionSpace(self.mesh, VE)
    T = FunctionSpace(self.mesh, TE)
    W = BlockFunctionSpace([V,V,P,P,T,P], restrict=[None,None,RD,PML,PML,I])
    print("[Function space] Number of dofs: {:d}".format(W.dim()))
    return W

  def set_problem(self):

    # Time stepping parameters
    dt    = self.TimeStepping.dt
    gamma = self.TimeStepping.gamma
    beta  = self.TimeStepping.beta

    # Function space for physical parameters
    P = FunctionSpace(self.mesh, 'DG', 0)

    # Porous medium parameters
    rho_s = self.parameters.rho_s  # solid densitiy
    rho_f = self.parameters.rho_f  # fluid density
    K_s   = self.parameters.K_s    # bulk moduli of fluid
    K_f   = self.parameters.K_f    # bulk moduli of the solid
    K_b   = self.parameters.K_b    # bulk moduli of the dry porous frame
    mu_b  = self.parameters.mu_b   # tortuosity
    kappa = self.parameters.kappa  # permeability of the fluid
    phi   = self.parameters.phi    # porosity
    alpha = self.parameters.alpha  # Biot-Willis coefficient
    rho   = self.parameters.rho    # effective density
    rho_w = self.parameters.rho_w
    eta   = self.parameters.eta    # shear moduli of the dry porous frame
    T     = self.parameters.T      # dynamic viscosity of the fluid
    M     = self.parameters.M      # fluid-solid coupling modulus
    File('output/debug/K_b.pvd') << K_b
    File('output/debug/mu_b.pvd') << mu_b
    File('output/debug/c1p.pvd') << self.parameters.c1p
    File('output/debug/cs.pvd') << self.parameters.cs

    # Stress and strain tensors
    stress = PoroelasticStress(mu_b=mu_b,K_b=K_b,alpha=alpha,M=M)
    strain = PoroelasticCompliance(mu_b=mu_b,K_b=K_b,alpha=alpha,M=M)

    # Attenuation and stretching functions and tensors
    alpha_1 = self.PMLFunctions[0]
    alpha_2 = self.PMLFunctions[1]
    beta_1 = self.PMLFunctions[2]
    beta_2 = self.PMLFunctions[3]
    a_ = alpha_1*alpha_2
    b_ = alpha_1*beta_2 + alpha_2*beta_1
    c_ = beta_1*beta_2

    Le = as_tensor([[alpha_2, 0],[0, alpha_1]])
    Lp = as_tensor([[beta_2, 0],[0, beta_1]])
    Lee = as_tensor([[alpha_1, 0],[0, alpha_2]])
    Lpp = as_tensor([[beta_1, 0],[0, beta_2]])

    # Create measure for the source term
    dx = Measure("dx")(subdomain_data=self.domains)
    ds = Measure("ds")(subdomain_data=self.boundaries)
    dS = Measure("dS")(subdomain_data=self.boundaries)

    # Project boundary conditions if bubble elements are being used
    if self.element=="P1b-P1":
      P = FiniteElement("CG", self.mesh.ufl_cell(), 1)
      B = FiniteElement("B", self.mesh.ufl_cell(), 3)
      V = FunctionSpace(self.mesh, MixedElement(2*[P+B]))
      if self.uD is not None:
        for i, ud in enumerate(self.uD):
          self.uD[i] = [project(ud[0], V, solver_type="gmres"), ud[1]]
      if self.wD is not None:
        for i, wd in enumerate(self.wD):
          self.wD[i] = [project(wd[0], V, solver_type="gmres"), wd[1]]

    # Boundary conditions
    bcs = []
    if self.uD is not None:
      for ud in self.uD:
        bcs.append(DirichletBC(self.W.sub(0), ud[0], self.boundaries, ud[1]))
    if self.wD is not None:
      for wd in self.wD:
        bcs.append(DirichletBC(self.W.sub(1), wd[0], self.boundaries, wd[1]))
    if self.pD is not None:
      for pd in self.pD:
        bcs.append(DirichletBC(self.W.sub(3), pd[0], self.boundaries, pd[1]))
    self.bcs = BlockDirichletBC(bcs)

    # Set up initial values
    self.U0 = BlockFunction(self.W)
    (self.u0, self.wu0, self.pu0, self.piu0, self.U0, _) = block_split(self.U0)

    self.V0 = BlockFunction(self.W)
    (self.v0, self.wv0, self.pv0, self.piv0, self.V0, _) = block_split(self.V0)

    self.A0 = BlockFunction(self.W)
    (self.a0, self.wa0, self.pa0, self.pia0, self.A0, _) = block_split(self.A0)

    # Function and test functions
    U_block = BlockTrialFunction(self.W)
    (u, w, p, pi, S, lp) = block_split(U_block)
    W_block = BlockTestFunction(self.W)
    (du, dw, dp, dpi, dT, ep) = block_split(W_block)

    # Newmark velocities and accelerations
    Nu = NewmarkHybrid(u0=self.u0, v0=self.v0, a0=self.a0, beta=beta, gamma=gamma, dt=dt)
    Np  = NewmarkHybrid(u0=self.pu0, v0=self.pv0, a0=self.pa0, beta=beta, gamma=gamma, dt=dt)
    Nw = NewmarkHybrid(u0=self.wu0, v0=self.wv0, a0=self.wa0, beta=beta, gamma=gamma, dt=dt)
    Npi = NewmarkHybrid(u0=self.piu0, v0=self.piv0, a0=self.pia0, beta=beta, gamma=gamma, dt=dt)  
    NS = NewmarkHybrid(u0=self.U0, v0=self.V0, a0=self.A0, beta=beta, gamma=gamma, dt=dt)

    # Normal vector and identity matrix
    n = FacetNormal(self.mesh)
    Id = Identity(2)

    # Variational problem
    a11 = rho*inner(Nu.ddotl(u), du)*dx(1) \
        + inner(stress.tensor(u), grad(du))*dx(1)
    a12 = rho_f*inner(Nw.ddotl(w), du)*dx(1)
    a13 = inner(-alpha*p*Id, grad(du))*dx(1)

    a11 += rho*inner(a_*Nu.ddotl(u) + b_*Nu.dotl(u) + c_*u, du)*dx(2)
    a12 += rho_f*inner(a_*Nw.ddotl(w) + b_*Nw.dotl(w) + c_*w, du)*dx(2)
    a14 = inner(-alpha*Npi.dotl(pi)*Id*Le - alpha*pi*Id*Lp, grad(du))*dx(2)
    a15 = inner(NS.dotl(S)*Le + S*Lp, grad(du))*dx(2)

    a21 = inner(rho_f*Nu.ddotl(u), dw)*dx(1)
    a22 = inner(rho_w*Nw.ddotl(w) + eta/kappa*Nw.dotl(w), dw)*dx(1)
    a23 = -p*div(dw)*dx(1) + dot(p*n, dw)*ds(3)

    a21 += inner(Lee*(rho_f*Nu.ddotl(u)), dw)*dx(2) \
         + inner(Lpp*(rho_f*Nu.dotl(u)), dw)*dx(2)
    a22 += inner(Lee*(rho_w*Nw.ddotl(w) + eta/kappa*Nw.dotl(w)), dw)*dx(2) \
         + inner(Lpp*(rho_w*Nw.dotl(w) + eta/kappa*w), dw)*dx(2)
    a24 = -Npi.dotl(pi)*div(dw)*dx(2)

    a31 = div(M*alpha*Nu.dotl(u))*dp*dx(1)
    a32 = div(M*Nw.dotl(w))*dp*dx(1)
    a33 = Np.dotl(p)*dp*dx(1)

    # a41 = div(Le*M*alpha*Nu.dotl(u))*dpi*dx(2) \
    #     + div(Lp*M*alpha*u)*dpi*dx(2)
    # a42 = div(Le*M*Nw.dotl(w))*dpi*dx(2) \
    #     + div(Lp*M*w)*dpi*dx(2)
    a41 = inner(Le, grad(M*alpha*Nu.dotl(u)))*dpi*dx(2) \
        + inner(Lp, grad(M*alpha*u))*dpi*dx(2)
    a42 = inner(Le, grad(M*Nw.dotl(w)))*dpi*dx(2) \
        + inner(Lp, grad(M*w))*dpi*dx(2)
    a44 = (a_*Npi.ddotl(pi) + b_*Npi.dotl(pi) + c_*pi)*dpi*dx(2)

    a51 = -0.5*inner(grad(u)*Lp + Lp*grad(u).T, dT)*dx(2) \
        - 0.5*inner(grad(Nu.dotl(u))*Le + Le*grad(Nu.dotl(u)).T, dT)*dx(2)
    a55 = inner(strain.tensor(a_*NS.ddotl(S) + b_*NS.dotl(S) + c_*S), dT)*dx(2)

    L1 = rho*inner(Nu.ddotr(u), du)*dx(1) \
       + rho_f*inner(Nw.ddotr(w), du)*dx(1) \
       + inner(self.g, du)*ds(3)
    L1 += rho*inner(a_*Nu.ddotr(u) + b_*Nu.dotr(u), du)*dx(2) \
       + rho_f*inner(a_*Nw.ddotr(w) + b_*Nw.dotr(w), du)*dx(2) \
       + inner((NS.dotr(S) - alpha*Npi.dotr(pi)*Id)*Le, grad(du))*dx(2)
    L2 = inner(rho_f*Nu.ddotr(u) + rho_w*Nw.ddotr(w) + eta/kappa*Nw.dotr(w), dw)*dx(1)
    L2 += -Npi.dotr(pi)*div(dw)*dx(2) \
       + inner(Lee*(rho_f*Nu.ddotr(u) + rho_w*Nw.ddotr(w) + eta/kappa*Nw.dotr(w)), dw)*dx(2) \
       + inner(Lpp*(rho_f*Nu.dotr(u) + rho_w*Nw.dotr(w)), dw)*dx(2)
    L3 = div(M*(Nw.dotr(w) + alpha*Nu.dotr(u)))*dp*dx(1) \
       + Np.dotr(p)*dp*dx(1)
    # L4 = div(Le*M*(alpha*Nu.dotr(u) + Nw.dotr(w)))*dpi*dx(2) \
    #    + (a_*Npi.ddotr(pi) + b_*Npi.dotr(pi))*dpi*dx(2)
    L4 = inner(Le, grad(M*(alpha*Nu.dotr(u) + Nw.dotr(w))))*dpi*dx(2) \
       + (a_*Npi.ddotr(pi) + b_*Npi.dotr(pi))*dpi*dx(2)
    L5 = inner(strain.tensor(a_*NS.ddotr(S) + b_*NS.dotr(S)), dT)*dx(2) \
       - 0.5*inner(grad(Nu.dotr(u))*Le + Le*grad(Nu.dotr(u)).T, dT)*dx(2)

    a36 = -dot(avg(lp), dp('+'))*dS(2)
    a46 = dot(avg(lp), dpi('-'))*dS(2)
    a63 = -dot(p('+'), avg(ep))*dS(2)
    a64 = dot(Npi.dotl(pi('-')), avg(ep))*dS(2)
    L6 = dot(Npi.dotr(pi('-')), avg(ep))*dS(2)


    # Left and right hand sides
    self.rhs = [L1,L2,L3,L4,L5,L6]
    self.lhs = [[a11, a12, a13, a14, a15, 0],
                [a21, a22, a23, a24, 0,   0],
                [a31, a32, a33, 0,   0,   a36],
                [a41, a42, 0,   a44, 0,   a46],
                [a51, 0,   0,   0,   a55, 0],
                [0,   0,   a63, a64, 0,   0]]

    # Assemble LHS and apply boundary conditions
    self.A = block_assemble(self.lhs)
    self.bcs.apply(self.A)

    # Solver
    if self.solver_type=="iterative":
      self.solver = PETScKrylovSolver()
      PETScOptions.set("ksp_type", "fgmres")
      PETScOptions.set("ksp_max_it", 1500)
      PETScOptions.set("ksp_monitor_true_residual")
      PETScOptions.set("ksp_rtol", "1.0e-7")
      PETScOptions.set("ksp_atol", "1.0e-9")
      # PETScOptions.set("ksp_rtol", "1.0e-15")
      # PETScOptions.set("ksp_atol", "1.0e-17")
      PETScOptions.set("ksp_initial_guess_nonzero", True)
      if self.use_preconditioner==True:
          # # PETScOptions.set("pctype", "fieldsplit")
          # # PETScOptions.set("pc_fieldsplit_type", "schur")
          # # PETScOptions.set("pc_fieldsplit_schur_fact_type", "full")
          # # PETScOptions.set('pc_fieldsplit_detect_saddle_point')
          # # PETScOptions.set("fieldsplit_0_ksp_type", "cg")
          # # PETScOptions.set("fieldsplit_0_pc_type", "ilu")
          # # PETScOptions.set("fieldsplit_0_ksp_rtol", 1e-12)
          # # PETScOptions.set("fieldsplit_1_ksp_type", "cg")
          # # PETScOptions.set("fieldsplit_1_pc_type", "none")
          # # PETScOptions.set("fieldsplit_1_ksp_rtol", 1e-12)
          PETScOptions.set("pctype", "fieldsplit")
          PETScOptions.set("pc_fieldsplit_type", "schur")
          PETScOptions.set("pc_fieldsplit_schur_fact_type", "full")
          PETScOptions.set('pc_fieldsplit_detect_saddle_point')
          PETScOptions.set("fieldsplit_0_ksp_type", "cg")
          PETScOptions.set("fieldsplit_0_pc_type", "ilu")
          PETScOptions.set("fieldsplit_0_ksp_rtol", 1e-12)
          PETScOptions.set("fieldsplit_1_ksp_type", "cg")
          PETScOptions.set("fieldsplit_1_ksp_rtol", 1e-12)
          PETScOptions.set("pc_fieldsplit_schur_precondition", "selfp")
          PETScOptions.set("fieldsplit_1_pc_type", "hypre")
          PETScOptions.set("fieldsplit_1_pc_hypre_type", "euclid")
      self.solver.set_from_options()
      self.solver.set_operator(self.A)
    else:
      self.solver = PETScLUSolver('mumps')
      PETScOptions.set("reuse_factorization",True)
      self.solver.set_from_options()
      self.solver.set_operator(self.A)

  def solve(self, output_suffix=''):
    # MPI comm and rank
    comm = self.mesh.mpi_comm()
    rank = MPI.rank(comm)

    # Output files
    if self.XDMF:
      tmp = "output/xdmf/{:s}_{:s}{:s}.xdmf"
      XDMFu = XDMFFile(comm, tmp.format('u','RD',output_suffix))
      XDMFw = XDMFFile(comm, tmp.format('w','RD',output_suffix))
      XDMFpa = XDMFFile(comm, tmp.format('p','RD',output_suffix))
      XDMFpb = XDMFFile(comm, tmp.format('p','PML',output_suffix))
      XDMFu.parameters["rewrite_function_mesh"] = False
      XDMFw.parameters["rewrite_function_mesh"] = False
      XDMFpa.parameters["rewrite_function_mesh"] = False
      XDMFpb.parameters["rewrite_function_mesh"] = False
    if self.HDF5:
      tmp = "output/h5/{:s}_{:s}{:s}.h5"
      HDF5u = HDF5File(comm, tmp.format('u','RD',output_suffix), "w")
      HDF5w = HDF5File(comm, tmp.format('w','RD',output_suffix), "w")
      HDF5pa = HDF5File(comm, tmp.format('p','RD',output_suffix), "w")
      HDF5pb = HDF5File(comm, tmp.format('p','PML',output_suffix), "w")
    if self.PVD:
      tmp = "output/pvd/{:s}_{:s}{:s}.pvd"
      PVDu = File(tmp.format('u','RD',output_suffix),"compressed")
      PVDw = File(tmp.format('w','RD',output_suffix),"compressed")
      PVDpa = File(tmp.format('p','RD',output_suffix),"compressed")
      PVDpb = File(tmp.format('p','PML',output_suffix),"compressed")

    # Discrete solution
    Uh = BlockFunction(self.W)

    # Time stepping parameters
    t     = self.TimeStepping.t
    t_end = self.TimeStepping.t_end
    dt    = self.TimeStepping.dt
    gamma = self.TimeStepping.gamma
    beta  = self.TimeStepping.beta

    # Solving loop
    step = 0
    while t < t_end - 0.5*dt:

      # Update time
      t += float(dt)
      step += 1

      if rank == 0:
        print('\n\rtime: {:.4f} (Progress: {:.2f}%)'.format(t, 100*t/t_end),flush=True)

      # Update source term
      self.g.t = t

      # Assemble rhs and apply boundary condition
      b = block_assemble(self.rhs)
      self.bcs.apply(b)

      # Solve system
      self.solver.solve(Uh.block_vector(), b)
      Uh.block_vector().block_function().apply("to subfunctions")    
      (u, w, p, pi, S, _) = block_split(Uh)

      # Update previous time step
      u.rename('u','u')
      w.rename('w','w')
      p.rename('p','p')
      # update(Uh, self.U0, self.V0, self.A0, beta, gamma, dt)
      update(u, self.u0, self.v0, self.a0, beta, gamma, dt)
      update(w, self.wu0, self.wv0, self.wa0, beta, gamma, dt)
      update(p, self.pu0, self.pv0, self.pa0, beta, gamma, dt)
      update(pi, self.piu0, self.piv0, self.pia0, beta, gamma, dt)
      update(S, self.U0, self.V0, self.A0, beta, gamma, dt)

      # Text file for debugging
      if rank==0:
        with open('timesteps_hybrid.txt', 'a') as debug:
          debug.write('Saved timestep {:d} (time {:.4f} s)\n'.format(step,t))

      # Get current pressure on the PML layer
      # (_, _, _, _, _, pv0, _, _, _) = block_split(self.V0)
      self.piv0.rename('p','p')

      # Save solution to XDMF file format
      if self.XDMF:
        if (step % self.save_every_n == 0) or (step == 1):
          XDMFu.write(u, t)
          XDMFw.write(w, t)
          XDMFpa.write(p, t)
          XDMFpb.write(self.piv0, t)
      if self.HDF5:
        HDF5u.write(u, "u_RD", t)
        HDF5w.write(w, "w_RD", t)
        HDF5pa.write(p, "p_RD", t)
        HDF5pb.write(self.piv0, "p_PML", t)
      if self.PVD:
        if (step % self.save_every_n == 0) or (step == 1):
          PVDu << (u, t)
          PVDw << (w, t)
          PVDpa << (p, t)
          PVDpb << (self.piv0, t)

    # Close file
    if self.XDMF:
      XDMFu.close()
      XDMFw.close()
      XDMFpa.close()
      XDMFpb.close()
    if self.HDF5:
      HDF5u.close()
      HDF5w.close()
      HDF5pa.close()
      HDF5pb.close()
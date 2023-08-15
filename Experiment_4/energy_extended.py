# Copyright (C) 2019 Hernan Mella
import sys
from subprocess import run

import yaml

# setting path
sys.path.append('../')

import numpy as np
from dolfin import *
from mpi4py import MPI as pyMPI
from PML_base.ConstitutiveRelations import PoroelasticStress
from PML_base.Newmark import Newmark, update
from PML_base.PoroelasticParameters import Set4
from PML_base.TimeStepping import TimeStepping3


def mpi4py_comm(comm):
    '''Get mpi4py communicator'''
    try:
        return comm.tompi4py()
    except AttributeError:
        return comm
    
def peval(f, x):
    '''Parallel synced eval'''
    try:
        yloc = f(x)
    except RuntimeError:
        yloc = np.inf*np.ones(f.value_shape())

    comm = mpi4py_comm(f.function_space().mesh().mpi_comm())
    yglob = np.zeros_like(yloc)
    comm.Allreduce(yloc, yglob, op=pyMPI.MIN)

    return yglob

if __name__=="__main__":

  # MPI parameters
  comm = MPI.comm_world
  rank = MPI.rank(comm)

  # Set log level for parallel
  set_log_level(LogLevel.ERROR)
  if rank == 0: set_log_level(LogLevel.PROGRESS)

  # Import yaml parameters
  with open('PARAMETERS.yaml') as file:
      yparams = yaml.load(file, Loader=yaml.FullLoader)  

  # Import mesh
  mesh = Mesh()
  with XDMFFile("mesh/mesh_extended.xdmf") as infile:
      infile.read(mesh)

  # Markers for subdomains
  mvc = MeshValueCollection("size_t", mesh, mesh.geometry().dim())
  with XDMFFile("mesh/mf_extended.xdmf") as infile:
      infile.read(mvc, "markers")
  mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

  # Markers for layers
  mvc = MeshValueCollection("size_t", mesh, mesh.geometry().dim())
  with XDMFFile("mesh/mf_extended_1.xdmf") as infile:
      infile.read(mvc, "markers")
  mf_1 = cpp.mesh.MeshFunctionSizet(mesh, mvc)

  # Markers for Boundary conditions
  mvc = MeshValueCollection("size_t", mesh, mesh.geometry().dim()-1)
  with XDMFFile("mesh/ff_extended.xdmf") as infile:
      infile.read(mvc, "markers")
  ff = cpp.mesh.MeshFunctionSizet(mesh, mvc)

  # Create function spaces
  FE1 = FiniteElement("CG", mesh.ufl_cell(), 1)
  VE1 = VectorElement("CG", mesh.ufl_cell(), 2, dim=2)
  VE2 = VectorElement("CG", mesh.ufl_cell(), 2, dim=2)
  W = FunctionSpace(mesh, MixedElement([VE1,VE2,FE1]))

  # Time stepping parameters
  time_stepping = TimeStepping3(t_end=2.0)
  dt      = time_stepping.dt
  t       = time_stepping.t
  t_end   = time_stepping.t_end
  gamma   = time_stepping.gamma
  beta    = time_stepping.beta
  omega_p = 2*15*DOLFIN_PI

  # Function space for physical parameters
  P = FunctionSpace(mesh, 'DG', 0)

  # Get cell indices for each strata
  layer_1 = np.where(mf_1.array()[:] == 3)[0]
  layers = [layer_1]

  # Porous medium parameters
  params = Set4(V=P, layers=layers)
  rho_s = params.rho_s  # solid densitiy
  rho_f = params.rho_f  # fluid density
  K_s   = params.K_s    # bulk moduli of fluid
  K_f   = params.K_f    # bulk moduli of the solid
  K_b   = params.K_b    # bulk moduli of the dry porous frame
  mu_b  = params.mu_b   # tortuosity
  kappa = params.kappa  # permeability of the fluid
  phi   = params.phi    # porosity
  alpha = params.alpha  # Biot-Willis coefficient
  rho   = params.rho    # effective density
  rho_w = params.rho_w
  eta   = params.eta    # shear moduli of the dry porous frame
  T     = params.T      # dynamic viscosity of the fluid
  M     = params.M      # fluid-solid coupling modulus

  # Create measure for the source term
  dx  = Measure("dx", domain=mesh, subdomain_data=mf)

  # Set up initial values
  U0 = Function(W)
  (u0, wu0, pu0) = split(U0)

  V0 = Function(W)
  (v0, wv0, pv0) = split(V0)

  A0 = Function(W)
  (a0, wa0, pa0) = split(A0)

  # Newmark accelerations and velocities
  Nu = Newmark(u0=u0, v0=v0, a0=a0, beta=beta, gamma=gamma, dt=dt)
  Nw = Newmark(u0=wu0, v0=wv0, a0=wa0, beta=beta, gamma=gamma, dt=dt)
  Np = Newmark(u0=pu0, v0=pv0, a0=pa0, beta=beta, gamma=gamma, dt=dt)  

  # Stress tensor
  stress = PoroelasticStress(mu_b=mu_b,K_b=K_b,alpha=alpha,M=M)

  # Current values
  U = Function(W)
  (u, w, p) = U.split(True)

  # Input files
  path = ""
  hdf5ua_file = HDF5File(mesh.mpi_comm(), path+"output/h5/u_extended.h5", "r")
  hdf5wa_file = HDF5File(mesh.mpi_comm(), path+"output/h5/w_extended.h5", "r")
  hdf5pa_file = HDF5File(mesh.mpi_comm(), path+"output/h5/p_extended.h5", "r")

  # Output files (if required)
  save_velocities = False
  if save_velocities:
    xdmfuv_file = XDMFFile(mesh.mpi_comm(), "output/xdmf/uv_extended.xdmf")
    xdmfwv_file = XDMFFile(mesh.mpi_comm(), "output/xdmf/wv_extended.xdmf")
    xdmfuv_file.parameters["rewrite_function_mesh"] = False
    xdmfwv_file.parameters["rewrite_function_mesh"] = False

  # Energy
  if rank==0:
    energy_file = open('energy and seismograms/energy_extended.txt','w')
    energy_file.close()

  # Get seismograms locations
  seismograms = []
  for i in range(0,10):
    try:
      seismograms.append(yparams[3]['s{:1.0f}'.format(i+1)])
    except:
      break

  # Find closest node to seismograms
  Vu, Vu_to_W = W.sub(0).collapse(True)
  Vw, Vw_to_W = W.sub(1).collapse(True)
  Vp, Vp_to_W = W.sub(2).collapse(True)  
  dof_coords_u = Vu.tabulate_dof_coordinates()
  dof_coords_p = Vp.tabulate_dof_coordinates()
  seismograms_dofs = []
  print(seismograms)
  for s in range(len(seismograms)):
    # Distance from dof coordinate to each seismogram
    dist_u = np.linalg.norm(dof_coords_u - seismograms[s], axis=1)
    dist_p = np.linalg.norm(dof_coords_p - seismograms[s], axis=1)

    # find nearest DOF:
    dof_u = np.argmin(dist_u)
    dof_p = np.argmin(dist_p)
    seismograms_dofs.append([dof_u,dof_u+1,dof_u,dof_u+1,dof_p])
    # print('dof {}, x = {}'.format(dof_u, dof_coords_u[dof_u]))
    # print(seismograms_dofs)

  # Create files for seismograms
  if dist_u==0:
    for i in range(len(seismograms)):
      f = open('energy and seismograms/s{:1.0f}_extended.txt'.format(i+1),'w')
      f.close()

  # Solving loop
  step = 0
  while t < t_end - 0.5*dt:

    # Update time
    t += float(dt)
    step += 1

    if rank == 0:
      print('\n\rtime: {:.3f} (Progress: {:.2f}%)'.format(t, 100*t/t_end),)

    # Read checkpoints
    hdf5ua_file.read(u, "u/vector_{:d}".format(step-1))
    hdf5wa_file.read(w, "w/vector_{:d}".format(step-1))
    hdf5pa_file.read(p, "p/vector_{:d}".format(step-1))
    assign(U.sub(0), u)
    assign(U.sub(1), w)
    assign(U.sub(2), p)

    # Assemble energy
    E = assemble(0.5*rho*inner(Nu.dot(u), Nu.dot(u))*dx(1) \
        + 0.5*rho*inner(stress.tensor(u), grad(u))*dx(1) \
        + 0.5*rho_w*inner(Nw.dot(w), Nw.dot(w))*dx(1) \
        + 1/(2*M)*p*p*dx(1) \
        + rho_f*inner(Nu.dot(u), Nw.dot(w))*dx(1))

    # Export energy
    if rank==0:
      print(E)
      with open('energy and seismograms/energy_extended.txt','a') as f:
        f.write('{:.18e} {:.18e}\n'.format(t,E))

    # Get and export seismograms information
    for s in range(len(seismograms)):
      dofs = seismograms_dofs[s]
      t0, t1 = u.vector()[dofs[0]], u.vector()[dofs[1]]
      t2, t3 = w.vector()[dofs[2]], w.vector()[dofs[3]]
      t4 = p.vector()[dofs[4]]
      if rank==0:
        with open('energy and seismograms/s{:1.0f}_extended.txt'.format(s+1),'a') as f:
          f.write('{:.18e} {:.18e} {:.18e} {:.18e} {:.18e} {:.18e}\n'.format(t,t0,t1,t2,t3,t4))

    # Update previous time step
    update(U, U0, V0, A0, beta, gamma, dt)

    # Save velocities
    if save_velocities and ((step % 20 == 0) or (step == 1)):
      (v0_, wv0_, pv0_) = V0.split(True)
      xdmfuv_file.write(v0_, t)
      xdmfwv_file.write(wv0_, t)

  # Close files
  hdf5ua_file.close()
  hdf5wa_file.close()
  hdf5pa_file.close()
  if save_velocities:
    xdmfuv_file.close()
    xdmfwv_file.close()

  # Remove HDF5 Files
  if True:
    run(['rm','-rf','output/h5/u_extended.h5'])
    run(['rm','-rf','output/h5/w_extended.h5'])
    run(['rm','-rf','output/h5/p_extended.h5'])
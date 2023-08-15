# Copyright (C) 2019 Hernan Mella

import sys

import yaml

# setting path
sys.path.append('../')

import numpy as np
from dolfin import *
from PML_base.PoroelasticParameters import Set1
from PML_base.PoroelasticProblems import ParaxialProblem
from PML_base.Source import neumann_ricker_pulse
from PML_base.TimeStepping import TimeStepping1

if __name__=="__main__":

  # Set backend to PETSC
  parameters["linear_algebra_backend"] = "PETSc"

  # MPI parameters
  comm = MPI.comm_world
  rank = MPI.rank(comm)

  # Set log level for parallel
  set_log_level(LogLevel.ERROR)
  if rank == 0: set_log_level(LogLevel.PROGRESS)

  # Import geometric parameters
  with open('PARAMETERS.yaml') as file:
      params = yaml.load(file, Loader=yaml.FullLoader)

  # Import mesh
  mesh = Mesh(comm)
  with XDMFFile(comm, "mesh/mesh_paraxial.xdmf") as infile:
      infile.read(mesh)

  # Markers for subdomains
  mvc = MeshValueCollection("size_t", mesh, mesh.geometry().dim())
  with XDMFFile(comm, "mesh/mf_paraxial.xdmf") as infile:
      infile.read(mvc, "markers")
  mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

  # Markers for Boundary conditions
  mvc = MeshValueCollection("size_t", mesh, mesh.geometry().dim()-1)
  with XDMFFile(comm, "mesh/ff_paraxial.xdmf") as infile:
      infile.read(mvc, "markers")
  ff = cpp.mesh.MeshFunctionSizet(mesh, mvc)

  # Function space for physical parameters
  P = FunctionSpace(mesh, 'DG', 0)

  # Markers for layers
  mvc = MeshValueCollection("size_t", mesh, mesh.geometry().dim())
  with XDMFFile("mesh/mf_paraxial_1.xdmf") as infile:
      infile.read(mvc, "markers")
  mf1 = cpp.mesh.MeshFunctionSizet(mesh, mvc)

  # Get cell indices for each strata
  layer_1 = np.where(mf1.array()[:] == 3)[0]
  layers = [layer_1]

  # Porous medium parameters
  parameters = Set1(V=P, layers=layers)

  # Source term
  fr = params[3]['fr']
  A  = params[3]['source_amplitude']
  g = neumann_ricker_pulse(A=A, t=0.0, omega=2*fr*DOLFIN_PI, degree=2)

  # Define problem
  problem = ParaxialProblem(mesh=mesh, parameters=parameters, 
                          g=g, TimeStepping=TimeStepping1(),
                          domains=mf,
                          boundaries=ff,
                          solver_type='direct')

  # Solve the problem
  #problem.solve()
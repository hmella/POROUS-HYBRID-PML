# Copyright (C) 2019 Hernan Mella
import sys

import yaml

# setting path
sys.path.append('../')

import numpy as np
from dolfin import *
from PML_base.ConstitutiveRelations import PoroelasticStress
from PML_base.Newmark import Newmark, update
from PML_base.PoroelasticParameters import Set3
from PML_base.Source import neumann_ricker_pulse
from PML_base.TimeStepping import TimeStepping3
from PML_base.PoroelasticProblems import Problem

if __name__=="__main__":

  # Set backend to PETSC
  parameters["linear_algebra_backend"] = "PETSc"

  # MPI parameters
  comm = MPI.comm_world
  rank = MPI.rank(comm)

  # Set log level for parallel
  set_log_level(LogLevel.ERROR)
  if rank == 0: set_log_level(LogLevel.PROGRESS)

  # Import YAML parameters
  with open('PARAMETERS.yaml') as file:
      params = yaml.load(file, Loader=yaml.FullLoader)

  # Import mesh
  mesh = Mesh()
  with XDMFFile("mesh/mesh_extended.xdmf") as infile:
      infile.read(mesh)

  # Markers for subdomains
  mvc = MeshValueCollection("size_t", mesh, mesh.geometry().dim())
  with XDMFFile("mesh/mf_extended_1.xdmf") as infile:
      infile.read(mvc, "markers")
  mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

  # Markers for Boundary conditions
  mvc = MeshValueCollection("size_t", mesh, mesh.geometry().dim()-1)
  with XDMFFile("mesh/ff_extended.xdmf") as infile:
      infile.read(mvc, "markers")
  ff = cpp.mesh.MeshFunctionSizet(mesh, mvc)

  # Function space for physical parameters
  P = FunctionSpace(mesh, 'DG', 0)

  # Get cell indices for each strata
  layer_1 = np.where(mf.array()[:] == 3)[0]
  layer_2 = np.where(mf.array()[:] == 4)[0]
  layer_3 = np.where(mf.array()[:] == 5)[0]
  layer_4 = np.where(mf.array()[:] == 6)[0]
  layers = [layer_1,layer_2,layer_3,layer_4]

  # Porous medium parameters
  parameters = Set3(V=P, layers=layers)

  # Set up boundary condition
  uD = Constant(("0.0", "0.0"))
  wD = Constant(("0.0", "0.0"))

  # Source term
  fr = params[3]['fr']
  A  = params[3]['source_amplitude']
  g = neumann_ricker_pulse(A=A, t=0.0, omega=2*fr*DOLFIN_PI, degree=2)

  # Define problem
  problem = Problem(mesh=mesh, parameters=parameters, 
                          g=g, TimeStepping=TimeStepping3(t_end=2.0),
                          domains=mf,
                          boundaries=ff,
                          uD=[[uD, 1],],
                          wD=[[wD, 1],])

  # Solve the problem
  problem.solve()
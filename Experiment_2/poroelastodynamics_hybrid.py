# Copyright (C) 2019 Hernan Mella

import sys

import yaml

# setting path
sys.path.append('../')

import numpy as np
from fenics import *
from PML_base.Domains import map_facets
from PML_base.Geometry import Domain
from PML_base.PoroelasticParameters import Set3
from PML_base.PMLFun import PMLFunctions, PMLFunctionsNU
from PML_base.Source import neumann_ricker_pulse
from PML_base.TimeStepping import TimeStepping3
from PML_base.PoroelasticProblems import HybridPMLProblem

if __name__=="__main__":

  # Set backend to PETSC
  parameters["linear_algebra_backend"] = "PETSc"
  parameters["ghost_mode"] = "shared_facet"  # required by dS
  parameters["allow_extrapolation"] = True

  # MPI parameters
  comm = MPI.comm_world
  rank = MPI.rank(comm)

  # Set log level for parallel
  set_log_level(LogLevel.ERROR)
  if rank == 0: set_log_level(LogLevel.PROGRESS)

  # Import geometric parameters
  with open('PARAMETERS.yaml') as file:
      params = yaml.load(file, Loader=yaml.FullLoader)
  Lx        = params[0]['Lx']
  Ly        = params[0]['Ly']  
  LPML      = params[0]['LPML']
  source_center    = params[3]['source_center']
  b         = params[3]['b']

  # Geometry
  domain = Domain(Lx=Lx,Ly=Ly,LPML=LPML,source_center=source_center,b=b)

  # Import mesh
  mesh = Mesh(comm)
  with XDMFFile(comm,"mesh/mesh_hybrid.xdmf") as infile:
      infile.read(mesh)

  # Markers for subdomains
  mvc = MeshValueCollection("size_t", mesh, mesh.geometry().dim())
  with XDMFFile(comm,"mesh/mf_hybrid.xdmf") as infile:
      infile.read(mvc, "markers")
  mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

  # Markers for Boundary conditions
  mvc = MeshValueCollection("size_t", mesh, mesh.geometry().dim()-1)
  with XDMFFile(comm,"mesh/ff_hybrid.xdmf") as infile:
      infile.read(mvc, "markers")
  ff = cpp.mesh.MeshFunctionSizet(mesh, mvc)

  # Markers for layers
  mvc = MeshValueCollection("size_t", mesh, mesh.geometry().dim())
  with XDMFFile(comm,"mesh/mf_hybrid_1.xdmf") as infile:
      infile.read(mvc, "markers")
  mf1 = cpp.mesh.MeshFunctionSizet(mesh, mvc)

  # Get cell indices for each strata
  layer_1 = np.where(mf1.array()[:] == 3)[0]
  layer_2 = np.where(mf1.array()[:] == 4)[0]
  layer_3 = np.where(mf1.array()[:] == 5)[0]
  layer_4 = np.where(mf1.array()[:] == 6)[0]
  layers = [layer_1, layer_2, layer_3, layer_4]

  # Function space for physical parameters
  P = FunctionSpace(mesh, 'DG', 0)

  # Porous medium parameters
  parameters = Set3(V=P, layers=layers)

  # Set up boundary condition
  uD_PML_1 = Constant(("0.0", "0.0"))
  wD_PML_1 = Constant(("0.0", "0.0"))
  pD_PML_1 = Constant("0.0")

  # Solve the proble for different PML function
  for suffix in ['','_multiaxial']:
      # Attenuation and stretching functions and tensors
      if suffix == '':
        PML = PMLFunctions(domain,parameters,form=1,R=1e-4,m=3,alpha_max=5.0)
      elif suffix == '_multiaxial':
        PML = PMLFunctions(domain,parameters,form=1,R=1e-4,m=3,alpha_max=0.0,PMLType='MPML',Pxy=0.01,Pyx=0.01)
      elif suffix == '_optimal':
        PML = PMLFunctions(domain,parameters,form=1,m=1,alpha_max=0.0,PMLType='OPML')

      # Source term
      fr = params[3]['fr']
      A  = params[3]['source_amplitude']
      g = neumann_ricker_pulse(A=A, t=0.0, omega=2*fr*DOLFIN_PI, degree=2)

      # Define problem
      problem = HybridPMLProblem(mesh=mesh,
                              domains=mf,
                              boundaries=ff,
                              parameters=parameters, 
                              g=g, TimeStepping=TimeStepping3(),
                              uD=[[uD_PML_1, 1],],
                              wD=[[wD_PML_1, 1],],
                              pD=[[pD_PML_1, 1],],
                              PMLFunctions=PML)

      # Solve the problem
      problem.solve(output_suffix=suffix)
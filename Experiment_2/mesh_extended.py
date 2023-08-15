# Copyright (C) 2019 Hernan Mella
import sys

import yaml

# setting path
sys.path.append('../')


import meshio
import numpy as np
from PML_base import Geometry

if __name__ == '__main__':

  # Import geometric parameters
  with open('PARAMETERS.yaml') as file:
      params = yaml.load(file, Loader=yaml.FullLoader)
  Lx        = params[2]['Lx']
  Ly        = params[2]['Ly']
  LPML      = params[0]['LPML']
  R_ext     = params[2]['R_ext']
  R_int     = params[2]['R_int']
  meshtype  = params[2]['meshtype']
  source_center    = params[3]['source_center']
  b         = params[3]['b']
  cl_base   = params[3]['cl_base']
  cl_coarse = params[3]['cl_coarse']
  cl_source = params[3]['cl_source']

  # Create domain and generate mesh
  domain = Geometry.Domain(Lx=Lx, Ly=Ly, LPML=LPML, R_ext=R_ext, R_int=R_int,
                source_center=source_center, b=b, cl_source=cl_source,
                cl_base=cl_base, cl_coarse=cl_coarse, meshtype=meshtype,
                flag=0)
  domain.generate_mesh()

  # Read msh mesh
  path = "mesh/"

  # Read mesh
  geometry = meshio.read(path+"mesh.msh")

  # Indices to obtain lines and triangles markers
  idx_lines     = [i for (i,c) in enumerate(geometry.cells) if c.type=='line']
  idx_triangles = [i for (i,c) in enumerate(geometry.cells) if c.type=='triangle']

  # Mesh points
  points = geometry.points[:,0:2]

  # Mesh cells and phyiscal data
  cells = np.vstack([c.data for c in geometry.cells if c.type=='triangle'])
  cell_data = np.hstack([geometry.cell_data['gmsh:physical'][i] for i in idx_triangles])

  # Separate cell_data into RD and PML region markers and stratas
  domains = [[3,8,13,18],[1,2,4,5,6,7,9,10,11,12,14,15,16,17]]
  cell_data_domains = np.array(cell_data)
  for idx, domain in enumerate(domains):
    for domain_i in domain:
      cell_data_domains[np.where(cell_data == domain_i)] = idx + 1

  stratas = [[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15],[16,17,18]]
  cell_data_stratas = np.array(cell_data)
  for idx, domain in enumerate(stratas):
    for domain_i in domain:
      cell_data_stratas[np.where(cell_data == domain_i)] = 3 + idx

  # Mesh facets and phyiscal data
  lines = np.vstack([c.data for c in geometry.cells if c.type=='line'])
  line_data = np.hstack([geometry.cell_data['gmsh:physical'][i] for i in idx_lines])

  # Export in xdmf
  meshio.write(path+"mesh_extended.xdmf", meshio.Mesh(points=points,
              cells={"triangle": cells}))

  meshio.write(path+"mf_extended.xdmf", meshio.Mesh(points=points,
              cells={"triangle": cells},
              cell_data={"markers": [cell_data_domains]}))

  meshio.write(path+"mf_extended_1.xdmf", meshio.Mesh(points=points,
              cells={"triangle": cells},
              cell_data={"markers": [cell_data_stratas]}))

  meshio.write(path+"ff_extended.xdmf", meshio.Mesh(points=points,
              cells={"line": lines},
              cell_data={"markers": [line_data]}))

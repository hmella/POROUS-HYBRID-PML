# Copyright (C) 2019 Hernan Mella

import os
from subprocess import run


class Domain:
  def __init__(self, Lx=10, Ly=10, Lz=-1, LPML=1,
               R_ext=100, R_int=50,
               b=1, source_center=[0,0,0],
               cl_base=5, cl_source=1, cl_coarse=5,
               cl_fine=1.0, cl_PML=5, cl_str_1=5.0,
               cl_str_2=5.0, cl_str_3=5.0, cl_str_4=5.0,
               cl_str_5=5.0, meshtype='GEO_1_HYBRID',
               flag=0):
    self.Lx = Lx
    self.Ly = Ly
    self.Lz = Lz
    self.LPML = LPML
    self.R_ext = R_ext
    self.R_int = R_int
    self.source_center = source_center
    self.b = b
    self.cl_base = cl_base
    self.cl_coarse = cl_coarse    
    self.cl_fine = cl_fine
    self.cl_PML = cl_PML
    self.cl_source = cl_source    
    self.cl_str_1 = cl_str_1    
    self.cl_str_2 = cl_str_2    
    self.cl_str_3 = cl_str_3    
    self.cl_str_4 = cl_str_4 
    self.cl_str_5 = cl_str_5
    self.meshtype = meshtype
    self.flag = flag

  def generate_mesh(self):

    # Current working directory
    cwd = os.getcwd()

    # Create mesh folder
    folder = os.path.dirname('{:s}/mesh/'.format(cwd))
    if not os.path.exists(folder):
      os.makedirs(folder)

    # Split the path in head and tail pair
    cwd_ht = os.path.split(cwd)

    # write gmsh file
    file = open('{:s}/PML_base/gmsh_geometries/dimensions'.format(cwd_ht[0]), 'w')
    dimensions = 'Lx = {:f};            // width of the interior domain'.format(self.Lx)
    dimensions += '\nLy = {:f};         // height of the interior domain'.format(self.Ly)
    dimensions += '\nLz = {:f};         // height of the interior domain'.format(self.Lz)
    dimensions += '\nLPML = {:f};       // width of the PML layer'.format(self.LPML)
    dimensions += '\nR_ext = {:f};         // external radius'.format(self.R_ext)
    dimensions += '\nR_int = {:f};         // internal radius'.format(self.R_int)
    dimensions += '\ncl_base = {:f};'.format(self.cl_base)
    dimensions += '\ncl_coarse = {:f};'.format(self.cl_coarse)
    dimensions += '\ncl_fine = {:f};'.format(self.cl_fine)
    dimensions += '\ncl_PML = {:f};'.format(self.cl_PML)
    dimensions += '\ncl_source = {:f};'.format(self.cl_source)
    dimensions += '\ncl_str_1 = {:f};'.format(self.cl_str_1)
    dimensions += '\ncl_str_2 = {:f};'.format(self.cl_str_2)
    dimensions += '\ncl_str_3 = {:f};'.format(self.cl_str_3)
    dimensions += '\ncl_str_4 = {:f};'.format(self.cl_str_4)
    dimensions += '\ncl_str_5 = {:f};'.format(self.cl_str_5)
    dimensions += '\nb = {:f};          // width of the Neumann boundary condition'.format(self.b)
    dimensions += '\nxc = {:f};         // center of the circle (x-coord)'.format(self.source_center[0])
    dimensions += '\nyc = {:f};         // center of the circle (y-coord)'.format(self.source_center[1])
    dimensions += '\nflag = {:d};         // height of the interior domain'.format(self.flag)
    print(dimensions)
    file.write(dimensions)
    file.close()

    # Generate PML meshes depending on the input parameters
    path_geo_file = '{:s}/PML_base/gmsh_geometries/{:s}.geo'.format(cwd_ht[0],self.meshtype)
    path_msh_file = 'mesh/mesh.msh'
    if self.Lz > 0:
      run(['gmsh',path_geo_file,'-3','-o',path_msh_file])
    else:
      run(['gmsh',path_geo_file,'-2','-o',path_msh_file])

    return True

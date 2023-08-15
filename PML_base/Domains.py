# Copyright (C) 2019 Hernan Mella

from dolfin import DOLFIN_EPS, SubDomain, SubsetIterator, MeshFunction, near


# Map facets between parent and child meshes to apply boundary conditions
def map_facets(parent_mesh,child_mesh,markers):

  # Mapped facets
  dim = child_mesh.geometry().dim()-1
  child_markers = MeshFunction("size_t", child_mesh, dim, 0)

  # Boundary facets in child mesh
  child_tmp = MeshFunction("size_t", child_mesh, dim, 0)
  Boundary().mark(child_tmp, 1)

  # Boundary facets in parent mesh
  parent_tmp = MeshFunction("size_t", parent_mesh, dim, 0)
  Boundary().mark(parent_tmp, 1)

  # Maximun label in parent facet functions
  labels = [*set(markers.array())]

  for label in labels[:-1]:
    for f_parent in SubsetIterator(markers, label):
      for f_child in SubsetIterator(child_tmp, 1):
        if f_parent.midpoint().distance(f_child.midpoint()) <= DOLFIN_EPS:
          child_markers[f_child.index()] = markers[f_parent.index()]

  return child_markers



class Dirichlet(SubDomain):
  def __init__(self,Lx,Ly,Lpml,**kwargs):
    self.Lx = Lx
    self.Ly = Ly
    self.Lpml = Lpml
    super().__init__(**kwargs)

  def inside(self, x, on_boundary):
    return x[0] <= -0.5*self.Lx - self.Lpml + DOLFIN_EPS \
        or x[0] >= 0.5*self.Lx + self.Lpml - DOLFIN_EPS \
        or x[1] <= -self.Ly-self.Lpml + DOLFIN_EPS \
        or x[1] >= self.Lpml - DOLFIN_EPS and on_boundary

class Neumann(SubDomain):
  def __init__(self,xc,width,**kwargs):
    self.xc = xc
    self.width = width
    super().__init__(**kwargs)

  def inside(self, x, on_boundary):
    halfwidth = 0.5*self.width + 1e-6*self.width
    return near(x[1], 0) and (-halfwidth + self.xc <= x[0] <= halfwidth + self.xc)

class NotNeumann(SubDomain):
  def __init__(self,xc,width,**kwargs):
    self.xc = xc
    self.width = width
    super().__init__(**kwargs)

  def inside(self, x, on_boundary):
    halfwidth = 0.5*self.width - 1e-6*self.width
    return near(x[1], 0) and ((x[0] < -halfwidth + self.xc) or (x[0] > halfwidth + self.xc))

class Interface(SubDomain):
  def __init__(self,Lx,Ly,**kwargs):
    self.Lx = Lx
    self.Ly = Ly
    super().__init__(**kwargs)

  def inside(self, x, on_boundary):
    return (near(x[0],-0.5*self.Lx) and x[1] >= -self.Ly - DOLFIN_EPS) \
        or (near(x[0],0.5*self.Lx) and x[1] >= -self.Ly - DOLFIN_EPS) \
        or (near(x[1],-self.Ly) and x[0] >= -0.5*self.Lx - DOLFIN_EPS and x[0] <= 0.5*self.Lx + DOLFIN_EPS)

# Source domain
class Circle(SubDomain):
  def __init__(self,c,rd,**kwargs):
    self.c  = c
    self.rd = rd
    super().__init__(**kwargs)

  def inside(self, x, on_boundary):
    return pow(x[0] - self.c[0], 2) + pow(x[1] - self.c[1], 2) <= pow(self.rd, 2) + 1.0e-6

# Regular domain
class RegularDomain(SubDomain):
  def __init__(self,Lx,Ly,**kwargs):
    self.Lx = Lx
    self.Ly = Ly
    super().__init__(**kwargs)

  def inside(self, x, on_boundary):
    return abs(x[0]) <= 0.5*self.Lx and x[1] >= -self.Ly

# Boundaries
class Boundary(SubDomain):
  def __init__(self,**kwargs):
    super().__init__(**kwargs)

  def inside(self, x, on_boundary):
    return on_boundary

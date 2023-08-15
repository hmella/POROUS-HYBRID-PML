# Copyright (C) 2019 Hernan Mella
from dolfin import Function
import ufl
from multiphenics import BlockFunction


# Update previous time step using Newmark scheme
def update(u, u0, v0, a0, beta, gamma, dt):
  """Update fields at the end of each time step."""

  # Get vectors (references)
  if isinstance(u, Function):
    u_vec, u0_vec  = u.vector(), u0.vector()
    v0_vec, a0_vec = v0.vector(), a0.vector()
  elif isinstance(u, BlockFunction):
    u_vec, u0_vec  = u.block_vector(), u0.block_vector()
    v0_vec, a0_vec = v0.block_vector(), a0.block_vector()

  # Update acceleration and velocity
  # a = 1/(beta*dt*dt)*(u - u0 - v0*dt) - (1-2*beta)/(2*beta)*a0
  a_vec = (1.0/(beta*dt*dt))*(u_vec - u0_vec - v0_vec*dt) - (1.0-2.0*beta)/(2*beta)*a0_vec

  # v = v0 + (1-gamma)*dt*a0 + gamma*dt*a
  v_vec = v0_vec + (1.0 - gamma)*dt*a0_vec + gamma*dt*a_vec

  # Update (u0 <- u)
  if isinstance(u, Function):
    v0.vector()[:], a0.vector()[:] = v_vec, a_vec
    u0.vector()[:] = u.vector()
  elif isinstance(u, BlockFunction):
    v0.block_vector()[:], a0.block_vector()[:] = v_vec, a_vec
    u0.block_vector()[:] = u.block_vector()

# Newmark time integration
class Newmark:
  def __init__(self, u0=[], v0=[], a0=[], beta=0.25, gamma=0.5, dt=0.001):
    self.u0 = u0
    self.v0 = v0
    self.a0 = a0
    self.beta = beta
    self.gamma = gamma
    self.dt = dt

  # Acceleration
  def ddot(self, u):
    return (u - self.u0 - self.dt*self.v0)/(self.beta*self.dt*self.dt) \
          - self.a0*(1.0 - 2.0*self.beta)/(2.0*self.beta)
  
  # Velocity
  def dot(self, u):
    return self.gamma*(u - self.u0)/(self.beta*self.dt) \
          + self.v0*(1 - self.gamma/self.beta) \
          + self.a0*self.dt*(1.0 - self.gamma/(2*self.beta))


# Newmark time integration for the hybrid PML
class NewmarkHybrid:
  def __init__(self, u0=[], v0=[], a0=[], beta=0.25, gamma=0.5, dt=0.001):
    self.u0 = u0
    self.v0 = v0
    self.a0 = a0
    self.beta = beta
    self.gamma = gamma
    self.dt = dt

  # Acceleration
  def ddotl(self, u):
    left  = u/(self.beta*self.dt*self.dt)
    return left

  def ddotr(self, u):
    right = (self.u0 + self.dt*self.v0)/(self.beta*self.dt*self.dt) \
          + self.a0*(1.0 - 2.0*self.beta)/(2.0*self.beta)
    return right

  # Velocity
  def dotl(self, u):
    left  = self.gamma*u/(self.beta*self.dt)
    return left

  def dotr(self, u):
    right = self.gamma*self.u0/(self.beta*self.dt) \
          - self.v0*(1 - self.gamma/self.beta) \
          - self.a0*self.dt*(1.0 - self.gamma/(2*self.beta))
    return right
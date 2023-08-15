# Copyright (C) 2019 Hernan Mella
import numpy as np
from dolfin import DOLFIN_PI, UserExpression
from ufl import exp


# Neumann explosive Ricker source
class neumann_ricker_pulse(UserExpression):
  def __init__(self, t, omega, A=1.0e+3, gdim=2, **kwargs):
    self.A     = A
    self.t     = t
    self.omega = omega
    self.cutoff = 6.0*np.sqrt(6.0)/omega
    self.gdim = gdim
    super().__init__(**kwargs)

  def eval(self, values, x):
    # Ricker pulse parameters
    u  = self.omega*self.t - 3*pow(6,0.5)
    Tp = ((0.25*pow(u, 2) - 0.5)*exp(-0.25*pow(u, 2)) \
         -13*exp(-13.5))/(0.5 + 13*exp(-13.5))

    values[0] = 0.0
    values[1] = 0.0
    if self.t <= self.cutoff:
      values[self.gdim-1] = self.A*Tp
    else:
      values[self.gdim-1] = 0.0

  def value_shape(self):
    return (self.gdim,)
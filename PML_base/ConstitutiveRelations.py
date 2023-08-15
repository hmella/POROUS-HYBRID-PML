from dolfin import Identity, grad, sym, tr

# Strain tensor
def eps(u):
  return (grad(u) + grad(u).T)/2

# Stress tensor
class PoroelasticStress:
  def __init__(self,mu_b,K_b,alpha,M):
    self.mu_b = mu_b
    self.lmbda_b = K_b - 2/3*mu_b + alpha**2*M

  def tensor(self,u):
    gdim = self.mu_b.function_space().mesh().geometry().dim()
    return 2.0*self.mu_b*sym(grad(u)) + self.lmbda_b*tr(sym(grad(u)))*Identity(gdim)

# Small deformations tensor
class PoroelasticCompliance:
  def __init__(self,mu_b,K_b,alpha,M):
    self.mu_b = mu_b
    self.lmbda_b = K_b - 2/3*mu_b + alpha**2*M

  def tensor(self,sigma):
    gdim = self.mu_b.function_space().mesh().geometry().dim()
    return sigma/(2.0*self.mu_b) \
      - self.lmbda_b/(4.0*self.mu_b*(self.lmbda_b + self.mu_b))*tr(sigma)*Identity(gdim)

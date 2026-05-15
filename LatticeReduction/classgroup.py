from sage.all import *


class ClassGroup:
  def __init__(self, d):
    self.d = d
  
  def inner_product(self, v, u):
    return v[0] * u[0] + abs(self.d) * v[1] * u[1]
  
  def norm(self, u):
    return self.inner_product(u, u)  #i.e. u1^2 + |d| u2^2 NOTE: squared norm







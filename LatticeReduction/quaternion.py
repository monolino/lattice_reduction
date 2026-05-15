from sage.all import *
'''
implicit basis (1,i,j,k) where
i^2 = -1
j^2 = -p
k = ij = -ji
'''
class QuaternionAlgebra:
  def __init__(self, p):
    self.p = p
  
  def nrd(self, x):
    return x[0]**2 + x[1]**2 + x[2]**2* self.p + x[3]**2* self.p
  
  def nrd_ideal(self, I):
    return gcd([self.nrd(x) for x in I])

  def print_el(self, x):
    terms = []
    if x[0] != 0:
      terms.append(str(x[0]))
    if x[1] != 0:
      terms.append(f"{x[1]}i")
    if x[2] != 0:
      terms.append(f"{x[2]}j")
    if x[3] != 0:
      terms.append(f"{x[3]}k")
    if not terms:
      print("0")
    else:
      s = " + ".join(terms).replace("+ -", "- ")
      print(s)

  
  

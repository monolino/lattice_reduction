from sage.all import *

'''
Class group Q(sqrt{d}) for d squarefree
element v of the class group is: v = v1 + i sqrt{d} v2
can be represented as a vector (v1,v2) in Q^2

||u|| = u1^2 + d u2^2       (-> field norm)
<v,u> = v1 u1 + d v2 u2 = u^T G v, with
G = [
  1, 0
  0, d]
'''
def nearest_integer(q):
  return (q + QQ(1)/2).floor()

def inner_product(v, u, d):
  return v[0] * u[0] + d * v[1] * u[1]

def norm(u, d):
  return inner_product(u, u, d)  #i.e. u1^2 + d u2^2

def lattice_reduction_2dim(v, u, d):
  #Assumption that ||v|| <= ||u||
  if norm(v,d) > norm(u,d):
    v, u = u, v 

  #iteration count L
  L = 0

  while norm(v,d) < norm(u,d):
    q = inner_product(u,v,d) / norm(v,d)
    q = nearest_integer(q)
    r = (u[0] - q * v[0], u[1] - q * v[1])  #u - q * v  faster than using vector(QQ) operations
    u, v = v, r
    L += 1
  test_minimal_basis(u,v,d)
  return (u,v), L

def test_minimal_basis(v, u, d):
  assert norm(v,d) <= norm(u,d) and abs(inner_product(v,u,d)) <= norm(v,d)/2, f"The basis ({v}, {u}) is not reduced"

if __name__ == "__main__":
  d = QQ(5)
  u = (QQ(1),QQ(0))
  v = (QQ(3),QQ(1))
  (vm,um), L = lattice_reduction_2dim(u,v,d)
  print(f"Reduced basis: {vm}, {um}, L={L}")



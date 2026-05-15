from sage.all import *
'''
Needed an Class element with methods
inner_product(v,u) and
norm(u)
'''
def nearest_integer(q):
  return (q + QQ(1)/2).floor()

#NOTE: Important!! the norm is the squared norm here.
def lattice_reduction_2dim(Obj, v, u):
  assert all(isinstance(x, (int, Integer, Rational)) for x in v + u), "v and u should lie in QQ^2"

  #Assumption that ||v|| <= ||u||
  if Obj.norm(v) > Obj.norm(u):
    print('Case B: swap')
    v, u = u, v 

  #iteration count L
  L = 0

  while Obj.norm(v) < Obj.norm(u):
    q = Obj.inner_product(u,v) / Obj.norm(v)
    q = nearest_integer(q)
    r = (u[0] - q * v[0], u[1] - q * v[1])   #u - q * v  faster than using vector(QQ) operations
    u, v = v, r
    L += 1
  test_minimal_basis(u,v)
  return (u,v), L

def test_minimal_basis(v, u):
  assert all(isinstance(x, (int, Integer, Rational)) for x in v + u), "v and u should lie in QQ^2"
  assert (Obj.norm(v) <= Obj.norm(u) and abs(Obj.inner_product(v,u)) <= Obj.norm(v)/2) or (Obj.norm(u) <= Obj.norm(v) and abs(Obj.inner_product(u,v)) <= Obj.norm(u)/2 ), f"The basis ({v}, {u}) is not reduced"


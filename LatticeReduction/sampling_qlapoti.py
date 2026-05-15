from sage.all import *
from L_rdc_algorithm import *
from quaternion import QuaternionAlgebra
from euclideanspace import EuclideanSpace
import random

p = 5 #TODO what is p in our quaternion algebra B(-1,p)?
B = QuaternionAlgebra(p)

def small_generator(n, QA=B):
  '''I is left ideal, n is the integer n=nrd(I) and b is the basis of I
  returns alpha
  '''
  while True:
    coeffs = [random.randint(1, 10**5) for _ in range(4)]
    if gcd(2*coeffs[0], n) == 1 and gcd(QA.nrd(coeffs), n**2) == 1:
      return coeffs

def lattice_vectors(n, QA=B):
  coeffs = small_generator(n, QA)
  u = (QQ(n), 0) #n is gcd so integer
  v = (QQ(-2*coeffs[1]), QQ(2* coeffs[0])) # coeffs are all integer
  return u,v

def histogram(n, QA=B, trials=1000):
  #NOTE cannot take v as the key since the same lattice has different basis v.
  return

if __name__ == "__main__":
  n = 1001
  u,v = lattice_vectors(n, B)
  E = EuclideanSpace(2)
  print(f"u: {u}, v: {v}")
  try:
    test_minimal_basis(E, v, u)
    print("The basis is minimal")
  except AssertionError:
    print("The basis is not minimal")

  

  


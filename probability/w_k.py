from sage.all import *
from DominantEigenvalue.DominantEigenvalue import *
from SubdominantEigenvalue.SubdominantEigenvalue import *
from ConstantC.constant import *

w_1_val = 0.28986813369645287294
w_2_val = 0.04848080144946363270
w_3_val = 0.01027816477906659643

precision = 53

def w_1():
  return pi.n(precision)**2/3-3

def w_2(N=1000):
  val = -5 + 2*pi.n(precision)**2/3 - 2*zeta(3)
  for n in range(1, N+1):
    val += 2*((-1)**n *(n+1)*zeta(n+4)*(zeta(n+2)-1))
  return val

def w_k(k,m=64,M=1000):
  lambda_dom = p(m)
  mu = subdominant_eigenvalue_upperbound(lam = lambda_dom, M=M)
  c = constant_c(m)
  return c * lambda_dom ** k * (1 + (abs(mu)/abs(lambda_dom))**k)

if __name__ == "__main__":
  #print(f"w_1: {w_1().n(precision)}")
  #print(f"w_2: {w_2().n(precision)}")
  k = 2
  M=1000
  w_k_val = w_k(k)
  print(f"w_{k}: {w_k_val.n(precision)}")
  print(f"paper w_{k}: {w_2.n(precision)}")
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
  print(f"lambda_dom: {lambda_dom},\nmu: {mu},\nc: {c}")
  print(f"symbolic w_{k}: {c*lambda_dom**k} + {c*lambda_dom**k * (abs(mu)/abs(lambda_dom))**k} * CONST")
  return c * lambda_dom ** k * (1 + (abs(mu)/abs(lambda_dom))**k)

def find_k_for_eps(epsilon, m=64, M=1000):
  lambda_dom = 0.199456936618669145744216645005185152214178590575950745618862326297299505461605521321064234246552894856623201238356834575285212183630482315506366412927878972060137032746047371063233 #p(m)
  mu = 0.0817848113104025438918947757350936382209595595733334792637958012116893651486785522213993792 #subdominant_eigenvalue_upperbound(lam = lambda_dom, M=M)
  c = 0.45696720642905160712078870828525952183190644320652401921920 #constant_c(m)
  k = 1
  while True:
    w_k_val = c*lambda_dom**k
    if w_k_val < epsilon:
      print(f"w_{k} = {w_k_val} + {w_k_val * (abs(mu)/abs(lambda_dom))**k} * CONST")
      return k
    k += 1

if __name__ == "__main__":
  #print(f"w_1: {w_1().n(precision)}")
  #print(f"w_2: {w_2().n(precision)}")
  k = 3
  M=1000
  w_k_val = w_k(k, m=8, M=M)
  print(f"w_{k}: {w_k_val.n(precision)}")
  #print(f"paper w_{k}: {w_2_val}")

  epsilon = 0.0001
  #k = find_k_for(epsilon, m=64, M=1000)
  #print(f"Smallest k such that w_k < {epsilon}: {k}")
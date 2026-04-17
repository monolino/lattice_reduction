from sage.all import *
from pathlib import Path

precision = 600
CF = ComplexField(precision)

#their values
tr = 0.1444623962461608158824990905254832038136
trG_sq = 0.046471825642727939835279753170
lam = 0.199458818343767
mu = 0.08178022638533167

def tau_4(m):
  return CF(1)/2 * (m + (m**2 + 4)**(1/2))

def tau_4_two_m(m1,m2):
  #|h|=2 i.e. G^2
  return (CF(1)/2 * ((m1 * m2 + 2) + sqrt((m1 * m2 + 2)**2 - 4))).n(precision) 

def trace_4():
  val = CF(0)
  for m in range(1, M+1):
    val += tau_4(m)**(-4) / (1 + tau_4(m)**(-2))
  print(f"Trace of G: {val}")
  print(f"Trace of G squared: {val**2}")
  return val

def trace_4_G_sq():
  #Tr(G^2)
  val = CF(0)
  print(M)
  for m1 in range(1, M+1):
    for m2 in range(1, M+1):
      val += tau_4_two_m(m1,m2)**(-4) / (1 - tau_4_two_m(m1,m2)**(-2)) #what is correct
  print(f"Trace of G squared: {val}")
  return val

def get_lambda():
  base = Path(__file__).resolve().parent
  path = base.parent / "DominantEigenvalue" / f"lambda_{precision}"
  return CF(path.read_text())


def subdominant_eigenvalue_upperbound():
  return real_part((trace_4_G_sq() - get_lambda()**2)**(CF(1)/2))

if __name__ == "__main__":
  M = 1000
  mu = subdominant_eigenvalue_upperbound()
  print(f"mu = {mu}")
  print(f"Precision: {mu.precision()} bits")

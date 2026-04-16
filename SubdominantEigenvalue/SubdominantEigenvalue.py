from sage.all import *
from pathlib import Path

CF = ComplexField(500)

def tau_4(m):
  return CF(1)/2 * (m + (m**2 + 4)**(1/2))

def trace_4():
  val = CF(0)
  for m in range(1, M+1):
    val += tau_4(m)**(-4) / (1 + tau_4(m)**(-2))
  print(f"Trace of G: {val}")
  print(f"Trace of G squared: {val**2}")
  return val

def get_lambda():
  base = Path(__file__).resolve().parent
  path = base.parent / "DominantEigenvalue" / "lambda"
  return CF(path.read_text())


def subdominant_eigenvalue_upperbound(m):
  return real_part((trace_4()**2 - get_lambda()**2)**(1/2))

if __name__ == "__main__":
  M = 1000
  print("Subdominant eigenvalue upper bound: ", subdominant_eigenvalue_upperbound(M))
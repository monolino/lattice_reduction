from sage import *
from DominantEigenvalue.DominantEigenvalue import *

a= QQ(1)/2 #1/2 as element in Q
x = var('x')
precision = 200
CF = ComplexField(precision)


def u(t):
  return 1/(1+t)**2

#u(t) represented via the same basis (x-a)^i
def coeff_vector_u(m):
    # Taylor expansion of u(t) around t = a   
  t = var('t')
  
  expr = 1/(1+t)**2
  coeffs = []
  for k in range(m):
    ck = diff(expr, t, k)(t=a) / factorial(k)
    coeffs.append(CF(ck))
  return vector(CF, coeffs)




def constant_c(m):
  A = matrix_T_m(m)
  data_left = A.eigenvectors_left()
  data_right = A.eigenvectors_right()
  lambda_dom_r, ev_right, mult_r = max(data_right, key=lambda data: abs(data[0])) #right eigenvector corresponding to the dominant eigenvalue
  lambda_dom_l, ev_left, mult_l = max(data_left, key=lambda data: abs(data[0])) #left eigenvector corresponding to the dominant eigenvalue
  
  ev_right = ev_right[0]
  ev_left = ev_left[0]

  #scale right eigenvector such that f(0) = 1 TODO: do this normalization also for dominant eigenvector.
  ev_right = ev_right / ev_right[0]

  #normalize the eigenvectors
  norm = ev_left.dot_product(ev_right)
  ev_left = ev_left / norm

  #f*[u] = ev_levt * [u]
  return real_part(ev_left.dot_product(coeff_vector_u(m)))

if __name__ == "__main__":
  m = 16
  c = constant_c(m)
  print(f"Constant c for m={m}: {c}")
  

from sage.all import *
from DominantEigenvalue.DominantEigenvalue import matrix_T_m

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

def dual_f_star(x, m):
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

  
  if abs(imag_part(ev_left.dot_product(x))) > 1e-10:
    print("Warning: significant imaginary part:", val)


  return real_part(ev_left.dot_product(x))

def eigenvector_f(x, m):
  A = matrix_T_m(m)
  data_left = A.eigenvectors_left()
  data_right = A.eigenvectors_right()
  lambda_dom_r, ev_right, mult_r = max(data_right, key=lambda data: abs(data[0])) #right eigenvector corresponding to the dominant eigenvalue
  lambda_dom_l, ev_left, mult_l = max(data_left, key=lambda data: abs(data[0])) #left eigenvector corresponding to the dominant eigenvalue
  
  ev_right = ev_right[0]
  ev_left = ev_left[0]

  #scale right eigenvector such that f(0) = 1 TODO: do this normalization also for dominant eigenvector.
  ev_right = ev_right / ev_right[0]

  if abs(imag_part(ev_right.dot_product(x))) > 1e-10:
    print("Warning: significant imaginary part:", val)

  return real_part(ev_right.dot_product(x))

def Constant_C_4(k, m, pdf, grid_res=100):
  def integral(m, pdf, grid_res):
    #grid of points since only in 2D
    x = np.linspace(0, 1, grid_res)
    y = np.linspace(-0.5, 0.5, grid_res)
    X, Y = np.meshgrid(x, y)
    #circle mask
    mask = (X - 0.5)**2 + Y**2 <= 0.25

    #space between points
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    integrand_vals = np.zeros_like(X) #zero array for integrand values
    for i in range(grid_res):
      for j in range(grid_res):
        point = np.array([X[i, j], Y[i, j]])
        f_val = eigenvector_f(point, m) #value of f
        pdf_val = pdf(point) #value of pdf
        integrand_vals[i, j] = f_val * pdf_val
    
    integrand_vals[~mask] = 0 # set values outside the disk to zero

    return np.sum(integrand_vals) * dx * dy #sum over all values and multiply by area of each cell
    #NOTE when one wants to improve precision look at np.sum(pdf_values, dtype=np.float128) or use math.fsum
  
  C_4 = integral(m, pdf, grid_res)
  return C_4

if __name__ == "__main__":
  m = 16
  #c = constant_c(m)
  #print(f"Constant c for m={m}: {c}")

  e0 = vector([1] + [0]*(m-1)) #constant 1 is (1,0,...,0) in the basis (x-a)^i
  m = 16
  print('f*[1] = ', dual_f_star(e0, m))
  print('f[1] = ', eigenvector_f(e0, m))

  print('Constant C_4 = ', Constant_C_4(k=4, m=m, pdf=?))
  

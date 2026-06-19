from sage.all import *
from DominantEigenvalue.DominantEigenvalue import matrix_T_m
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

a= QQ(1)/2 #1/2 as element in Q
x = var('x')
precision = 200
CF = ComplexField(precision)

def compute_eigenvectors(m, plot=False):
  A = matrix_T_m(m)
  data_left = A.eigenvectors_left()
  data_right = A.eigenvectors_right()

  lambda_dom_r, ev_right, mult_r = max(data_right, key=lambda data: abs(data[0].real())) #right eigenvector corresponding to the dominant eigenvalue
  lambda_dom_l, ev_left, mult_l = max(data_left, key=lambda data: abs(data[0].real())) #left eigenvector corresponding to the dominant eigenvalue

  assert abs(lambda_dom_r - lambda_dom_l) < 1e-10, "Dominant eigenvalues for left and right eigenvectors do not match!"
  
  ev_right = ev_right[0]
  ev_left = ev_left[0]

  #real vectors
  ev_right = vector([x.real() for x in ev_right])
  ev_left  = vector([x.real() for x in ev_left])

  f0 = sum(ev_right[i] * (-0.5)**i for i in range(len(ev_right)))
  assert f0 != 0, "f(0) is zero, cannot normalize eigenvector!"
  print(f"f(0) before normalization: {f0}")
  
  #scale right eigenvector such that f(0) = 1
  ev_right = ev_right / f0

  if plot:
    x_vals = np.linspace(0, 1, 100)
    f_vals = [sum(ev_right[i] * (x - 0.5)**i for i in range(len(ev_right))) for x in x_vals]
    
    plt.plot(x_vals, f_vals, label='Eigenvector f(x)')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title('Right Eigenvector as a Function of x')
    plt.legend()
    plt.grid()
    plt.show()

  #normalize the eigenvectors such that <ev_left, ev_right> = 1 i.e dual
  norm = ev_left.dot_product(ev_right)
  print(f"<f*,f> = {norm}")
  assert norm != 0, "Left and right eigenvectors are orthogonal, cannot normalize!"

  ev_left = ev_left / norm
  
  return ev_left, ev_right

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

def dual_f_star(x, m, ev_left=None):
  if ev_left is None:
    ev_left, _ = compute_eigenvectors(m)
  
  x_sage = vector(CC, x)
    
  if abs(imag_part(ev_left.dot_product(x_sage))) > 1e-10:
    print("Warning: significant imaginary part:",ev_left.dot_product(x_sage))

  return real_part(ev_left.dot_product(x_sage))

def eigenvector_f(w, m, ev_right=None):
  if ev_right is None:
    _, ev_right = compute_eigenvectors(m)

  val = 0
  for i in range(m):
    val += ev_right[i] * (w - 1/2)**i

  if abs(imag_part(val)) < 10**(-precision/2):
    return real_part(val)

  return val

def Constant_C_4(k, m, grid_res=300, plot=False):
  # C_4 = \int_D \int_[0,1] (y - y^2) f_4(x + i(1-2y)t) dy dA(w) where w = x + it
  def integral(m, grid_res):
    #grid of points for the integral dA(w)
    x = np.linspace(0, 1, grid_res)
    t = np.linspace(-0.5, 0.5, grid_res)
    X, T = np.meshgrid(x, t)
    #circle mask
    mask = ((X - 0.5)**2 + T**2) <= 0.25
    #space between points
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    #grid for the integral dy
    Y = np.linspace(0, 1, grid_res)
    dy = Y[1] - Y[0]

    #only compute ev_right once not every time in the loop
    _, ev_right = compute_eigenvectors(m)

    total = 0.0 + 1j*0.0

    I = np.zeros_like(X, dtype=float) #for plot

    for i in range(grid_res):
      for j in range(grid_res):
        if not mask[i, j]:
          continue #skip points outside the disk
        inner_sum = 0.0 + 1j*0.0
        for y_k in Y:
          s = (1 - 2*y_k)
          w = X[i,j] + 1j * s * T[i,j]
          f_val = eigenvector_f(w, m, ev_right) #value of f_4 at the point x + i(1-2y)t
          inner_sum += f_val * (y_k - y_k**2) * dy
        I[i,j] = (inner_sum).real #for plot
        total += inner_sum
    
    I[~mask] = np.nan #set points outside the disk to nan for better plotting
    return 6 * total * dx * dt, (X,T,I,dx,dt)

  C_4, (X,T,I,dx,dt) = integral(m, grid_res)

  #Plot
  if plot: 
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')

    if True:
      ax.plot_surface(X, T, I, cmap='viridis')
      surf = ax.plot_surface(X, T, I, cmap='viridis')
      fig.colorbar(surf, ax=ax, shrink=0.6, aspect=10)
    else:
      ax.bar3d(
        X.flatten(),
        T.flatten(),
        np.zeros_like(I).flatten(),
        dx,
        dt,
        I.flatten(),
        shade=True
      )
    
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('I(x,t)')
    ax.set_title('Averaged Integrand for $C_4$')

    plt.tight_layout()
    plt.show()

  #return result
  if abs(imag_part(C_4)) > 1e-10:
    print("Warning: significant imaginary part in C_4:", C_4)
  print("C_4 (complex):", C_4)
  with open(f"C_4_value_{m}.txt", "w") as f:
    f.write(f"C_4 = {real_part(C_4)}\n")
  return real_part(C_4)

if __name__ == "__main__":
  m = 32
  #c = constant_c(m)
  #print(f"Constant c for m={m}: {c}")

  e0 = vector([1] + [0]*(m-1)) #constant 1 is (1,0,...,0) in the basis (x-a)^i
  m = 64
  #print('f*[1] = ', dual_f_star(e0, m))
  #print('f[1] = ', eigenvector_f(0.2 + 0j, m))

  
  print('Constant C_4 = ', Constant_C_4(k=4, m=m, grid_res=400, plot=False))

  #print(f"eigenvectors for m={m}: ", compute_eigenvectors(m, plot=True))
  

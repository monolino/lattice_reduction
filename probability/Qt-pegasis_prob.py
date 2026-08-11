from sage.all import *
from ConstantC.constant import dual_f_star, Constant_C_4
from DominantEigenvalue.DominantEigenvalue import matrix_T_m
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def get_samples(filename):
  samples = []
  
  with open(filename, 'r') as f:
    for line in f:
      line = line.strip().strip('()') #of the form (x1, x2)
      x1, x2 = map(float, line.split(','))
      samples.append(x1 + 1j * x2)
  samples = np.array(samples)
  return samples

#========================= PDF function and F(w,/overline{w}) = pdf(w) =========================#

def pdf(z):
  # d=0
  #d = 1e-9
  #d = 2.5027632403871487e-05 #mean
  d = 1.0620646970966825e-09 #lower bound
  #z is a complex number
  return np.exp(- (z.imag - d) ** 2 / ( 2* 0.05 ** 2)) #sigma=0.05

def F(u,v):
  return  np.exp(-((v-u)/(2 * 1j))**2 * (1/(2*0.05**2))) #sigma=0.05

def normalization_constant_N(pdf,m, samples=100, plot=False):
  #vector with basis (x-1/2)^i
  #compute the coefficients by numerical integration
  
  #sample points in D
  x = np.linspace(0, 1, samples)
  t = np.linspace(-0.5, 0.5, samples)
  X, T = np.meshgrid(x, t)

  mask = (X - 0.5)**2 + T**2 <= 0.25
  Z = X + 1j*T

  dx = x[1] - x[0]
  dt = t[1] - t[0]

  
  points = Z[mask]
  values = pdf(points)

  values_plot = pdf(Z)
  values_plot[~mask] = np.nan


  if plot:
    
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(X, T, values_plot, cmap='viridis')

    fig.colorbar(surf, ax=ax, shrink=0.6)

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$t$')
    ax.set_zlabel(r'$\rho(x+it)$')
    ax.set_title('PDF over the disk $D$')

    ax.set_box_aspect([1,1,0.5])

    plt.tight_layout()
    plt.show()


  #normalize pdf so that the integral over D is 1
  Z_norm = np.sum(values) * dx * dt

  if Z_norm == 0:
    raise ValueError("Integral of pdf over D is zero, cannot normalize.")

  print("Normalization constant N = ", (1/Z_norm), Z_norm)
  return 1/Z_norm

#KDE
def K_H_kde(u,v):
  sigma_x = 0.2
  sigma_y = 0.02
  return np.exp(-1/2 * ((u + v)**2 / (4*sigma_x) + (u - v)**2 / (-4*sigma_y)))

def pdf_kde(x, sample_points):
  n = len(sample_points)
  sol = 0
  for i in range(n):
    z = x - sample_points[i]
    sol += K_H_kde(z, z.conjugate())
  return 1/n * sol

def F_kde(u,v, sample_points):
  n = len(sample_points)
  sol = 0 
  for i in range(n):
    sol += K_H_kde(u - sample_points[i],v - sample_points[i].conjugate())
  return 1/n * sol

def f_kde(u, sample_points):
  return F_kde(u,u,sample_points)

def f_two_variables(x,y,sample_points):
  u = x + 1j* y
  return f_kde(u, sample_points)

'''helper functions'''
def basis(z, m):
  terms = []
  for i in range(m):
    terms.append((z-1/2)**i)
  return np.array(terms)

def basis_functions(m):
  return np.array( lambda x, i=i: (x - 1/2) **i for i in range(m)) # does not work like intended? maybe make a matrix out of it and also needs a two-dimensional basis function ?


def f_to_vector(function, m, sample_points, samples=100, A=None, f_vals=None):
  
  def get_A_f_points(sample_points, samples, m):
    '''sparse grid'''
    x = np.linspace(0, 1, samples)
    y = np.linspace(0, 1, samples)
    X, Y = np.meshgrid(x, y)
    mask = (X - 0.5)**2 + Y**2 <= 0.25
    Z = X + 1j * Y

    grid_points = Z[mask]

    f_vals = np.array([f_kde(z, sample_points) for z in grid_points])
    print(f_vals)

    A = np.vstack([basis(z, m) for z in grid_points])
    return A, f_vals, grid_points
  
  if A is None and f_vals is None:
    A, f_vals, _ = get_A_f_points(sample_points, samples, m)
    with open(f"Matrix_A_fVal_m_{m}.txt", "w") as f:
      f.write(f"{A}\n{f_vals}\n")

  'solve Ac = f'
  coeffs, *_ = np.linalg.lstsq(A, f_vals, rcond=None)
  print(coeffs)
  return coeffs

def f_to_vector_any_function(function, m, samples=100,A=None, f_vals=None):
  
  def get_A_f_points(samples, m):
    '''sparse grid'''
    x = np.linspace(0, 1, samples)
    y = np.linspace(0, 1, samples)
    X, Y = np.meshgrid(x, y)
    mask = (X - 0.5)**2 + Y**2 <= 0.25
    Z = X + 1j * Y

    grid_points = Z[mask]

    f_vals = np.array([function(z) for z in grid_points])
    print(f_vals)

    A = np.vstack([basis(z, m) for z in grid_points])
    return A, f_vals, grid_points
  
  if A is None and f_vals is None:
    A, f_vals, _ = get_A_f_points(samples, m)
    with open(f"Matrix_A_fVal_m_{m}_pdf.txt", "w") as f:
      f.write(f"{A}\n{f_vals}\n")

  'solve Ac = f'
  coeffs, *_ = np.linalg.lstsq(A, f_vals, rcond=None)
  print(coeffs)
  return coeffs



#============================= Probability parts =============================#

def spectral_gap(m):
  lam_4 = 0.1994588183437672601918456859798790
  mu_4 = 0.081780
  T = matrix_T_m(m, s=6)
  lam_6, _, _ = max(T.eigenvectors_right(), key=lambda data: abs(data[0].real()))
  lam_6 = lam_6.real()
  print(f"m={m}: lambda_4 = {lam_4}, lambda_6 = {lam_6}, mu_4 = {mu_4}")
  return 1/lam_4 * max(lam_6, mu_4)

def norm_F(m, F, sample_points,samples=200):
  x = np.linspace(0, 1, samples)
  y = np.linspace(-0.5, 0.5, samples)
  X, Y = np.meshgrid(x, y)
  a = np.linspace(0, 1, samples)
  b = np.linspace(-0.5, 0.5, samples)
  A, B = np.meshgrid(a, b)

  mask_U = (X - 0.5)**2 + Y**2 <= 0.25
  U = X + 1j*Y
  mask_V = (A - 0.5)**2 + B**2 <= 0.25
  V = A + 1j*B

  u_points = U[mask_U]
  v_points = V[mask_V]

  
  n = len(u_points)

  idx_u = np.random.choice(n, 1000)
  idx_v = np.random.choice(n, 1000)

  


  #N = normalization_constant_N(pdf, m)
  values = F(u_points[idx_u], u_points[idx_v], sample_points)
  max_val = np.max(np.abs(values))
  print(f"||F|| = {max_val}")
  return max_val

def probability_estimate(k, m, pdf, grid_res=100, C_4=None, plot=False):
  lam = 0.1994588183437672601918456859798790
  N = normalization_constant_N(pdf, m, samples=grid_res, plot=plot)
  f =  vector(CC, [N] + [0]*(m-1)) #F(u,u) = f(u) = N * e^0 = N where N is normalization constant of the pdf
  dual_f_p = dual_f_star(f, m)
  if C_4 is None:
    C_4 = Constant_C_4(k, m, grid_res)
  return lam ** k * dual_f_p * C_4  #+ small term


def prob_estimate_kde(k, m, sample_points, C_4=None, grid_res=100,multi=1):
  lam = 0.1994588183437672601918456859798790
  if C_4 is None:
    C_4 = Constant_C_4(k, m, grid_res)
  f_vector = f_to_vector(f_kde, m , sample_points)
  dualval = dual_f_star(f_vector, m)
  proba = lam**k * dualval * C_4
  with open(f"prob_estimate_multi_{multi}_m_{m}.txt", "w") as f:
    f.write(f"k= {k}\nlam= {lam}\nC_4= {C_4}\nf*[f]= {dualval}\nP[L>=k+1]= {proba}")
  return proba

def f_dual_of_pdf(m, pdf):
  #d = 1e-9
  #d = 2.5027632403871487e-05 #mean
  #d = 1.0620646970966825e-09 #lower bound
  N = normalization_constant_N(pdf, m)
  normalized_pdf = lambda z: N * pdf(z)
  pdf_vector = f_to_vector_any_function(normalized_pdf, m)
  f_dual = dual_f_star(pdf_vector, m)
  return f_dual



if __name__ == "__main__":
  #k = 30
  m = 128
  #with open(f"C_4_value_{m}.txt", "r") as f:
  #  C_4 = float(f.read().strip().split()[2])
  print("START")
  #prob = probability_estimate(k=k, m=16, pdf=pdf, C_4=C_4, grid_res=300, plot=False)
  #normF = norm_F(m, pdf)
  #gap = spectral_gap(m)**k
  #print(f"P[L >= {k} + 1] = {prob} * (1 + O({normF * gap}))")
  #print(f"error = {prob * normF * gap}")
  #norm_F(m, pdf)
  #print(spectral_gap(m))

  #multi=10
  #m = 64

  #sample_points = get_samples("probability/points_-21474836479_samples_10000_multi_10.txt")

  #prob = prob_estimate_kde(k,m,sample_points, multi=multi)
  #normF = norm_F(m, F_kde, sample_points=sample_points)
  #gap = spectral_gap(m)**k
  #print(f"P[L >= {k} + 1] = {prob} * (1 + O({normF * gap}))")
  #print(f"error = {prob * normF * gap}")

  print(f_dual_of_pdf(m, pdf))

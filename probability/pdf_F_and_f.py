from sage.all import *
from ConstantC.constant import dual_f_star, Constant_C_4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def pdf(z):
  d=0 #concentrated at y=d
  #z is a complex number
  return np.exp(- ((z.conjugate() - z)/(2 * 1j)+d) ** 2 * 1 / (  2* 0.05 ** 2))  #sigma=0.05

def plot_pdf_p(pdf,m, samples=100):
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
  Z_norm = np.sum(values) * dx * dt

  values_plot = pdf(Z) * (1/Z_norm)
  values_plot[~mask] = np.nan

    
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

  print("Normalization constant N = ", (1/Z_norm))
  print("Integral of pdf over D = ", Z_norm)

def F(u,v):
  return np.exp(- ((v - u)/(2 * 1j)) ** 2 * 1 / (  2* 0.05 ** 2))  #sigma=0.05

def f(u):
  return F(u,u)

def plot_f(m, grid_res=100):
  #sample points in D
  x = np.linspace(0, 1, grid_res)
  t = np.linspace(-0.5, 0.5, grid_res)
  X, T = np.meshgrid(x, t)

  mask = (X - 0.5)**2 + T**2 <= 0.25
  Z = X + 1j*T

  dx = x[1] - x[0]
  dt = t[1] - t[0]

  values_plot = f(Z)
  values_plot[~mask] = np.nan

    
  fig = plt.figure(figsize=(8,6))
  ax = fig.add_subplot(111, projection='3d')

  surf = ax.plot_surface(X, T, values_plot, cmap='viridis')
  fig.colorbar(surf, ax=ax, shrink=0.6)

  ax.set_xlabel(r'$x$')
  ax.set_ylabel(r'$t$')
  ax.set_zlabel(r'$f(x+it)$')
  ax.set_title('Diagonal function f(u) = F(u,u) over the disk $D$')

  ax.set_box_aspect([1,1,0.5])
  plt.tight_layout()
  plt.show()

if __name__ == "__main__":
  m = 16
  plot_pdf_p(pdf, m, samples=300)
  plot_f(m, grid_res=300)
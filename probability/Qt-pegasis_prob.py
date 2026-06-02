from sage.all import *
from ConstantC.constant import dual_f_star, Constant_C_4
import numpy as np

def pdf(x):
  #x is np.array and second argument is imaginary part of z
  return np.exp(- x[1] ** 2 / ( 2* 0.01 ** 2)) #sigma=0.01

def probability_estimate(k, m, pdf, grid_res=100):
  lam = 0.1994588183437672601918456859798790
  dual_f_1 = dual_f_star(vector([1] + [0]*(m-1)), m)
  C_4 = Constant_C_4(k, m, pdf, grid_res)
  return lam ** k * dual_f_1 * C_4 #+ small term

if __name__ == "__main__":
  k = 80
  print("what")
  print(f"P[L >= {k} + 1] = ", probability_estimate(k=k, m=16, pdf=pdf))
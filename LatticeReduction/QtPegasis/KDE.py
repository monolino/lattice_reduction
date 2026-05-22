import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import sqrtm #square root of a matrix

#Kernel function x = np.array([x1, x2]) is a 2D vector, H is a 2x2 positive definite matrix (covariance matrix)
def K(x, H_det, H_inv):
  return (2*np.pi)**(-1) * H_det**(-1/2) * np.exp(-0.5 * x.T @ H_inv @ x)

#H = np.array([[1, 0], [0, 1]]) #TODO find good H
#H_det = np.linalg.det(H)
#H_inv = np.linalg.inv(H)

def pdf(x, H_det, H_inv, samples):
  return 1/len(samples) * sum(K(x - sample, H_det, H_inv) for sample in samples)

def get_samples(filename):
  samples = []
  
  with open(filename, 'r') as f:
    for line in f:
      line = line.strip().strip('()') #of the form (x1, x2)
      x1, x2 = map(float, line.split(','))
      samples.append(np.array([x1, x2]))
  samples = np.array(samples)
  samples_scaled = samples.copy()
  samples_scaled[:,1] *= 1e49
  return samples_scaled

def get_H(samples): #TODO find good H, for now use Scott's rule
  n = len(samples)
  Sigma = np.cov(samples.T)
  H = n**(-1/6) * sqrtm(Sigma) * 0.99 #H^ = n^*(-1/(d+4)) * Sigma^1/2
  return H, np.linalg.det(H), np.linalg.inv(H)

def pdf_1D(x, samples, h):
  x_vals = samples[:,0]
  return 1/ (len(samples) * h) * sum(np.exp(-0.5 * ((x - x_i)/h)**2) for x_i in x_vals)

def get_h(samples, factor=1.0):
  n = len(samples)
  sigma = np.std(samples[:,0])
  h = factor * 1.06 * sigma * n**(-1/5) #Silverman's rule of thumb for 1D
  return h

def plot_pdf(samples, H_det, H_inv):
  x = np.linspace(-1, 1, 100)
  y = np.linspace(-1, 1, 100)
  X, Y = np.meshgrid(x, y)
  Z = np.zeros_like(X)

  for i in range(X.shape[0]):
    for j in range(X.shape[1]):
      Z[i, j] = pdf(np.array([X[i, j], Y[i, j]]), H_det, H_inv, samples)

  plt.figure(figsize=(6,6))
  plt.contourf(X, Y, Z, levels=50, cmap='viridis')
  plt.colorbar(label='Estimated PDF')
  plt.title('Estimated PDF from Samples (y scaled by 1e49)')
  plt.xlabel('x1')
  plt.ylabel('x2')
  plt.xlim(-1, 1)
  plt.ylim(-1, 1)
  plt.gca().set_aspect('equal')
  plt.show()

def plot_pdf_1D(samples, h):
  x = np.linspace(-1, 1, 100)
  y = [pdf_1D(x_i, samples, h) for x_i in x]

  plt.figure(figsize=(6,4))
  plt.plot(x, y, label='Estimated PDF')
  plt.title('Estimated PDF from Samples (1D)')
  plt.xlabel('x')
  plt.ylabel('Density')
  plt.legend()
  plt.show()

if __name__ == "__main__":
  D = 5 * 2**32 - 1
  num_samples = 10000
  samples = get_samples(f"points_{-D}_samples_{num_samples}.txt")
  #samples = get_samples('test_points.txt')
  #H, H_det, H_inv = get_H(samples)
  #print(f"H:\n{H}\nH_det: {H_det}\nH_inv:\n{H_inv}")
  #plot_pdf(samples, H_det, H_inv)

  factor = 0.5
  h = get_h(samples, factor)
  plot_pdf_1D(samples, h)
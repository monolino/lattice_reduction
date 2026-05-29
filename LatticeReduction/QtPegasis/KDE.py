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
  return 1/ (len(samples) * h * np.sqrt(2*np.pi)) * sum(np.exp(-0.5 * ((x - x_i)/h)**2) for x_i in x_vals)

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
  x = np.linspace(0, 1, 100)
  y = [pdf_1D(x_i, samples, h) for x_i in x]

  plt.figure(figsize=(6,4))
  plt.plot(x, y, label='Estimated PDF')
  plt.title('Estimated PDF from Samples (1D)')
  plt.xlabel('x')
  plt.ylabel('Density')
  plt.legend()
  plt.savefig("pdf_1D.png")
  plt.show()

def gaussian_2D(x, mean, cov_inv, cov_det):
  n = 2 #2D
  diff = x - mean #as np.array
  exponent = -0.5 * diff.T @ cov_inv @ diff
  return (2 * np.pi) ** (-n / 2) * cov_det ** (-0.5) * np.exp(exponent)

def holomorphic_gaussian_2D(z, sigma, holomorphic=True):
  if holomorphic:
    return np.exp(- z ** 2 / (  sigma ** 2))
  return np.exp(- np.imag(z) ** 2 / ( 2* sigma ** 2))


def plot_gaussian_2D(mean, cov):
  x = np.linspace(0, 1, 100)
  y = np.linspace(-0.5, 0.5, 100)
  X, Y = np.meshgrid(x, y)
  Z = np.zeros_like(X)
  mask =  (X - 0.5)**2 + Y**2 <= 0.25 #disk circle equation

  cov_inv = np.linalg.inv(cov)
  cov_det = np.linalg.det(cov)

  for i in range(X.shape[0]):
    for j in range(X.shape[1]):
      Z[i, j] = gaussian_2D(np.array([X[i, j], Y[i, j]]), mean, cov_inv, cov_det)

  Z[~mask] = np.nan # only inside disk

  plt.figure(figsize=(6,6))

  theta = np.linspace(0, 2*np.pi, 300)
  circle_x = 0.5 + 0.5 * np.cos(theta)
  circle_y = 0.5 * np.sin(theta)


  plt.plot(circle_x, circle_y, 'r', label='Re(1/z) = 1')
  
  plt.contourf(X, Y, Z, levels=50, cmap='viridis')
  plt.colorbar(label='Gaussian PDF')
  plt.title('2D Gaussian PDF')
  plt.xlabel('x1')
  plt.ylabel('x2')
  plt.xlim(-0.5, 1.5)
  plt.ylim(-1, 1)
  plt.gca().set_aspect('equal')
  plt.show()

def plot_holomorphic_gaussian_2D(sigma, holomorphic=True):
  x = np.linspace(0, 1, 400)
  y = np.linspace(-0.5, 0.5, 400)
  X, Y = np.meshgrid(x, y)

  # build complex grid
  Z_complex = X + 1j * Y

  # evaluate function
  sigma = 0.03
  Z = holomorphic_gaussian_2D(Z_complex, sigma, holomorphic=holomorphic)

  # mask (your circle!)
  mask = (X - 0.5)**2 + Y**2 <= 0.25
  Z[~mask] = np.nan


  # plot
  plt.figure(figsize=(6,6))
  if not holomorphic:
    plt.contourf(X, Y, Z, levels=50, cmap='viridis')
    plt.clim(0, 1)
  else:
    Z_abs = np.abs(Z)
    plt.contourf(X, Y, Z_abs, levels=100, cmap='viridis')
    plt.clim(0, 1)
  
  plt.colorbar()
  #plt.contour(X, Y, Z, levels=[0.05, 0.2, 0.5], colors='white')

  # draw circle boundary
  theta = np.linspace(0, 2*np.pi, 300)
  circle_x = 0.5 + 0.5*np.cos(theta)
  circle_y = 0.5*np.sin(theta)
  plt.plot(circle_x, circle_y, 'r')

  plt.gca().set_aspect('equal')
  plt.xlabel("Re(z)")
  plt.ylabel("Im(z)")
  plt.xlim(-0.5, 1.5)
  plt.ylim(-1, 1)
  
  
  plt.title(
    "Gaussian-like concentration near Im(z)=0\n"
    + ("f(z) = exp(-z^2 / sigma^2), sigma=" + str(sigma) if holomorphic 
    else "f(z) = exp(-Im(z)^2 / (2*sigma^2)), sigma=" + str(sigma))
  )

  plt.show()

if __name__ == "__main__":
  D = 5 * 2**32 - 1
  #num_samples = 1000
  #samples = get_samples(f"points_{-D}_samples_{num_samples}.txt")
  #samples = get_samples('test_points.txt')
  #H, H_det, H_inv = get_H(samples)
  #print(f"H:\n{H}\nH_det: {H_det}\nH_inv:\n{H_inv}")
  #plot_pdf(samples, H_det, H_inv)

  #factor = 1.0
  #h = get_h(samples, factor)
  #plot_pdf_1D(samples, h)
  mean = np.array([0.5, 0]) 
  
  cov = [
    [0.16, 0],    # σx = 0.2  (wide)
    [0, 0.0004]   # σy = 0.02 (very thin)
  ] 

  plot_gaussian_2D(mean, cov)
  #plot_holomorphic_gaussian_2D(sigma=0.01, holomorphic=False)
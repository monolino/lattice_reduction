from sage.all import *
from L_rdc_algorithm import *
from quaternion import QuaternionAlgebra
from euclideanspace import EuclideanSpace
import random
import matplotlib.pyplot as plt

def small_generator(n, QA):
  '''I is left ideal, n is the integer n=nrd(I) and b is the basis of I
  returns alpha
  '''
  while True:
    coeffs = [random.randint(1, 10**5) for _ in range(4)]
    if gcd(2*coeffs[0], n) == 1 and gcd(QA.nrd(coeffs), n**2) == 1:
      return coeffs

def lattice_vectors(n, QA):
  coeffs = small_generator(n, QA)
  u = (QQ(n), 0) #n is gcd so integer
  v = (QQ(-2*coeffs[1]), QQ(2* coeffs[0])) # coeffs are all integer
  return u,v

def histogram(n, QA, num_samples=1000, log=False):
  #NOTE cannot take v as the key since the same lattice has different basis v.
  E = EuclideanSpace(2)

  p = QA.p
  if p == 5 * 2**248 - 1:
    p_str_title = "5 * 2^248 - 1"
    p_str = "p1"
  elif p == 65 * 2**376 - 1:
    p_str_title = "65 * 2^376 - 1"
    p_str = "p2"
  elif p == 27 * 2**500 - 1:
    p_str_title = "27 * 2^500 - 1"
    p_str = "p3"

  histogram = {}
  for _ in range(num_samples):
    v, u = lattice_vectors(n, QA)
 
    (vm, um), L = lattice_reduction_2dim(E, v, u) #from here now euclidean norm.
    if L not in histogram:
      histogram[L] = 1
    else:
      histogram[L] += 1
  
  # Sort by L
  L_values = sorted(histogram.keys())
  counts = [histogram[L] for L in L_values]

  #log
  if log:
    histogram_file = f"histogram_{p_str}_{num_samples}.txt"
    with open(histogram_file, "w") as f:
      for L in L_values:
        f.write(f"L = {L}: {histogram[L]}\n")

  

  plt.figure()
  bars = plt.bar(L_values, counts)

  #Show count on top of bars
  for bar in bars:
    height = bar.get_height()
    plt.text(
      bar.get_x() + bar.get_width()/2,
      height,
      str(int(height)),
      ha='center',
      va='bottom',
      fontsize=8
    )
  
  plt.yscale('log')

  plt.xlabel("Number of reduction steps L")
  plt.ylabel("Frequency")
  
  plt.title(
    f"Histogram of lattice reduction steps\n"
    f"Quaternion algebra B(-1,{p_str_title}), samples = {num_samples}"
  )

  plt.savefig(f"histogram_{p_str}_samples_{num_samples}.png")

  return histogram
  

if __name__ == "__main__":
  p1 = 5 * 2**248 - 1
  p2 = 65 * 2**376 - 1
  p3 = 27 * 2**500 - 1
  B = QuaternionAlgebra(p1)
  n = 1001
  u,v = lattice_vectors(n, B)
  E = EuclideanSpace(2)
  print(f"u: {u}, v: {v}")

  histo = histogram(n, B, num_samples=10000)
  print(histo)
  

  

  


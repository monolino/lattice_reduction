from sage.all import *
import random as pyrandom #Sage also has a random
from L_rdc_qtpegasis import *
import matplotlib.pyplot as plt
import numpy as np

#how to hash into class groups Construction 2 p.11 of https://eprint.iacr.org/2024/034.pdf


#------------Prime generation----------------#
#code from here https://github.com/seresistvanandras/hashingToClassGroups/blob/main/ClassGroupPlayground%2025.ipynb
def random_prime_length(n):
    prime_candidate = 2
    while True:
            prime_candidate = getLowLevelPrime(n)
            if not isMillerRabinPassed(prime_candidate):
                continue
            else:
                break
    return prime_candidate

# Pre generated primes
first_primes_list = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                     31, 37, 41, 43, 47, 53, 59, 61, 67,
                     71, 73, 79, 83, 89, 97, 101, 103,
                     107, 109, 113, 127, 131, 137, 139,
                     149, 151, 157, 163, 167, 173, 179,
                     181, 191, 193, 197, 199, 211, 223,
                     227, 229, 233, 239, 241, 251, 257,
                     263, 269, 271, 277, 281, 283, 293,
                     307, 311, 313, 317, 331, 337, 347, 349]

def nBitRandom(n):
  return pyrandom.randrange(2**(n-1)+1, 2**n - 1)

def getLowLevelPrime(n):
    '''Generate a prime candidate divisible 
    by first primes'''
    while True:
        # Obtain a random number
        pc = nBitRandom(n)
 
        # Test divisibility by pre-generated
        # primes
        for divisor in first_primes_list:
            if pc % divisor == 0 and divisor**2 <= pc:
                break
        else:
            return pc

def isMillerRabinPassed(mrc):
    '''Run 30 iterations of Rabin Miller Primality test'''
    maxDivisionsByTwo = 0
    ec = mrc-1
    while ec % 2 == 0:
        ec >>= 1
        maxDivisionsByTwo += 1
    assert(2**maxDivisionsByTwo * ec == mrc-1)
 
    def trialComposite(round_tester):
        if pow(round_tester, ec, mrc) == 1:
            return False
        for i in range(maxDivisionsByTwo):
            if pow(round_tester, 2**i * ec, mrc) == mrc-1:
                return False
        return True
 
    # Set number of trials here
    numberOfRabinTrials = 30
    for i in range(numberOfRabinTrials):
        round_tester = pyrandom.randrange(2, mrc)
        if trialComposite(round_tester):
            return False
    return True

def bit_length_n(D):
  '''input discriminant D approx 2^n, output n/2 -1'''
  n = abs(D).bit_length() # D approx 2^n
  k = (n // 2) - 1
  return k

def sage_random_prime(n):
  p = random_prime(2**n - 1) # lbound=2**(n-1)+1
  return p

#------------Sampling class group----------------#
'''what is the discriminant D?'''
''' qt-pegasus paper uses D, see the primes they denote and that are actually the d from sqrt -d hmm TODO unsure if that means the discriminant is negative or positive??
d1 = 3 * 11 * 2 ** 503 - 1
d2 = 3 * 5 * 2 ** 1004 - 1
d3 = 3 * 3 * 2 ** 1551 - 1
d4 = 3 * 17 * 2 ** 2026 - 1
d5 = 3 * 3 * 7 * 2 ** 4084 - 1
'''




#Wesolowski’s hash-to-class group construction (not hash but random sampling)
def random_class_group_element(D, multi=1):
  assert D % 4 == 1 or D % 4 == 0, "D should be congruent to 1 mod 4 or 0 mod 4"
  n = (bit_length_n(D)) * multi
  #rejection sampling
  while True:
    p = sage_random_prime(n)
    if legendre_symbol(D, p) == 1: #p splits in Q(sqrt{D})
      D_mod_p = D % p
      b = modular_sqrt(D_mod_p, p) #solutions b, p - b in mod p

      #valid binary quadratic form with disciminant D, (p,b,c)
      numerator = b*b - D
      if numerator % (4*p) != 0:
        continue
      c = numerator // (4*p)
      if c == 0 or b*b - 4*p*c != D:
        continue

      break

  #not hashing (s = p mod 3), but randomly choosing
  sign = pyrandom.randint(0,1)
  if sign:
    return  (p, b)
  return  (p, -b)

def legendre_symbol(a, p):
  """ Compute the Legendre symbol a|p using
      Euler's criterion. p is a prime, a is
      relatively prime to p (if p divides
      a, then a|p = 0)
      Returns 1 if a has a square root modulo
      p, -1 otherwise.
  """
  ls = pow(a, (p - 1) // 2, p)
  return -1 if ls == p - 1 else ls

def modular_sqrt(a, p):
  """ Find a quadratic residue (mod p) of 'a'. p
      must be an odd prime.
      Solve the congruence of the form:
        x^2 = a (mod p)
      And returns x. Note that p - x is also a root.
      0 is returned is no square root exists for
      these a and p.
      The Tonelli-Shanks algorithm is used (except
      for some simple cases in which the solution
      is known from an identity). This algorithm
      runs in polynomial time (unless the
      generalized Riemann hypothesis is false).
  """
  if legendre_symbol(a, p) != 1:
    return 0
  elif a == 0:
    return 0
  elif p == 2:
    return p
  elif p % 4 == 3:
    return pow(a, (p + 1) // 4, p) #pow(base, exp, mod)
  
  s = p - 1
  e = 0
  while s % 2 == 0:
    s //= 2
    e += 1
  
  n = 2
  while legendre_symbol(n, p) != -1:
    n += 1

  x = pow(a, (s + 1) // 2, p)
  b = pow(a, s, p)
  g = pow(n, s, p)
  r = e

  while True:
    t = b
    m = 0
    for m in range(r):
      if t == 1:
        break
      t = pow(t, 2, p)

    if m == 0:
      return x

    gs = pow(g, 2 ** (r - m - 1), p)
    g = (gs * gs) % p
    x = (x * gs) % p
    b = (b * g) % p
    r = m


#------------from class group to vector----------------#
#NOTE: working over basis (1,sqrt{D}) not (1,i)
def class_group_element_to_vector_in_QQ(element):
  p, b = element
  v = (QQ(p), QQ(0))
  u = (-QQ(b)/2, QQ(1)/2)
  assert all(isinstance(x, (int, Integer, Rational)) for x in v + u), "v and u should lie in QQ^2"
  return v, u

#----------------Histogram----------------#
def histogram_of_lattice_reduction(D, num_samples, log=False, multi=10):
  if -D == 3 * 11 * 2 ** 503 - 1:
    d_str = "d1"
  elif -D == 3 * 5 * 2 ** 1004 - 1:
    d_str = "d2"
  elif -D == 3 * 3 * 2 ** 1551 - 1:
    d_str = "d3"
  elif -D == 3 * 17 * 2 ** 2026 - 1:
    d_str = "d4"
  elif -D == 3 * 3 * 7 * 2 ** 4084 - 1:
    d_str = "d5"
  elif -D == 5 * 2**32 - 1:
    d_str = "5 * 2**32 - 1"
  else:  d_str = f"D={D}"

  histogram = {}

  not_minimal = 0

  for _ in range(num_samples):
    class_group_element = random_class_group_element(D, multi=multi)
    
    if False:
      cglog = f"classgroupelementslog_{d_str}_{num_samples}_multi_{multi}.txt"
      with open(cglog, "a") as f:
        f.write(f"p = {class_group_element[0]}, b = {class_group_element[1]}\n")

    v, u = class_group_element_to_vector_in_QQ(class_group_element)

    
    try:
      test_minimal_basis(v, u, D)
      is_minimal = True
    except AssertionError:
      is_minimal = False
      not_minimal += 1


    (vm, um), L = lattice_reduction_2dim(v, u, D)
    if L not in histogram:
      histogram[L] = 1
    histogram[L] += 1

  # Sort by L
  L_values = sorted(histogram.keys())
  counts = [histogram[L] for L in L_values]

  #log
  if log:
    histogram_file = f"histogram_{d_str}_{num_samples}_multi_{multi}.txt"
    with open(histogram_file, "w") as f:
      for L in L_values:
        f.write(f"L = {L}: {histogram[L]}\n")

  

  plt.figure()
  plt.bar(L_values, counts)

  plt.xlabel("Number of reduction steps L")
  plt.ylabel("Frequency")
  

  plt.title(
    f"Histogram of lattice reduction steps\n"
    f"Discriminant D = {d_str}, samples = {num_samples}"
  )

  plt.savefig(f"histogram_{d_str}_samples_{num_samples}_multi_{multi}.png")
  plt.show()

  print(f"Not minimal basis count: {not_minimal} out of {num_samples}")
  return histogram

def plot_z_in_disk(D, num_samples=1000, multi=1, log=False):  
  minimal_count = 0
  swapped = 0
  points = []
  for _ in range(num_samples):
    class_group_element = random_class_group_element(D, multi=multi)
    v, u = class_group_element_to_vector_in_QQ(class_group_element)
    try :
      test_minimal_basis(v, u, D)
      minimal_count += 1
    except AssertionError:
      pass
    
    #NOTE z needs to be in B = { z | 0 <= Re(z) <= 1} and in the disk D = { z | Re(1/z) >= 1}
    #z = v/u
    if ((v[0]*u[0] + abs(D)*v[1]*u[1]) / (v[0]**2 + abs(D)*v[1]**2)) < 1: #Re(1/z) < 1 i.e z is outside the disk
      swapped += 1
      v, u = u, v 

    z = ((v[0]*u[0] + abs(D)*v[1]*u[1]) / (u[0]**2 + abs(D)*u[1]**2), (abs(D))**(0.5)*(v[1]*u[0] - v[0]*u[1]) / (u[0]**2 + abs(D)*u[1]**2))

    # z is now such that Re(1/z) >= 1
    if z[0] < 0 or z[0] > 1:
     z = (z[0] - math.floor(z[0]), z[1]) #shift to get back in the fundamental domain
      
    assert z[0] >= 0 and z[0] <= 1, f"Re(z) should be in [0,1], got {z[0]}"
    points.append(z)

  print(f"Minimal basis count: {minimal_count} out of {num_samples}")
  print(f"Swapped: {swapped} out of {num_samples}")
  xs, ys = zip(*points)

  plt.figure(figsize=(6,6))

  #log points
  if log:
    with open(f"points_{D}_samples_{num_samples}_multi_{multi}.txt", "w") as f:
      for z in points:
        f.write(f"({float(z[0])}, {float(z[1])})\n")
  # plot the points
  plt.scatter(xs, ys, s=3, alpha=1, color='black')
  

  # draw the circle (x-1/2)^2 + y^2 = (1/2)^2
  theta = np.linspace(0, 2*np.pi, 300)
  circle_x = 0.5 + 0.5 * np.cos(theta)
  circle_y = 0.5 * np.sin(theta)


  plt.plot(circle_x, circle_y, 'r', label='Re(1/z) = 1')
  
  #plt.axhline(0)
  #plt.axvline(0)
  
  plt.xlim(-0.5, 1.5)
  plt.ylim(-1, 1)

  #plt.gca().set_aspect('equal', adjustable='box')
  plt.gca().set_aspect('equal')


  plt.title(
    f"Points z = v/u in the complex plane\n"
    f"Class group with discriminant D = {D}, samples = {num_samples}"
  )
  plt.suptitle(f"Swapped: {swapped} out of {num_samples}")
  plt.xlabel("Re(z)")
  plt.ylabel("Im(z)")
  plt.legend()

  plt.savefig(f"Disk_{D}_samples_{num_samples}_multi_{multi}.png")

  plt.show()

def multiple_histograms(D, num_samples, multi_list):
  plt.figure()
  for multi in multi_list:
    histogram = {}
    for _ in range(num_samples):
      class_group_element = random_class_group_element(D, multi=multi)
      v, u = class_group_element_to_vector_in_QQ(class_group_element)
      (vm, um), L = lattice_reduction_2dim(v, u, D)
      if L not in histogram:
        histogram[L] = 1
      histogram[L] += 1
    
    L_values = sorted(histogram.keys())
    counts = [histogram[L] for L in L_values]
    
    
    plt.bar(
    L_values,
    counts,
    alpha=0.5,                #transparency
    label=f"exp times {multi}" 
    )
    print(f"multi={multi} ploted")

  plt.xlabel("Number of reduction steps L")
  plt.ylabel("Frequency")

  plt.title(
    f"Histogram of lattice reduction steps\n"
    f"Discriminant D = 5 * 2**32 - 1, samples = {num_samples}"
  )
  plt.legend()
  plt.tight_layout()
  plt.savefig(f"multiple_histogram_samples_{num_samples}.png")
  plt.show()
  return plt

  
if __name__ == "__main__":
  d1 = 3 * 11 * 2 ** 503 - 1
  d2 = 3 * 5 * 2 ** 1004 - 1
  d3 = 3 * 3 * 2 ** 1551 - 1
  d4 = 3 * 17 * 2 ** 2026 - 1
  d5 = 3 * 3 * 7 * 2 ** 4084 - 1

  D = 5 * 2**32 - 1

  #NOTE: the discriminant D is -d, so -d1, -d2, -d3, -d4, -d5
  if False:
    class_group_element = random_class_group_element(-d1)
    print(f"Random class group element: {class_group_element}")
    v, u = class_group_element_to_vector_in_QQ(class_group_element)
    print(f"Corresponding vector in QQ^2: v={v}, u={u}")

  multi = 10
  #histogram = histogram_of_lattice_reduction(-D, 10000, log=True, multi=multi)
  plot_z_in_disk(-D, num_samples=10000, log=True, multi=multi)
  #multiple_histograms(-D, num_samples=1000, multi_list=[10, 11, 12, 13])

  if False:
    p,b = class_group_element
    q = - (p * b // 2)
    print(f"q = {q}")


 
  


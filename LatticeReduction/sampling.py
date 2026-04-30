from sage.all import *
import random

#how to hash into class groups Construction 2 p.11 of https://eprint.iacr.org/2024/034.pdf


#------------Prime generation----------------#
#code from here https://github.com/seresistvanandras/hashingToClassGroups/blob/main/ClassGroupPlayground%2025.ipynb
def random_prime(n):
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
  return random.randrange(2**(n-1)+1, 2**n - 1)

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
        round_tester = random.randrange(2, mrc)
        if trialComposite(round_tester):
            return False
    return True

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
def random_class_group_element(D, n):
  #rejection sampling
  while True:
    p = random_prime(n)
    if legendre_symbol(D, p) == 1: #p splits in Q(sqrt{D})
      D_mod_p = D % p
      b = modular_sqrt(D_mod_p, p)

      #valid binary quadratic form with disciminant D, (p,b,c)
      numerator = b*b - D
      if numerator % (4*p) != 0:
        continue
      c = numerator // (4*p)
      if c == 0 or b*b - 4*p*c != D:
        continue

      break

  #not hashing (s = p mod 3), but randomly choosing
  sign = randint(0,1) 
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
  v = (QQ(p), 0)
  u = (-QQ(b)/2, QQ(1)/2)
  return v, u
  
  
if __name__ == "__main__":
  D = -2545453607 #cool discriminant (?)
  n = 16 #depend on discriminant see paper hashing
  class_group_element = random_class_group_element(D,n)
  print(f"Random class group element: {class_group_element}")
  v, u = class_group_element_to_vector_in_QQ(class_group_element)
  print(f"Corresponding vector in QQ^2: v={v}, u={u}")

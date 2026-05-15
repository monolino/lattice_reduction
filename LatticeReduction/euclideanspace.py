class EuclideanSpace:
    def __init__(self, n=2):
      self.n = n
    
    def inner_product(self, v, u):
      res = 0
      for i in range(self.n):
        res += v[i] * u[i]
      return res
    
    def norm(self, u):
      return self.inner_product(u,u)

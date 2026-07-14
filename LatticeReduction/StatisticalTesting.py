import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

set_p = np.linspace(0.01, 0.99, 99)


def read_data(filename):
  data = []
  with open(filename, "r") as f:
    for line in f:
      _, y = line.split(":")
      y = float(y.strip())
      data.append(y)
  data = np.array(data)

  #normalize such that sum is 1
  data /= np.sum(data)
  print("Sum of normalized data:", np.sum(data))
  return sorted(data)


def get_q_p(data, p):
  sum_pk = 0
  k = 0
  while sum_pk < p:
    sum_pk += data[k]
    k += 1
  return k

def quantiles(data, p_set):
  q_p = []
  for p in p_set:
    q_p.append(get_q_p(data, p))
  return np.array(q_p)

if __name__ == "__main__":
  filename1 = "Total_Length_Lineintersection_d_2.5027632403871487e-05.txt"
  filename2 = "Total_Length_Lineintersection_d_1e-09.txt"
  filename4 = "Total_Length_Lineintersection_d_1.0620646970966825e-09.txt"
  filename3 = "QtPegasis/histogram_5 * 2**32 - 1_10000_multi_3 copy.txt"
  data = read_data(filename1)
  print(data)
  q_p = quantiles(data, set_p)
  print("Quantiles for p values:", q_p)
  stats.probplot(q_p, dist="norm", plot=plt)
  plt.title("Normal Q-Q Plot")
  plt.show()
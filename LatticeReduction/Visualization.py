from sage.all import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects



class Disk:
  def __init__(self, center, radius, k):
    self.center = center
    self.radius = (radius)
    self.k = k
    assert center[1] == 0
    assert k >= 0
  
  def intersection_line(self, d):
    #d = y-value of the line i.e L_d : y = d
    #return [a,b] interval (only x values since y = d trivial)
    if abs(d) > self.radius:
      return None #no intersection with disk
    
    sqrt_term =  np.sqrt(self.radius**2 - d**2)

    a_x = float(self.center[0] - sqrt_term)
    b_x = float(self.center[0] + sqrt_term)
    return (a_x, b_x)

  def plot(self, ax, color='black', alpha=0.5, label=None):
    theta = np.linspace(0, 2*np.pi, 300)
    center_x, center_y = self.center
    
    x = center_x + self.radius * np.cos(theta)
    y = center_y + self.radius * np.sin(theta)

    ax.fill(x, y, color=color, alpha=alpha, label=label)


class Domain_D_k:
  def __init__(self, k):
    self.k = k
  
  def __str__(self):
    return f"D_{self.k}"

  
  def __repr__(self):
    return f"D_{self.k}"


  @staticmethod
  def apply_word(word, z):
    for m in reversed(word):
      z = 1 / (m + z)
    return z

  @staticmethod
  def generate_words( k, m_max):
    if k == 0: return [[]]
    smaller = Domain_D_k.generate_words(k-1, m_max)
    words = []
    for w in smaller:
      for m in range(1, m_max+1):
        words.append(w + [m])
    return words

  def next_word(self, word, d, m_max):
    w = word[:]  # copy

    for i in range(len(w) - 1, -1, -1):
      
      if w[i] == m_max:
        w[i] = 1 #reset
      else:
        w[i] += 1
        word_length = self.check_word_radius(w,d)
        if word_length is None:
          w[i] = 1
        else: return word_length #found valid next word
    return None #came back to word [1,1,1,...,1] so with this m_max no further disks

  def check_word_radius(self, w, d):
    a, b = self.apply_word(w, 0), self.apply_word(w, 1)
    radius = abs(a - b) / 2
    center = abs(a + b) / 2
    #log
    with open(f"Word_log_k_{self.k}.txt", "a") as f:
      f.write(f"[{min(a,b):.6f}, {max(a,b):.6f}], radius = {radius:.6f}, word = {w}\n")
    disk = Disk(center, radius, self.k)
    interval = disk.intersection_line(d)
    if interval is None: return None
    length_interval = interval[1] - interval[0]
    return length_interval

  def next_word_2(self, word, d):
    w = word[:]
    w_new = w+[1]
    res = self.check_word_radius(w_new, d)
    #if new word intersects with line return, word, length of intersection, disk
    if res is not None: return w_new, res
    #if new word does not intersect with line, try to find next word
    w_new_new = w[:]
    w_new_new[-1] += 1
    res = self.check_word_radius(w_new_new, d)
    if res is not None: return w_new_new, res
    #if new word does not intersect with line, try to find next word
    w_new_new_new = w[:-1]
    #check if finished i.e again []
    if w_new_new_new == []: return None
    w_new_new_new[-1] += 1
    res = self.check_word_radius(w_new_new_new, d)
    if res is not None: return w_new_new_new, res
    return self.next_word_2(w_new_new_new, d) #recursively find next word

        
  def intervals_domain_m(self, m_max):
    #from right to left ! e.g ((1/2,1), (1/3,1/2), (1/4,1/3), ..., (1/(m+1), 1/m))
    #domains_set = []
    #for i in range(1,m+1):
    #  domains_set.append((1/(i+1), 1/i))
    #return tuple(domains_set)
    words = Domain_D_k.generate_words(self.k, m_max)
    intervals = []

    for w in words:
      a = Domain_D_k.apply_word(w, 0)
      b = Domain_D_k.apply_word(w, 1)
      intervals.append((min(a, b), max(a, b)))
    return intervals

  def union_disks(self, m_max):
    intervals = self.intervals_domain_m( m_max)
    return [Disk(center = ((a+b)/2, 0), radius= (b-a)/2, k=self.k) for (a,b) in intervals]

  def plot_domain(self, ax, color='black', d=None, label=None):
    if d is not None:
      m_max = self.find_m(d)
    else: m_max = 20
    set_of_disks = self.union_disks(m_max)

    for i,d in enumerate(set_of_disks):
      if i == 0:
        d.plot(ax, color=color, alpha=0.5, label=label)
      else:
        d.plot(ax, color=color, alpha=0.5)
  
  def biggest_disk_in_domain(self):
    return self.union_disks(1)[0]

  def find_m(self, d):
    m = 1
    while True:
      word = [m for _ in range(self.k)]
      a, b = self.apply_word(word,0), self.apply_word(word,1)
      radius = abs( a - b)/2
      if radius < d:
        return m

def plot_until_k(max_k):
  fig, ax = plt.subplots()

  colors = plt.cm.viridis_r

  #Plot D_0
  D_0 = Disk(center=(0.5, 0), radius=0.5, k=0)
  D_0.plot(ax, color=colors(0), alpha=1, label=f"k={0}")

  #Plot D_k
  for k in range(1, max_k+1):
    color = colors(k / (max_k + 1))
    D_k = Domain_D_k(k)
    D_k.plot_domain(ax, color=color, label=f"k={k}")
    print(f"plot for k {k}")
  
  ax.legend(loc="upper right")
  ax.set_aspect('equal')
  ax.set_xlim(0, 1)
  ax.set_ylim(-0.5, 0.5)

  plt.savefig(f"Domains_k_{max_k}.png")
  plt.show()

def line_through_domains(d):
  if d > 0.5:
    raise Warning(f'The line y = {d} does not intersect with the disk of diameter [0,1]')
  elif d == 0.5:
    raise Warning(f'The line y = {d} does only intersect the disk at one point')
  set_of_domains = []
  k = 0
  D_k = Domain_D_k(k) #initial D_0
  while abs(d) < D_k.biggest_disk_in_domain().radius:
    set_of_domains.append(D_k)
    k += 1
    D_k = Domain_D_k(k)
  return set_of_domains

#---Helper function---#
def interval_dividing(dictionary ,list1, k, max_k):
  if k == max_k:
    return dictionary
  changed_dictionary = dictionary.copy()
  list2 = dictionary[k+1]
  new_list_1 = list1.copy()
  for (a, b) in list2:
    updated = []
    for (c, d) in new_list_1:
      if a >= c and b <= d:
        updated.append((c,a))
        updated.append((b,d))
      else:
        updated.append((c,d))
    new_list_1 = updated
      
  changed_dictionary[k] = new_list_1

  return interval_dividing(changed_dictionary, list2, k+1, max_k)

def line_intervals_domains(d,store_intervals=True):
  'Intervals = {k : [(a,b), (c,d), (e,f), ...], etc.}'
  Intervals_raw = {} if store_intervals else None
  Length_k = {}
  set_of_domains = line_through_domains(d)
  print(set_of_domains) 
  max_k = max(domain.k for domain in set_of_domains)

  for D_k in set_of_domains:
    length_k = 0
    m_max = D_k.find_m(d)
   #TODO
    unions_disks = D_k.union_disks(m_max)

    list_intervals = []
    for disks in unions_disks:
      interval_disk = disks.intersection_line(d) #(a,b)
      if interval_disk is not None:
        if store_intervals: list_intervals.append(interval_disk) #[(a,b)]
        length = interval_disk[1] - interval_disk[0]
        assert length >= 0
        length_k += length


    if store_intervals: Intervals_raw[D_k.k] = list_intervals
    Length_k[D_k.k] = length_k
    print(f"k={D_k.k}, {len(list_intervals)}")
  if store_intervals: Intervals_final = interval_dividing(Intervals_raw, Intervals_raw[0], 0, max_k)
  else: Intervals_final = {}

  return Intervals_final, Length_k

def compute_length_list_intervals(li):
  length = 0
  for intervals in li:
    minus = intervals[1] - intervals[0]
    assert minus >= 0
    length += minus
  return length

def plot_line_domains(d, plot_domains=False):
  Intervals,_ = line_intervals_domains(d)
  fig, ax = plt.subplots()

  ax.axhline(y=d, color='black', linestyle='--', linewidth=1) #line L_d

  colors = plt.cm.viridis_r

  max_k = max(Intervals.keys())

  for k, intervals in Intervals.items():
    color = colors(k / (max_k + 1))

    #plot the domains
    if plot_domains:
      D_k = Domain_D_k(k)
      color = colors(k / (max_k + 1))
      D_k.plot_domain(ax, color=color, d=d)

    #plot the intervals
    for i, (a, b) in enumerate(intervals):
      label = f"k={k}" if i == 0 else None
      # plot interval as thick line segment
      ax.plot([a, b], [d, d], linewidth=6, color=color, label=label)

      # compute length
      length = b - a

      # label at midpoint
      mid = (a + b) / 2
      text = ax.text(mid, d + 0.02*(k+1), f"{length:.3f}", ha='center', color=color, fontsize=8)
      

      text.set_path_effects([
        path_effects.Stroke(linewidth=2, foreground='white'),
        path_effects.Normal()
      ])
  y_offset = -0.15 #below the plot
  for k, intervals in Intervals.items():
    total_length = compute_length_list_intervals(intervals)

    ax.text(0.02, y_offset - 0.07*k,
      f"k={k}: length = {total_length:.3f}",
      transform=ax.transAxes,
      fontsize=9,
      verticalalignment='top',
      )



  ax.legend()
  ax.set_title(f"Line intersection with domains (y = {d})")
  ax.set_xlim(0, 1)
  ax.set_ylim(d - 0.2, d + 0.2)
  ax.set_aspect('equal')

  plt.show()

def total_length_intervals(d):
  Domains = line_through_domains(d)
  total_length = {}

  for D_k in Domains: #D_0, D_1, D_2 , ...
    k = D_k.k
    m_max = 5
    word = tuple(1 for _ in range(k)) #word that corresponds to largest disk
    length_k = 0
    radius = 1
    while radius >= d: #disks largest, ..., smallest that intersects
      a, b = apply_word(word, 0), apply_word(word, 1) #this is the interval which is the diameter of the disk
      a, b = min(a,b), max(a,b)
      radius = (b-a)/2
      disk = Disk(center = ((a+b)/2, 0), radius=radius, k=k)
      interval_disk = disk.intersection_line(d)
      if interval_disk is not None:
        length_k += interval_disk[1] - interval_disk[0]
      word = next_word(word, m_max)
    total_length[k] = length_k
    if k != 0:
      total_length[k-1] -= length_k 

def histogram_for_k(d, plot=True):
  _, Length_dict = line_intervals_domains(d, store_intervals=False)
    
  if not plot:
    return Length_dict

  x = list(Length_dict.keys())
  y = list(Length_dict.values())
  
  plt.bar(x, y)
  plt.xlabel("k")
  plt.ylabel("Length")
  plt.title("Histogram of lengths by k")

  plt.savefig(f"Histogram_lines_length_d{d}")
  plt.show()

def plot_all_disk_intersecting_line(d):
  k = self.k
  word = [1] * (k - 1) + [0] #next_word computes [1,1,1,...,1] the actual first word

  while new_word != None:
    new_word, length, disk = self.next_word(word, d, m_max)
    if new_word is not None:
      word = new_word
      print(f"word = {word}, length = {length}, disk = {disk}")
    else:
      print("No more words found")
      break






def apply_word(word, z):
  for m in reversed(word):
    z = 1 / (m + z)
  return z

def check_word_radius(w, d):
  a, b = apply_word(w, 0), apply_word(w, 1)
  radius = abs(a - b) / 2
  center_x = abs(a + b) / 2
  disk = Disk((center_x, 0), radius, k=len(w))
  interval = disk.intersection_line(d)
  if interval is None: return None
  length_interval = interval[1] - interval[0]
  return w, length_interval, disk


def next_dfs_word(word, d):
  #RETURNS NONE if done otherwise returns (word, (length of intersection, disk))
  assert d != 0, "d cannot be 0, because infinite tree"
  w = word[:]
  child = w+[1]
  res = check_word_radius(child, d)
  #if new word intersects with line return, word, length of intersection, disk
  if res is not None: return child, res
  #if new word does not intersect with line, try to find next word
  sibling = w[:]
  if sibling == []: return None
  sibling[-1] += 1
  res = check_word_radius(sibling, d)
  if res is not None: return sibling, res
  #if new word does not intersect with line, try to find next word
  parent_sibling = w[:-1]
  #check if finished i.e again []
  if parent_sibling == []: return None
  parent_sibling[-1] += 1
  res = check_word_radius(parent_sibling, d)
  if res is not None: return parent_sibling, res
  return next_dfs_word(parent_sibling, d) #recursively find next word


def plot_all_disk_intersecting_line(d, plot=True):
  word = []
  length_k = {}
  fig, ax = plt.subplots()
  colors = plt.cm.viridis_r

  #plot the line y = d
  ax.axhline(y=d, color='black', linestyle='--', linewidth=1, label=f"y = {d}") #line L_d
  #plot unit disk
  D_0 = Disk(center=(0.5, 0), radius=0.5, k=0)
  D_0.plot(ax, color=colors(0), alpha=1, label=f"k={0}")
  length_k[0] = D_0.intersection_line(d)[1] - D_0.intersection_line(d)[0] if D_0.intersection_line(d) is not None else 0


  while True:
    res = next_dfs_word(word, d)
    if res is None: break
    word, (_, length, disk) = res
    k = len(word)
    length_k[k] = length_k.get(k, 0) + length
    color = colors(k / 10)  # Adjust the denominator for color scaling as needed
    if word == [1] * k: label = f"k={k}"
    else: label = None
    if plot: disk.plot(ax, color=color, alpha=0.5, label=label)

  #recompute lengths for each k
  total_length = {}
  for k in length_k.keys():
    total_length[k] = length_k[k] - length_k.get(k+1, 0)

  if not plot:
    return total_length
   
  #write the results
  text = "\n".join(
    f"k={k}: length = {total_length[k]:.3f}"
    for k in sorted(total_length)
  )

  fig.text(0.15, 0.02, text, fontsize=9,va='bottom')
  plt.subplots_adjust(bottom=0.3)
  plt.title(f"Disks intersecting line y = {d}")
  ax.set_aspect('equal')
  ax.set_xlim(0, 1)
  ax.set_ylim(-0.5, 0.5)
  ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1))
  plt.savefig(f"Disks_intersecting_line_d_{d}.png")
  plt.show()

def total_length_intersection(d):
  word = []
  length_k = {}
  D_0 = Disk(center=(0.5, 0), radius=0.5, k=0)
  length_k[0] = D_0.intersection_line(d)[1] - D_0.intersection_line(d)[0] if D_0.intersection_line(d) is not None else 0
  while True:
    res = next_dfs_word(word, d)
    if res is None: break
    word, (_, length, disk) = res
    k = len(word)
    length_k[k] = length_k.get(k, 0) + length

  #recompute lengths for each k
  total_length = {}
  for k in length_k.keys():
    total_length[k] = length_k[k] - length_k.get(k+1, 0)

  for k in sorted(total_length):
    print(f"k={k}: length = {total_length[k]:.3f}")

  return total_length


def plot_histogram_for_thesis():
  k = np.arange(20)
  #lengths = [0.008,0.045,0.125,0.215,0.248,0.198,0.109,0.041,0.010,0.001,0] # mean value of lengths for d = 2.5027632403871487e-05
  #lengths = [0.000,0.001,0.003,0.012,0.032,0.067,0.113,0.155,0.175,0.163,0.127,0.081,0.043,0.019,0.007,0.002,0.000,0.000,0.000,0.000] # infimum value of lengths for d = 1.0620646970966825e-09
  lengths = [0.000,0.000,0.002,0.006,0.019,0.043,0.081,0.124,0.158,0.168,0.151,0.114,0.072,0.038,0.017,0.006,0.002,0.000,0.000,0.000] # lowerbound value of lengths for d = 1.8626451493176932e-10

  #k = (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
  #lengths = [9, 85, 443, 1193, 2192, 2641, 2112, 1048, 258, 27, 3] # sampled lattices for multi=3


  plt.figure(figsize=(8,5))

  plt.bar(k, lengths, width=0.85, color="#1f77b4", edgecolor="white")

  plt.title(r"Lower bound value ($d = 1.86\cdot 10^{-10}$)", fontsize=16)
  plt.xlabel("Number of reduction steps $k$", fontsize=14)
  plt.ylabel(r"Number of Lattices", fontsize=14)
  plt.xticks(k)
  plt.grid(axis='y', alpha=0.3)
  plt.tight_layout()
  #plt.savefig("histogram_reduction_steps_meanvalue.png", dpi=300)
  plt.show()

def plot_skew_vs_minimal_basis():

  b1 = np.array([3,0])
  b2 = np.array([1,3])

  pts = []
  for m in range(0,6):
    for n in range(0,6):
        pts.append(m*b1+n*b2)

  fig, ax = plt.subplots(figsize=(4,4))
  pts = np.array(pts)

  ax.scatter(pts[:,0], pts[:,1], color='black', s=10)
  b1 = np.array([3,0])
  b2 = np.array([10,3])
  ax.arrow(0, 0, b1[0], b1[1],color='black',width=0.03,head_width=0.25,head_length=0.35,length_includes_head=True)
  ax.arrow(0, 0, b2[0], b2[1],color='black',width=0.03,head_width=0.25,head_length=0.35,length_includes_head=True)
  ax.set_aspect('equal')

  ax.axis('off')
  plt.show()


if __name__ == "__main__":
  #plot_until_k(5)
  D_1 = Domain_D_k(0)
  #print(D_1.biggest_disk_in_domain().center)
  d = 0.05#2e-10#5e-44

  d = 0.00714 #2857142857143
  #print( line_intervals_domains(d))
  #plot_line_domains(0.04, plot_domains=True)
  #histogram_for_k(d)

  #d=0
  #d = 1e-9
  d = 2.5027632403871487e-05 #mean
  #d = 1.0620646970966825e-09 #lower bound
  #d = 1.8626451493176932e-10 #proved lower bound

  #plot_all_disk_intersecting_line(d=1.0620646970966825e-09, plot=False)
  #total_length_intersection(d)

  #plot_histogram_for_thesis()
  plot_skew_vs_minimal_basis()






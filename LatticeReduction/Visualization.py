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
    a, b = apply_word(w, 0), apply_word(w, 1)
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
    m = 5 #can be adjustable when one knows the d will be small so in many disks. Then take bigger m to start with
    
    while True:
      union_disks = self.union_disks(m)

      if (disks.radius < abs(d) for disks in union_disks):
        return m
      m += 5

  def blabla(self, d, m_max):
    k = self.k
    word = [1] * (k - 1) + [0] #next_word computes [1,1,1,...,1] the actual first word
    word, length = self.next_word(word, d, m_max) if not None else return 0
    while next_word is not None:
      word, newlength = self.next_word(word, d, m_max) if not None else word, newlength = None, 0
      length += newlength



    

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





if __name__ == "__main__":
  #plot_until_k(5)
  D_1 = Domain_D_k(0)
  #print(D_1.biggest_disk_in_domain().center)
  d = 0.05#2e-10#5e-44
  #print( line_intervals_domains(d))
  #plot_line_domains(d, plot_domains=False)
  histogram_for_k(d)
  
    






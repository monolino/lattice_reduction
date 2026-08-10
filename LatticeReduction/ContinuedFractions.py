from itertools import product
import matplotlib.pyplot as plt
import numpy as np

def apply_word(word, z):
  for m in reversed(word):
    z = 1 / (m + z)
  return z

def radius_of_word(word):
  a, b = apply_word(word,0), apply_word(word,1)
  radius = abs( a - b)/2
  center = abs(a+b)/2
  return radius, center

def sort_results(results):
  return sorted(results, key=lambda x: x["radius"])

def continued_fractions(m_max, k):
  results = []
  for word in product(range(1, m_max + 1), repeat=k):
    radius, center = radius_of_word(word)
    results.append({
      "word": list(word),
      "radius": radius,
      "center": center
    })
  return results

def interval(word):
  a, b = apply_word(word, 0), apply_word(word, 1)
  return min(a,b), max(a,b)

def plot_disks():
  alpha = 1
  fig, ax = plt.subplots()

  #colors = plt.cm.viridis_r

  #Fundamental Disk
  center=(0.5, 0)
  radius=0.5
  theta = np.linspace(0, 2*np.pi, 300)
  center_x, center_y = center
    
  x = center_x + radius * np.cos(theta)
  y = center_y + radius * np.sin(theta)

  ax.fill(x, y, color="#FFF7AE", alpha=alpha, label=f"[]")

  #Disks [m]
  for m in range(1,7):
    radius, center = radius_of_word([m])
    theta = np.linspace(0, 2*np.pi, 300)
    center_x, center_y = center, 0
    
    x = center_x + radius * np.cos(theta)
    y = center_y + radius * np.sin(theta)

    ax.fill(x, y, color="#B8E6B8", alpha=alpha, label=f"[{m}]")

  #Disks [1,m]
  for m in range(1,7):
    radius, center = radius_of_word([1,m])
    theta = np.linspace(0, 2*np.pi, 300)
    center_x, center_y = center, 0
    
    x = center_x + radius * np.cos(theta)
    y = center_y + radius * np.sin(theta)

    ax.fill(x, y, color="#BFD7FF", alpha=alpha, label=f"[1,{m}]")\

  
  #Disks [2,m]
  for m in range(1,7):
    radius, center = radius_of_word([2,m])
    theta = np.linspace(0, 2*np.pi, 300)
    center_x, center_y = center, 0
    
    x = center_x + radius * np.cos(theta)
    y = center_y + radius * np.sin(theta)

    ax.fill(x, y, color="grey", alpha=alpha, label=f"[2,{m}]")
  
  ax.legend(loc="upper left", bbox_to_anchor=(1.05, 1))
  ax.set_aspect('equal')
  ax.set_xlim(0, 1)
  ax.set_ylim(-0.5, 0.5)

  plt.savefig(f"Disks_corresponding_to_words.png")
  plt.show()

def max_value_in_specific_word(word, d):
  #word = [m_1, m_2, ..., m_{k-1}, m_k]
  #returns value of how many disks one needs to consider in the disk [m_1, m_2, ..., m_{k-1}]
  index = len(word) - 1  # m[index] = m_k
  while True:
    radius, _ = radius_of_word(word)
    if radius < d:
      return word[index] - 1
    word[index] += 1

def length_of_intersection_line(radius, d):
  if radius < d:
    return 0
  else:
    return float(2 * np.sqrt(radius**2 - d**2))

def Total_Length_Lineintersection(d):
  #NOTE improvement possible return radius in next_word already when computed
  stored_raw_results = {}
  #Fundamental Disk
  length = length_of_intersection_line(0.5, d)
  stored_raw_results[0] = length
  #==================================================================
  def next_word(word,d):
    if word == []: #START at root node
      return [1]

    #valid child
    child = word + [1]
    if radius_of_word(child)[0] >= d: 
      return child

    #No valid child, go to next sibling
    w = word.copy() #just to be save
    w[-1] += 1
    while True: #next higher part of tree
      if radius_of_word(w)[0] >= d: #found valid sibling
        return w
          
      del w[-1] #remove depth
        
      if w == []: #Again at root node so finished
        return None

      w[-1] += 1 #goes to next sibling above
  #==================================================================

  #ROOT
  word = []
  while True:
    #CREATE NEXT WORD
    word = next_word(word,d)
    print(word)
    #IF again at ROOT then finished
    if word is None:
      break
    k = len(word)
    #COMPUTE LENGTH
    radius, _ = radius_of_word(word)
    length = length_of_intersection_line(radius, d)
    #ADD LENGTH
    if k not in stored_raw_results:
      stored_raw_results[k] = length
    else:
      stored_raw_results[k] += length

  total_length = {}
  for k in stored_raw_results:
    total_length[k] = stored_raw_results[k] - stored_raw_results[k+1] if k+1 in stored_raw_results else stored_raw_results[k]
  
  with open(f"Total_Length_Lineintersection_d_{d}.txt", "w") as f:
    for k in sorted(total_length.keys()):
      f.write(f"{k}: {total_length[k]}\n")
  return total_length

def make_histogram_from_file(filename):
  k = []
  values = []
  with open(filename, "r") as f:
    for line in f:
      x, y = line.split(":")
      k.append(int(x))
      values.append(float(y))
  plt.figure(figsize=(8, 4))
  
  plt.bar(k, values, width=1.0, align="center", edgecolor="black")

  plt.xlabel("k")
  plt.ylabel("Value")
  plt.title("Histogram")
  plt.xticks(k)

  plt.savefig(f"Histogram")
  plt.show()


def get_mean_y_value(filename):
  y_values = []
  with open(filename, "r") as f:
    for line in f:
      _, y = line.split(",")
      y, _ = y.split(")")
      val = abs(float(y.strip())) #because other wise it could be 0 the mean
      y_values.append(val)
  mean_y = sum(y_values) / len(y_values)
  return mean_y

def get_smallest_y_value(filename):
  y_values = []
  with open(filename, "r") as f:
    for line in f:
      _, y = line.split(",")
      y, _ = y.split(")")
      val = abs(float(y.strip()))
      y_values.append(val)
  return min(y_values)

#result = continued_fractions(20,1)

#sortedresult = sort_results(result)

#d = 1e-9
#d = 2.5027632403871487e-05 #mean
d = 1.0620646970966825e-09 #lower bound

#for r in sortedresult:
  #print(r["word"], r["radius"])

radius = 0.25 #m=1

#for m in range(2,21):
  #radius = radius * (m-1)/(m+1)
  #print(radius)


#plot_disks()

#print(max_value_in_specific_word([10000,1]))

#print(Total_Length_Lineintersection(d))

#make_histogram_from_file("Total_Length_Lineintersection_d_1e-09.txt")
make_histogram_from_file("Total_Length_Lineintersection_d_2.5027632403871487e-05.txt")

#print(get_mean_y_value("QtPegasis/points_-21474836479_samples_10000_multi_3.txt"))
#print(get_smallest_y_value("QtPegasis/points_-21474836479_samples_10000_multi_3.txt"))
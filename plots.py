from sage.all import *
from hurwitzZeta import *

def plot_lambdas(s,M):
  lambdas = []
  for m in [4,8,16]:
    lambdas.append(lambda_estimate(s,M,m))
  print(lambdas)
    
  P = list_plot(
    list(zip([4,8,16], lambdas)),
    plotjoined=True,
    marker='o',
    title= 'Lambdas for different m',
    axes_labels=(r"$m$", r"$\lambda_m$")
  ) 
  P.save("lambda_vs_m.png")
  print("Saved plot")
  return

def plot_eigenfunctions():
  p_4 = p(4)
  p_8 = p(8)
  p_16 = p(16)
  p_32 = p(32)
  p_64 = p(64)
  plt_4 = plot(
    p_4,
    (0, 1),
    color= 'blue',
    legend_label=r"$p_4$",
    plot_points=200,
    axes_labels=[r"$t$", r"$p_m(t)$"],
    title=r"$Eigenfunctions for different m$"
  )
  plt_8 = plot(
    p_8,
    (0, 1),
    color= 'red',
    legend_label=r"$p_8$",
    plot_points=200,
  )
  plt_16 = plot(
    p_16,
    (0, 1),
    color= 'orange',
    legend_label=r"$p_{16}$",
    plot_points=200,
  )
  plt = plt_4 + plt_8 + plt_16
  plt.save('eigenfunctions_vs_m.png')
  print("Saved plot")
  return

def plot_Gf_f(s,M,m):
  eig_function = lambda x: real_part(p(m)(x))
  def G_eval_div_f_numeric_m(t):
    return G_eval(eig_function, s, M, t)/eig_function(t)
  
  plt = plot(
    G_eval_div_f_numeric_m,
    (0, 1),
    plot_points=200,
    ymin=0.194,
    ymax=0.225,
    axes_labels=[r"$t$", r"$G[f](t)/f(t)$"],
    title=r"$M = 1000,\; s = 4$"
  )
  plt.save(f"G_plot_own_m{m}.png")
  print("Saved plot")

def plot_eigenfunction_m(m):
  eig_function = lambda x: real_part(p(m)(x))
  plt = plot(
    eig_function,
    (0, 1),
    plot_points=200,
    axes_labels=[r"$t$", r"$p_m(t)$"],
    title=r"$Eigenfunction$"
  )
  plt.save(f"eigenfunction_p_m{m}.png")
  print("Saved plot")
  return

if __name__ == "__main__" :
  #plot_lambdas(4,1000)
  #plot_eigenfunctions
  plot_Gf_f(4,1000,8)
  #plot_eigenfunction_m(8)

  '''
  def G_eval_div_f_numeric(t):
    g = lambda t: 1-87/50*(t-a)+ 391/200*(t-a)**2 - 1687/1000*(t-a)**3
    return G_eval(g, 4, 1000, t)/g(t)
  plt = plot(
    G_eval_div_f_numeric,
    (0, 1),
    plot_points=200,
    axes_labels=[r"$t$", r"$G[f](t)/f(t)$"],
    title=r"$M = 1000,\; s = 4$"
  )
  plt.save("G_plot.png")
  print("Saved plot to G_plot.png")
  '''
    
  '''
  p1 = plot(
    eigen_function,
    (0, 1),
    color='blue',
    thickness=2,
    legend_label=r"$f_{\mathrm{computed}}(x)$",
    plot_points=300
  )

  p2 = plot(
    their_eigfunc,
    (0, 1),
    color='red',
    linestyle='--',
    thickness=2,
    legend_label=r"$f_{\mathrm{paper}}(x)$",
    plot_points=300
  )

  plt_eig =  p1 + p2  
  plt_eig.save('eigenfunc_plot.png')
  '''
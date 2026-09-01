This is the code base for my master thesis "An Analysis of the Number of Reduction Steps for Special Two-Dimensional Lattices".

# Abstract
The Gaussian algorithm for lattice reduction in two dimensions is analyzed for special input lattices arising from ideals in imaginary quadratic fields. In particular, we examine the distribution of the number of reduction steps. For this, we examine the tail probability $\Pr[L\ge k+1]$, which quantifies the probability that the reduction algorithm requires more than $k$ reduction steps.


We analyze the geometric properties of such lattices arising from ideals and construct a probabilistic model for their initial distribution in a fundamental domain. Based on the transfer operator framework introduced by  Daudé, Flajolet, and Vallée [[DFV96](https://inria.hal.science/inria-00073892/file/RR-2798.pdf)], we derive numerical estimates for the reduction step probability of these special lattices. We further introduce a geometric line-intersection model that exploits the geometric properties observed in the initial distribution of the lattices. A numerical comparison with experimental data shows that this method provides a more accurate prediction of the observed probabilities.

# File Tree

```
.
├── README.md
├── ConstantC/
|    ├── constant.py
|    └── ...
├── DominantEigenvalue/
|    ├── DominantEigenvalue.py
|    ├── plots.py
|    └── ...
├── LatticeReduction/
|    ├── Log_and_plots/
|    |    └── ...
|    ├── QtPegasis/
|    |    ├── Log_and_plots/
|    |    |    └── ...
|    |    ├── KDE.py
|    |    ├── L_rdc_qtpegasis.py
|    |    └── sampling_qtpegasis.py
|    ├── classgroup.py
|    ├── ContinuedFractions.py
|    ├── euclideanspace.py
|    ├── L_rdc_algorithm.py
|    ├── quaternion.py
|    ├── sampling_qlapoti.py
|    ├── StatisticalTesting.py
|    └── Visualization.py
├── Log_and_plots/
|    └── ...
├── probability/
|    ├── pdf_F_and_f.py
|    ├── Qt-pegasis_prob_G.py
|    ├── Qt-pegasis_prob.py
|    ├── w_k.py
|    └── ...
└── SubdominantEigenvalue/
    └── SubdominantEigenvalue.py
```

# Content
The most important functions and files are listed below:\
The folders **Log_and_plots** contain all the created plots as well as some logged data and results.

**constant.py:**\
```compute_eigenvectors(m, plot=False)```
computes left and right eigenvector (normalized) of the matrix $T_m$.

```constant_c(m):```
computes $f_4[u]$, the constant used in the probability estimate using the uniform transfer operator.

```dual_f_star(x, m, ev_left=None):```
computes $f^*_4[x]$, the dual eigenfunctional of a coefficient vector $x$.

```eigenvector_f(w, m, ev_right=None):```
computes the value $f_4(w)$ for a number $w$.

```Constant_C_4(m, grid_res=300, plot=False):```
computes the constant $C_4$ used in the probability estimate using the general transfer operator.

**DominantEigenvalue.py**\
```G_eval(f, s, M, t):```
computes $\mathcal{G}_s[f](t)$.

```G_eval_div_f(f,s,M,t):```
computes $\frac{\mathcal{G}_s[f](t)}{f(t)}$.

```M_coeff(i,j, s=4):```
computes the matrix coefficients $M_{i,j}$.

```matrix_T_m(m, s=4):```
creates the matrix $T_m$ out of the coefficients $M_{i,j}$.

```alpha_beta(f,s,M,m,left=0.0,right=1.0, N=2000):```
computes the $\alpha$ and $\beta$ to get a bound on the dominant eigenvalue $\lambda_4$.

```lambda_estimate(s,M,m):```
computes the approximate eigenvalue $\lambda_4$ using the $\alpha$ and $\beta$.

```p(m,s=4):```
computes the function $p_m$ which is the truncated version of the eigenfunction $f_4$.


## LatticeReduction/
**ContinuedFractions.py**\
```apply_word(word, z):```
computes $h_{\mathbf{m}}(z)$ for the word $\mathbf{m}$.

```radius_of_word(word):```
computes the radius of the disk $h_{\mathbf{m}}(\mathcal{D})$ corresponding to word $\mathbf{m}$.

```continued_fractions(m_max, k):```
computes all the disks $h(\mathcal{D})$ for a fixed upper bounds on the $m_i$ in the word  $\mathbf{m} = [m_1, \dots, m_k]$.

```interval(word):```
computes the interval which is the diameter of the disk $h(\mathcal{D})$ corresponding to word $\mathbf{m}$.

```length_of_intersection_line(radius, d):```
computes the length of the intersection of a disk with radius $radius$ and the line $l:y=d$.

```Total_Length_Lineintersection(d):```
computes the total length of the intersection of the line $l:y=d$ with all possible domains $\mathcal{D}_k$.

**L_rdc_algorithm.py**\
```lattice_reduction_2dim(Obj, v, u):```
is the Standard Gaussian Algorithm in 2 dimensions, which is adjusted to the norms corresponding to the object $Obj$. This is used because we have the special norm from the quadratic field $\mathbb{Q}[\sqrt{\Delta}]$.



**sampling_qlaoti.py** contains the sampling of the Qlapoti construction, which is not finished, since the thesis didn't need it.



**StatisticalTesting.py** uses the ```stats.probplot``` method to compare the predicted distributions of number of reduction steps with the expected one.



**Visualization.py** \
```class Disk:```
is the class for the disks $h(\mathcal{D})$ and contains the necessary methods.

```class Domain_D_k:```
is the class for the domains $\mathcal{D}_k$ and contains the necessary methods.

```line_through_domains(d):```
returns all the domains $\mathcal{D}_k$ that the line $l:y=d$ intersects.

```line_intervals_domains(d,store_intervals=True):```
computes all the intervals $\mathcal{D}_k \cap l$ for every $k$.

```compute_length_list_intervals(li):```
computes the lengths for the $k$ using the list of intervals from above.

### QtPegasis/

**KDE.py**\
```K(x, H_det, H_inv):```
is the kernel function for the covariance $H$.

```pdf(x, H_det, H_inv, samples):```
is the KDE probability density function.

**sampling_qtpegasis.py**\
```random_class_group_element(D, multi=1):```
samples a random ideal class group element of the quadratic field $\mathbb{Q}[\sqrt{D}]$.

```histogram_of_lattice_reduction(D, num_samples, log=False, multi=10):```
applies the lattice reduction algorithm to the sampled lattices and counts the number of reducton steps. (Creates histogram out of it).

```plot_z_in_disk(D, num_samples=1000, multi=1, log=False): ```
plots the complex representation $z$ of the sampled lattices in the fundamental disk $\mathcal{D}$.

## SubdominantEigenvalue/

**SubdominantEigenvalue.py**\
```trace_4_G_sq(M)```
computes the trace of the uniform transfer operator squared $\mathcal{G}_4^2$.

```subdominant_eigenvalue_upperbound(lam=get_lambda(), M=1000):```
uses the trace value to compute the upper bound on the subdominant eigenvalue $\mu_4$.

## probability

**Qt-pegasis_prob.py**\
```normalization_constant_N(pdf, samples=100, plot=False):```
computes the normalization factor $\frac{1}{Z_d}$, which ensures that the probability density function is actually a probablity denstiy function.

```f_to_vector_any_function(function, m, samples=100,A=None, f_vals=None):```
computes the coefficient vector of a function with respect to the basis $\{(x-1/2)^i\}$.

```spectral_gap(m):```
computes the spectral gap used in the general transfer operator estimate. $1/lam_4  \max(lam_6, mu_4)$.

```Norm_F(F, pdf, samples=1000):```
computes the norm $\|F\|$ of the holomorphic extension $F$ of the probability density function $\rho$.

```probability_estimate(k, m, pdf, grid_res=100, C_4=None, plot=False):```
computes the leading term of the probablity estimate using the general transfer operator.

```f_dual_of_pdf(m, pdf):```
computes $f^*_4[\rho]$.

```prob(values):```
computes both the leading term and the error term of the asymptotical probablity estimate.


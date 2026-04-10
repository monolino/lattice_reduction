from sage.all import *

'''
Hurwitz zeta function: hurwitz_zeta(s,w)
Binomial: binomial(n,k)
'''

#variables
m, s, t, x = var('m s t x')
CF = ComplexField(1000) #100 bits precision
a = QQ(1)/2 #1/2 as element in Q


#transfer operator for uniform input density
def G_s(f,s):
    return sum( (1/(m+t)**s) * f(1/(m+t)), m, 1, infinity )

#transfer operator G_s with an upper sum bound M (not infinity as before)
def G_eval(f, s, M, t):
    return sum(
        (1/(m+t)**s) * f(1/(m+t)),
        m, 1, M
    )

'''
@param f: function f(t) which returns a value
@param s: the s in G_s
@param M: upper bound of the sum of G_s
@param t: number in [0,1] to evaluate

@output : evaluation of G_s(f)/f at t
'''
def G_eval_div_f(f,s,M,t):
    return G_eval(f, s, M, t)/f(t)


def f(t):
    return 1-87/50*(t-a)+ 391/200*(t-a)**2 - 1687/1000*(x-a)**3

def find_M(f):
    M = 1
    while abs(G_eval_div_f(f,4,M,0.5) - G_eval_div_f(f,4,1000,0.5)) > 0.000000000000000000000001:
        M += 1
    return M

def M_coeff(i,j):
    value = CF(0)
    for l in range(0,j+1):
        value += binomial(j,l)* binomial(i+l+3,i) * (-a)**(j-l) * hurwitz_zeta(4+l+i,a+1)
    value *= (-1)**i
    return value.n(100)

def matrix_T_m(m):
    return matrix(CF, m, m,lambda i, j: M_coeff(i, j))

def p_first(m, a=CC(1)/2):
    A = matrix_T_m(m)
    char = A.charpoly()
    #print('characteristic polynomial: ',char)
    eig_values = char.roots(ring=CF) #100 bit precision |lambda' - lambda| <= 2**(-100). list [(eigenvalue, multiplicity)]
    #print('eigenvalues of T_m', eig_values)

    lambda_dom = max([r[0] for r in eig_values], key=lambda x: abs(x))
    #print('dominant eigenvalue of T_m: ', lambda_dom)
    p = vector(CF, [1] * m) #initial vector [1,1,1,1,1,...]
    p /= p.norm() #normalize
    I = identity_matrix(CF, m) #identity matrix mxm
    '''
    Inverse-iteration:
    (A-lambda I)p_k+1 = p_k
    '''
    for k in range(30):   
        p = (A - lambda_dom*I).solve_right(p)
        p /= p.norm()
    #print('Eigenvector: ', p)
    def eig_function(x):
        return sum(p[i] * (x-a)**i for i in range(m))

    return eig_function

def alpha_beta(s,M,m):
    eigen_function = p(m)
    alpha = G_eval(eigen_function, s, M, 0)/eigen_function(0)
    beta = G_eval(eigen_function, s, M, 1)/eigen_function(1)
    return (alpha, beta)

def lambda_estimate(s,M,m):
    '''
    with noticing that Gf/f has a plateau over [0.1,0.8] i.e is nearly constant there and then take midpoint
    '''
    if m == 4:
        t_left, t_right = (0.1,0.8)
    if m == 8:
        t_left, t_right = (0.1,0.8)
    if m == 16:
        t_left, t_right = (0.1,0.1)
    eigen_function = p(m)
    mid = (t_left + t_right)/2
    #alpha, beta = (G_eval(eigen_function, s, M, 0.1)/eigen_function(0.1),G_eval(eigen_function, s, M, 0.8)/eigen_function(0.8))
    #error = beta-alpha
    return G_eval(eigen_function, s, M, mid)/eigen_function(mid)

def p(m):
    A = matrix_T_m(m)
    data = A.eigenvectors_right() #[eigval, [tuple], mult]
    dom = max(data, key=lambda data: abs(data[0])) #[dominanteigval, [dominant eigenvector], multiplicity]
    lambda_dom = dom[0]
    eigvec = dom[1][0]
    mult = dom[2]
    #print("lambda: ", lambda_dom, "\neigenvector: ", eigvec, "\nmult: ", mult)
    def eig_function(x):
        return sum(eigvec[i] * (x-a)**i for i in range(m))
    return eig_function
    

        


if __name__ == "__main__" :
    their_eigfunc = lambda t: 1-87/50*(t-a)+ 391/200*(t-a)**2 - 1687/1000*(t-a)**3

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

    print(alpha_beta(4,100000,64))
    print('lambda: ', lambda_estimate(4,1000,64))
    #print('alpha, beta = ', alpha_beta(4,1000,16))
    #print(matrix_T_m(16))

    
    #print(ef_test(16))
    #print(same_eigenvectors(ef_test(16),p(16)))
    #print(p(16)(1/2),"\n",ef_test(16)(1/2))
    #print(alpha_beta(4,1000,64))
    



    




from sage.all import *

'''
Hurwitz zeta function: hurwitz_zeta(s,w)
Binomial: binomial(n,k)
'''

#variables
m, s, t, x = var('m s t x')
precision = 1000
CF = ComplexField(precision) #100 bits precision
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
    t = CF(t)
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
    return value.n(precision)

def matrix_T_m(m):
    return matrix(CF, m, m,lambda i, j: M_coeff(i, j))


def G_div_f_numeric(f,s,M,t):
    t = CF(t)
    G_f = CF(0)
    for m in range(1, M + 1):
        G_f += (1/(m+t)**s) * f(1/(m+t))
    f = f(t)
    return G_f/f

def alpha_beta(f,s,M,m,left=0.0,right=1.0, N=2000):

    values= []
    for i in range(N+1):
        t = left + (right-left)*i/N #grid points on intervall
        values.append(real_part(G_div_f_numeric(f,s,M,t)))
    return min(values), max(values)

def lambda_estimate(s,M,m):
    '''
    with noticing that Gf/f has a plateau over [0.1,0.8] i.e is nearly constant there
    '''
    p = p(m)
    alpha,beta = alpha_beta(p,s,M,m,0.1,0.8)
    alpha = str(alpha)
    beta = str(beta)
    lambda_m = ""
    for x, y in zip(a, b):
        if x == y:
            prefix += x
        else:
            break
    return lambda_m.rstrip('.')

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
    
    m=32
    p = p(m)
    alpha, beta = alpha_beta(p,4,1000,m)
    print(alpha.n(100), " , ", beta.n(100))
    #print('lambda: ', lambda_estimate(4,1000,64))
    #print('alpha, beta = ', alpha_beta(4,1000,16))

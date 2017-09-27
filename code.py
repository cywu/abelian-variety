# This should be run by SageMath.
# For cusp resolution
def fracExpansion(w):
    #returns fractional expansion in the form of van der Geer p38 eq (3)
    #w is the w0 in loc. cit. satisfying w>1>w'>0
    l=[]
    
    while True:
        b=ceil(w)
        if (len(l)>0 and l[0]==b):
            break
        l.append(b)
        w=1/(b-w)
    
    return l

# Examples
fracExpansion((3+sqrt(5))/2)
fracExpansion((5+sqrt(13))/2)

# Some scripts for elliptic points
# Q(sqrt(5))
K.<a>=QuadraticField(5)
# unit
eps=(1+a)/2
epsConj = eps.galois_conjugate()
# the ideal $\mathfrak{p}
I=K.ideal(2)
# generator for isotropy group of each elliptic point for full level structure
gammas=[matrix(K,[[0,1],[-1,0]]),matrix(K,[[0,-epsConj],[-eps,0]]),matrix(K,[[0,1],[-1,1]]),matrix(K,[[0,-epsConj],[-eps,1]]),matrix(K,[[0,1],[-1,eps]]),matrix(K,[[1,-epsConj],[epsConj,-epsConj]])]
# representatives of residue field of $\mathfrak{p}
myRes=[0,1,eps,epsConj]
# now considering level structure
deltainf=matrix(K,[[0,1],[-1,0]])
deltas=[matrix(K,[[1,0],[x,1]]) for x in myRes];deltas.append(deltainf)
isotropyGen =[filter(lambda x: x[1][1][0] in I, [(delta,delta*gamma*delta^(-1)) for delta in deltas]) for gamma in gammas]

# Q(sqrt(13))
K.<a>=QuadraticField(13);
eps=(3+a)/2;
epsconj=eps.galois_conjugate();
O=K.maximal_order()
I=K.ideal(4+a)
delta=[Matrix(O,[[0,1],[-1,0]]),Matrix(O,[[1,0],[0,1]]),Matrix(O,[[1,0],[1,1]]),Matrix(O,[[1,0],[2,1]])]
# 2 torsion inequiv
# using Gundlach's method to find generators for isotropy group of inequivalent elliptic points for full level structure
A2=[Matrix(O,[[0,1],[-1,0]]),Matrix(O,[[0,-epsconj],[-eps,0]])]
# find generators for isotropy group of inequivalent elliptic points for $\Gamma_0(I)$ structure
B2=filter(lambda x: x[1][1][0] in I, [(g,g*gamma*g^(-1)) for g in delta for gamma in A2])
#for (3;1,1)-torsion inequiv
A311=[Matrix(O,[[ 0, 1],[-1, 1]]),Matrix(O,[[ 1/2*a + 3/2, 2],[-1/2*a - 5/2, -1/2*a - 1/2]])]
B311=filter(lambda x: x[1][1][0] in I, [(g,g*gamma*g^(-1)) for g in delta for gamma in A311])
#for (3;1,-1)-torsion inequiv
A31minus1=[Matrix(O,[[-1/2*a - 3/2, a + 1],[-1/2*a - 3/2, 1/2*a + 5/2]]),Matrix(O,[[-1, 1/2*a + 1/2],[-1/2*a + 1/2, 2]])]
B31minus1=filter(lambda x: x[1][1][0] in I, [(g,g*gamma*g^(-1)) for g in delta for gamma in A31minus1])

# Estimation of Chern numbers
def c1Estimate (D,n,p3=False,p2=False,usePrecisecp=False):
    # D, n as in paper
    # p3: p=3?
    # p2: p=2?
    # usePrecisecp: should we use a more precise formula to compute?
    expr1=n* D^1.5 / 180
    if (usePrecisecp):
        expr1+= precisecEstimate(D)
    else:
        expr1+= -D^0.5 / 2 * (3/2/pi^2 * log(D)^2 + 1.05 * log(D))
    expr2=0
    expr3=0
    if (p3):
        expr2=4/pi*sqrt(3*D)*log(3*D)
    else:
        expr2= 1/4/pi*sqrt(3*D)*log(3*D)
        
        
    if (p2):
        expr3= 3/pi*sqrt(4*D)*log(4*D) 
    return (expr1 - expr2 - expr3).n()

def c2Estimate(D,n):
    return n*D^1.5/360

# auxiliary functions
def sigma00(n):
    return len([d for d in range(1,n+1) if n%d == 0 if gcd(n/d,d) == 1])

def F(x):
    return sum([sigma00(n) * (2*sqrt(x^2 - n^2).n() /n + 1) for n in range(1,x+1)])

def logthing(x):
    return (6/pi^2*x*log(2*x)^2 + 2.1*x*log(2*x)).n()

def dd(x):
    return -F(x) +logthing(x)

def precisecEstimate (D):
    return -sum([sigma((D-x^2)/4,0) for x in range(-sqrt(D),sqrt(D)) if (D -x^2)%4==0])/2

def computen (D,p3=False,p2=False,usePrecisecp=False):
    denom1=D^1.5 / 180
    denom2=D^1.5 / 180 + D^1.5/360
    expr1= precisecEstimate(D) if usePrecisecp else -D^0.5 / 2 * (3/2/pi^2 * log(D)^2 + 1.05 * log(D))
    expr2= 4/pi*sqrt(3*D)*log(3*D) if p3 else 1/4/pi*sqrt(3*D)*log(3*D)
    expr3= 3/pi*sqrt(4*D)*log(4*D) if p2 else 0
    # first one ensures c_1^2 >0
    # second one ensures c_1^2 + c_2 >12       
        
    return  max((-(expr1 - expr2 - expr3)/denom1).n(), ((12-(expr1 - expr2 - expr3))/denom2).n())

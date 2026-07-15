msolve = True
try :
    Ideal([PolynomialRing(GF(3),'x',1).gens()[0]]).groebner_basis(algorithm="msolve")
except:
    print("msolve not found, using sagemath default GB solver.")
    msolve = False

### VOX Parameters


VOXparams =[ 
    (251,8,9,6,6,"I"),
    (251,4,5,13,6,"Ia"),
    (251,5,6,11,6,"Ib"),
    (251,6,7,9,6,"Ic"),
    (1021, 10, 11, 7, 7,"III"),
    (1021, 5, 6, 15, 7,"IIIa"),
    (1021, 6, 7, 13, 7,"IIIb"),
    (1021, 7, 8, 11, 7,"IIIc"),
    (4093, 12, 13, 8, 8,"V"),
    (4093, 6, 7, 17, 8,"Va"),
    (4093, 7, 8, 14, 8,"Vb"),
    (4093, 8, 9, 13, 8,"Vc")
]


### KeyGen
def complete_basis(B, E) :
    """
    This function takes as input a list of linearily independant vectors B in E and 
    completes them into a basis naively.
    If the list is empty, a random basis of E is returned
    """
    res = list(B) #Avoid modifying the input.
    if len(res) == 0 :
        b = E.random_element()
        while b.is_zero() :
            b = E.random_element()
        res.append(b)

    while len(res) != E.dimension() :
        b = E.random_element()
        if not(b in span(res)) :
            res.append(b)
    return res

def hpKeyGen(q,o,v,m,t ) :
    """ 
    m equations of a UOVhp(qovt) key, assuming odd characteristic.
    """
    n = o + v
    F = []
    
    #Hat plus equations
    for _ in range(t) : 
        f = matrix([[ GF(q).random_element() for _ in range(n)]for _ in range(n)] )
        F.append(f + f.transpose())
    #UOV equations
    for _ in range(m-t) :
        f = matrix([[ GF(q).random_element() for _ in range(n)]for _ in range(n)] )
        for i in range(o) :
            for j in range(o) :
                f[i,j] = 0
        F.append(f+f.transpose())

    A = matrix(complete_basis([], GF(q)^n))

    #Public key
    G = [A.transpose()*f*A for f in F]
    
    Sp = matrix([[GF(q).random_element() for _ in range(t) ]for _ in range(m)])
    
    G2 = [ G[i] + sum([Sp[i,j]*G[j] for j in range(t)]) for i in range(m)]
    
    return (A,Sp,F, G), G2







#### Cost estimation
def dreg_semi_reg(n,m) :
    """ 
    HS of a semi-regular quadratic system m in n variables.   
    """
    R.<t> = PowerSeriesRing(ZZ)
    h = (1-t^2)^(m)/(1-t)^(n) 

    L = list(h)
    for i in range(len(L)) :
        if L[i] <= 0 :
            return i
    return len(L)+1

def a_op(q) :
    """
    Approximated cost of arithmetic in Fq
    """
    return 2*log(q,2)^2 + log(q,2)


### IO ### 
def ToMSolve(F, finput="/tmp/in.ms"): #From msolve library interfaces
    """Convert a system of sage polynomials into a msolve input file.

    Inputs :
    F (list of polynomials): system of polynomial to solve
    finput (string): name of the msolve input file.

    """
    A = F[0].parent()
    assert all(A1 == A for A1 in map(parent,F)),\
            "The polynomials in the system must belong to the same polynomial ring."
    variables, char = A.variable_names(), A.characteristic()
    s = (", ".join(variables) + " \n"
            + str(char) + "\n")

    B = A.change_ring(order = 'degrevlex') 
    F2 = [ str(B(f)).replace(" ", "") for f in F ]
    if "0" in F2:
        F2.remove("0")
    s += ",\n".join(F2) + "\n"

    fd = open(finput, 'w')
    fd.write(s)
    fd.close()

if msolve :
    lim  = 60
else: 
    lim  = 50


print("We begin by estimating the attack cost for all relevant fields.")
exp_params = []
for q, o, v, c, t, name in VOXparams :
    print("### Parameter set " + name + " ###")
    for c1 in c.divisors()[::-1]:
        
        print('with l\' =',c1)
        c2 = c/c1
        N = (o+v)*c2 
        O = o*c2 
        print("nbr eqs",o*c,"dim",O,"kr",N-o*c, O-t)
        if N > o*c :  # Direct attack is not a key recovery attack in the chosen field.
            print("System is not overdetermined.")
            continue
        if O<=t :     # Direct attack fails in the chosen field.
            print("Variety is empty.")
            continue
        
        d = dreg_semi_reg(N-O+t,o*c)
        Cost =round(log(a_op(q**c2)*(o*d)*binomial(N-O+t+d-1, d)^2.81  , 2), 1)  
        
        
        print("Gate count:", Cost)
        print("Degree of regularity", d)

        if Cost < lim :
            exp_params.append([next_prime(q**c2), O, N-O, o*c,t, name, c1]) ### We use nextprime because msolve currently only supports prime fields.
            
        print()
        
### Experiments

print("We run experiments on instances with estimated costs below", lim)
if not msolve:
    print("Note that Sage's default GB algorithm is slower than msolve, which was used in experiments.")
import time
for q, o, v, m, t, name, c1 in exp_params:
    print("Parameter set", name, "with", "l' = "+str(c1))
    print("KeyGen")
    _, G = hpKeyGen(q,o,v,m,t)

    ## Modelisation
    R = PolynomialRing(GF(q), 'x', v+t)
    X = vector([1] + [0 for _ in range(o - t -1)] + list(R.gens()))

    eqs = [ X*g*X for g in G]
    I = Ideal(eqs)
    ToMSolve(eqs, "/tmp/VOX_"+name)
    print("Attack")
    t0 = time.time()
    if msolve :
        gb = I.groebner_basis(algorithm="msolve")
    else :
        gb = I.groebner_basis()
    dt = time.time() - t0
    print("GB computed in", round(dt,3),'s')
    print("Linear terms in the GB:")
    for i in gb:
        if i.degree()==1:
            print(i)  
    print(" ")
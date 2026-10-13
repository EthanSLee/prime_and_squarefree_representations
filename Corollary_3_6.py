# Preload bennett et al constants..

c_theta = {}

with open("cthetas.txt", "r") as file: # "cthetas.txt" is available at http://www.nt.math.ubc.ca/BeMaObRe/
    for line in file:
        if line.startswith("#"):
            continue
        columns = line.split()
        q = int(columns[0])
        c_theta[q] = float(columns[5])

am_comps = {
    1 : 0.328739,
    3 : 0.19272,
    5 : 0.24844,
    7 : 0.261257,
    15 : 0.143104,
    21 : 0.149224,
    35 : 0.16973,
    105 : 0.087477 
}

# Define functions..

def c(q,n=2*10**(12)):
    if q <= 100:
        return fdiv(1,q) + fdiv(1,420*log(n))
    if q <= 10**5:
        return fdiv(1,q) + fdiv(1,160*log(n))
    else:
        return 0.000041
def least_prime(m,n=2*10**(12)):
    primes = []
    indx = len(factorint(m))
    np = nextprime(max(factorint(m)))
    while indx < 4:
        primes.append(np)
        np = nextprime(np)
        indx += 1
    while am_comps[m] <= fdiv(log(2),n) + sum(c(p) for p in primes):
        primes = [nextprime(p) for p in primes]
    return primes

# Code to identify least conditions..

for m in [3,5,7,15,21,35,105]:
    print(f"If $m = {m}$, then \\eqref{{eqn:equivalent_problem}} holds when $q_{{j+1}} \\geq {least_prime(m)[0]}$.")

# Define more functions..

def CArtin():
    return 0.373955813619202288054728
def E(m,n):
    if m == 1:
        return float(fdiv(0.375,power(log(n),2)))
    if m > 1:
        return c_theta[m]
def one(m):
    if m == 1:
        return 1
    else:
        return 0
def inner_sum(d,k,n):
    s = 0
    for e in divisors(k):
        s += E(e*d**2,n)
    return s
def squarefree(a):
    for p, e in factorint(a).items():
        if e > 1:
            return False
    return True
def sum1_m(k,n=2*10**(12)):
    c = floor(10**(5/2)/sqrt(k))
    sum1 = 0
    d = 1
    while d <= c:
        if squarefree(d) == True and gcd(d,k) == 1:
            sum1 += inner_sum(d,k,n)
        d += 1
    return sum1
def sum2_m(k,n=2*10**(12),Z=10**5):
    c = floor(10**(5/2)/sqrt(k))
    sum2 = 0
    sum3 = 0
    d = int(c + 1)
    while d <= Z:
        if squarefree(d) == True and gcd(d,k) == 1:
            sum2 += fdiv(1,totient(d**2))
            sum3 += fdiv(1,totient(d**2)*(1 - fdiv(log(k*d**2),log(n))))
        d += 1
    return [sum2 + fdiv(4,Z), sum3]
def sum1_q(q,n=2*10**(12)):
    c1 = floor(sqrt(10**(5)/q))
    c2 = floor(10**(5/2)/q)
    sum1A, sum1B = 0, 0
    d = 1
    while d <= max(c1,c2):
        if squarefree(d) == True and gcd(d,q) == 1:
            if d <= c1:
                sum1A += E(q*d**2,n)
            if d <= c2:
                sum1B += E((q*d)**2,n)
        d += 1
    return [sum1A, sum1B]
def sum2_q(q,n=2*10**(12),Z=10**5):
    c1 = floor(10**(5/2)/q)
    c2 = floor(sqrt(10**(5)/q))
    sum2_A = 0
    sum2_B = 0
    sum3_A = 0
    sum3_B = 0
    d = 1
    while d <= Z:
        if squarefree(d) == True and gcd(d,q) == 1:
            if d >= int(c1 + 1):
                sum2_A += fdiv(1,totient(d**2))
                sum3_A += fdiv(1,totient(d**2))*power(1-fdiv(log(q*d**2),log(n)),-1)
            if d >= int(c2 + 1):
                sum2_B += fdiv(1,totient(d**2))
                sum3_B += fdiv(1,totient(d**2))*power(1-fdiv(2*log(q*d),log(n)),-1)
        d += 1
    sum2_A += fdiv(4,Z)
    sum2_B += fdiv(4,Z)
    return [sum2_A, sum2_B, sum3_A, sum3_B]

def CurlyE1(m,A,s1,s2,s3,n=2*10**(12),Z=10**5):
    prd1, prd2 = 1, 1
    for prm in factorint(m).keys():
        prd1 *= fdiv(prm-2,prm-1)
        prd2 *= fdiv(prm+1,prm)
    om = len(factorint(m))
    return fdiv(s1,log(n)) + fdiv(power(n,A-1)*log(m),sqrt(m)) + (s2 + 2*s3 + fdiv(8,Z*(1-2*A)))*prd1 + (fdiv(power(2,om)*power(n,A-1),sqrt(m)) + fdiv(prd2,fdiv(power(n,A),sqrt(m)) - 1) + fdiv(power(2,om),sqrt(n)))*log(n)
def CurlyE2(q,A2,s1A,s1B,s2A,s2B,s3A,s3B,n=2*10**(12),Z=10**5):
    c1 = floor(10**(5/2)/q)
    c2 = floor(sqrt(10**(5)/q))
    out = (q*power(n,-2*A2) + power(n,-A2) + fdiv(power(n,A2-1),q))*log(n)
    for i in range(1,3):
        if i == 1:
            out += fdiv(s1A,log(n)) + fdiv(s2A,totient(q)) + fdiv(2,totient(q))*(s3A + fdiv(4,Z*(1-2*A2)))
        if i == 2:
            out += fdiv(s1B,log(n)) + fdiv(s2B,totient(q**2)) + fdiv(2,totient(q**2))*(s3B + fdiv(4,Z*(1-2*A2)))
    return out

def LHS(m,k):
    prd1, sum1 = 1, 0
    for prm in factorint(k).keys():
        if m % prm == 0:
            prd1 *= fdiv(prm*(prm-2),prm*(prm-1)-1)
        else:
            sum1 += fdiv(prm-1,prm*(prm-1)-1)
    return (prd1 - sum1)*CArtin()
def RHS_As(m,As,k,s1,sum1_q_pairs,s2,sum2_q_pairs,s3,n=2*10**(12)):
    counter = 0
    out = fdiv(log(2),n) + CurlyE1(m,As[counter],s1,s2,s3,n)
    for q in factorint(k//m):
        counter += 1
        s1A, s1B = sum1_q_pairs[q]
        s2A, s2B, s3A, s3B = sum2_q_pairs[q]
        out += CurlyE2(q,As[counter],s1A,s1B,s2A,s2B,s3A,s3B,n)
    return float(out)

# More definitions..

from scipy.optimize import differential_evolution

def optimise_As(m,k,n=2*10**(12)):
    s1_m = sum1_m(m,n)
    s2_m, s3_m = sum2_m(m,n)
    sum1_q_pairs = {}
    sum2_q_pairs = {}
    for q in factorint(k//m):
        sum1_q_pairs[q] = sum1_q(q,n)
        sum2_q_pairs[q] = sum2_q(q,n)
    
    def objective(x):
        return RHS_As(m,x,k,s1_m,sum1_q_pairs,s2_m,sum2_q_pairs,s3_m,n)

    result = differential_evolution(objective, bounds=[(0.2, 0.495) for u in range(len(factorint(k//m))+1)], tol=1e-6, polish=True, seed=1)
    return result.x, RHS_As(m, result.x, k, s1_m, sum1_q_pairs, s2_m, sum2_q_pairs, s3_m, n)
  
def condition(m,k):
    prd1, sum1 = 1, 0
    for prm in factorint(k).keys():
        if m % prm == 0:
            prd1 *= fdiv(prm*(prm-2),prm*(prm-1)-1)
        else:
            sum1 += fdiv(prm-1,prm*(prm-1)-1)
    return prd1 > sum1
def optimise_A(m,k,n=2*10**12,Z=10**5):
    As_opt, maximum = optimise_As(m,k,n)
    return [round(float(A_opt),6) for A_opt in As_opt], float(LHS(m,k)), float(maximum), LHS(m,k) > maximum

# Code to produce Table 2..

cases = {
    1 : [3,[11,13]],
    2 : [5,[11]],
    3 : [7,[11]],
    4 : [15,[11,13]],
    5 : [21,[11,13]],
    6 : [35,[11]],
    7 : [105,[11,13]],
}

for c_no in range(1,8):
    M, options = cases[c_no]
    for o in options:
        k = M
        while len(factorint(k)) < 5:
            k *= o
            o = nextprime(o)
        oAs, left, right, condition_holds = optimise_A(M,k,2*10**12,10**5)
        A_optimal = oAs[0]
        Aqs_optimal = oAs[1:]
        print(f"${M}$ & ${k}$ & ${A_optimal}$ & {str(Aqs_optimal)[1:-2]} & ${left:.6f}$ & ${right:.6f}$ \\\\") 

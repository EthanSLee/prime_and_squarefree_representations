# Preload bennett et al constants..

c_theta = {}

with open("cthetas.txt", "r") as file: # "cthetas.txt" is available at http://www.nt.math.ubc.ca/BeMaObRe/
    for line in file:
        if line.startswith("#"):
            continue
        columns = line.split()
        q = int(columns[0])
        c_theta[q] = float(columns[5])

# Main code..

from sympy import nextprime, isprime, prevprime, factorint, primerange, floor, totient, gcd, divisors, mobius
import mpmath as mpm
from mpmath import mpf, sqrt, pi, log, exp, power, fdiv, ceil, floor, euler, sqrt
from numpy import linspace

mpm.mp.dps = 120

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
def sum1(k,n):
    c = floor(10**(5/2)/sqrt(k))
    sum1 = 0
    d = 1
    while d <= c:
        if squarefree(d) == True and gcd(d,k) == 1:
            sum1 += inner_sum(d,k,n)
        d += 1
    return sum1
def sum2(Z,k,n):
    c = floor(10**(5/2)/sqrt(k))
    sum2 = 0
    sum3 = 0
    d = int(c + 1)
    while d <= Z:
        if squarefree(d) == True and gcd(d,k) == 1:
            sum2 += fdiv(1,totient(d**2))
            sum3 += fdiv(1,totient(d**2)*(1-fdiv(log(k*d**2),log(n))))
        d += 1
    return sum2 + fdiv(4,Z), sum3 
def main_term(k,n):
    prd1 = 1
    for q in factorint(k):
        prd1 *= ((q-2)/(q-1))*(q**2 - q)/(q**2 - q - 1)
    return float(prd1*CArtin())
def lower_bound(A,k,n,sum1,sum2,sum3):
    Z = 10**5
    c = floor(10**(5/2)/sqrt(k))
    if k == 1:
        ok = 0
        prd1 = 1
        prd2 = 1
        prd3 = 1
    else:
        ok = 0
        prd1 = 1
        prd2 = 1
        prd3 = 1
        for q in factorint(k):
            prd1 *= (q**2 - 2*q)/(q**2 - q - 1)
            prd2 *= (q-2)/(q-1)
            prd3 *= (q+1)/q
            ok += 1
    lower_bound = prd1*CArtin() - fdiv(sum1,log(n)) - fdiv(power(n,A-1)*log(k),sqrt(k)) - (sum2 + 2*sum3 +  fdiv(8,Z*(1-2*A)))*prd2 - (fdiv(power(2,ok),sqrt(k))*power(n,A-1) + fdiv(prd3,fdiv(power(n,A),sqrt(k)) - 1) + fdiv(power(2,ok) + one(k),sqrt(n)))*log(n)
    return float(lower_bound)
def optimise_A(k,n):
    s1_comp = sum1(k,n)
    s2_comp, s3_comp = sum2(10**5,k,n)
    optimalA = 0.2
    for diffs in [0.01,0.001,0.0001,0.00001,0.000001]:
        while lower_bound(optimalA + diffs,k,n,s1_comp,s2_comp,s3_comp) > lower_bound(optimalA,k,n,s1_comp,s2_comp,s3_comp):
            optimalA += diffs
        optimalA -= diffs
    optimalA += 0.000001
    return optimalA, lower_bound(optimalA,k,n,s1_comp,s2_comp, s3_comp), s1_comp, s2_comp, s3_comp

for y in [1,3,5,7,3*5,3*7,5*7,3*5*7]:
    n0 = 2*10**(12)
    A_optimised, lb, s1_comp, s2_comp, s3_comp = optimise_A(y,n0)
    print(f"${y}$ & ${round(A_optimised,6)}$ & ${round(s1_comp + 0.00001,5)}$ & ${round(s2_comp + 0.00001,5)}$ & ${round(s3_comp + 0.00001,5)}$ & ${round(lb-0.000001,6)}$ \\\\")

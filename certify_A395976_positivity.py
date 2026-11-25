"""Finite rational certificates for the infinite positivity proof in the supplement.
Run with Python 3 and SymPy. All assertions use exact arithmetic.
"""
from fractions import Fraction as Q
import sympy as s
from math import factorial
from check_A395976_recurrence import calc

def base_certificate():
    d=calc(400)
    assert d[1]==Q(55,24) and d[4]==Q(-4579919,7962624)
    q=[None]+[d[j]*Q(4**j,factorial(j-1)) for j in range(1,401)]
    R=Q(4,5); U=Q(13,10); K=20
    assert all(x>0 for x in d[6:])
    assert q[K]>=1
    assert all(R*q[j-1]<=q[j]<=U*q[j-1] for j in range(K+1,401))
    assert all(abs(q[j])*R**(K-j)<=100*q[K] for j in range(1,K))
    print('Exact base: d_6,...,d_400 positive; ratio and early-index bounds verified.')

def polynomial_certificate():
    j,t=s.symbols('j t'); M=30; L=26; n0=400
    R=s.Rational(4,5); U=s.Rational(13,10)
    weights=[{0:1},{0:-3,2:-1,4:1},{0:3,2:-5,4:1,6:1},{0:-1,2:6,4:-13,6:12,8:-4}]
    zero=s.Poly(0,j,domain=s.QQ)
    coeff=[zero for k in range(M+1)]
    for r,w in enumerate(weights):
        log=[None]+[s.Poly(4*s.binomial(s.Rational(1,2),(k+1)//2)*(-r)**((k+1)//2) if k%2 else s.Rational(r**(k//2+1),(k//2)*(k//2+1))+j*s.Rational(r**(k//2),k),j,domain=s.QQ) for k in range(1,M+1)]
        E=[s.Poly(1,j,domain=s.QQ)]
        for k in range(1,M+1):
            E.append(sum((i*log[i]*E[k-i] for i in range(1,k+1)),zero).mul_ground(s.Rational(1,k)))
        for k in range(M+1): coeff[k]+=sum((v*E[k-i] for i,v in w.items() if i<=k),zero)
    assert all(p.is_zero for p in coeff[:4]); assert coeff[4].as_expr()==-4*j
    den=s.Poly((t+n0)*s.prod(t+n0-i for i in range(1,L+1)),t,domain=s.QQ)
    low=s.Poly(0,t,domain=s.QQ); high=low
    for l in range(1,L+1):
        p=s.Poly((4**l*coeff[l+4].as_expr()/4).subs(j,t+n0-l).expand(),t,domain=s.QQ)
        pos=all(x>=0 for x in p.all_coeffs()); neg=all(x<=0 for x in p.all_coeffs()); assert pos or neg
        fac=s.Poly(s.prod(t+n0-i for i in range(l+1,L+1)),t,domain=s.QQ)
        v=p*fac
        if pos: low+=v.mul_ground(U**(1-l)); high+=v.mul_ground(R**(1-l))
        else: low+=v.mul_ground(R**(1-l)); high+=v.mul_ground(U**(1-l))
    assert all(x>=0 for x in (low-den.mul_ground(s.Rational(81,100))).all_coeffs())
    assert all(x>=0 for x in (den.mul_ground(s.Rational(119,100))-high).all_coeffs())
    print('Exact polynomial certificate: 81/100 <= head <= 119/100 for all N>=400.')
    near=300*400*Q(3,8)**27/(1-Q(3,8))
    far=100*300*400*Q(3,5)**200/(1-Q(3,5))
    zero=20736*Q(9,80)**400
    assert near+far+zero<Q(1,1000)
    print('Exact rational tail bound: < 1/1000.')

if __name__=='__main__':
    polynomial_certificate()
    base_certificate()
    print('ALL FINITE CERTIFICATES PASSED.')

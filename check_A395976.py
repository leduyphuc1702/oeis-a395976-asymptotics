from fractions import Fraction as Q
from math import factorial,comb
import sympy as sp
import mpmath as mp
# Work with v=i*u and integrate via tilted Gaussian moments.
def add(a,b):
    out=a.copy()
    for k,v in b.items(): out[k]=out.get(k,0)+v
    return {k:v for k,v in out.items() if v}
def mul(a,b):
    out={}
    for k,v in a.items():
        for l,w in b.items(): out[k+l]=out.get(k+l,0)+v*w
    return out
def scale(a,s): return {k:v*s for k,v in a.items() if v*s}
def calc(N):
    b=[Q(1)]
    for j in range(1,N+1): b.append(b[-1]*Q((2*j-1)**2,16*j))
    B=[sum(b[j]*b[k-j] for j in range(k+1)) for k in range(N+1)]
    moments=[1,-2]
    for r in range(1,3*N): moments.append(-2*moments[-1]-r*moments[-2])
    L=[{}]+[{k+2:Q(1,factorial(k+2)),k+1:Q(4,2**(k+1)*factorial(k+1))} for k in range(1,N+1)]
    for m in range(1,(N+2)//4+1):
        k=4*m-2
        if k<=N: L[k][0]=Q(2*sp.bernoulli(2*m),2*m*(2*m-1))
    d=[Q(0) for _ in range(N+1)]
    for j in range(N+1):
        L1=add(L[1],{1:Q(-j-1,2)})
        exp=[{0:Q(1)}]
        d[j]+=B[j]
        for k in range(1,N-j+1):
            p={}
            for r in range(1,k+1): p=add(p,scale(mul(L1 if r==1 else L[r],exp[k-r]),Q(r,k)))
            exp.append(p)
            d[j+k]+=B[j]*sum(v*moments[i] for i,v in p.items())
    return d
if __name__=='__main__':
    import sys
    N=int(sys.argv[1]) if len(sys.argv)>1 else 20
    d=calc(N)
    expected=[Q(1),Q(55,24),Q(2581,1152),Q(255893,414720),Q(-4579919,7962624),Q(-1075995491,6688604160),Q(2191041054791,4815794995200),Q(4477083471901,23115815976960),Q(454915545479411,3170169048268800),Q(222043505853843095117,263631258054033408000),Q(2554196391424824098947,1610547321930095001600)]
    assert d[:min(N+1,len(expected))]==expected[:min(N+1,len(expected))]
    print('OEIS coefficients matched:',min(N+1,len(expected)))
    for n in [1,2,3,4,5,6,10,15,20,25,30,40,50]:
        if n<=N:
            q=d[n]*Q(4**n,factorial(n-1)); print(n,float(q)*float(mp.pi/(2*mp.e**2)), 'positive' if d[n]>0 else 'negative')

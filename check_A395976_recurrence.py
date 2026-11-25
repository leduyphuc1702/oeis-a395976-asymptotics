"""Exact independent recurrence from A336293; no floating point in sign tests."""
from fractions import Fraction as Q
from math import factorial, comb
import sys

def calc(N):
    M=N+4
    fact=[factorial(k) for k in range(M+1)]
    sc=[2**k*fact[k] for k in range(M+1)]
    choose=[[comb(k,i) for i in range(k+1)] for k in range(M+1)]
    # 2^k k! times coefficient of each log shift, excluding j-dependent part.
    # binomial(1/2,a) computed exactly.
    half=[Q(1)]
    for a in range(1,M+2): half.append(half[-1]*Q(3-2*a,2*a))
    weights={1:{0:-3,2:-1,4:1},2:{0:3,2:-5,4:1,6:1},3:{0:-1,2:6,4:-13,6:12,8:-4}}
    base={}; linear={}
    for r in (1,2,3):
        base[r]=[0]*(M+1); linear[r]=[0]*(M+1)
        for k in range(1,M+1):
            if k%2: v=4*half[(k+1)//2]*(-r)**((k+1)//2)*sc[k]
            else:
                a=k//2; v=Q(r**(a+1),a*(a+1))*sc[k]
                w=Q(r**a,2*a)*sc[k]; assert w.denominator==1; linear[r][k]=w.numerator
            assert v.denominator==1; base[r][k]=v.numerator
    acc=[Q(0) for _ in range(N+1)]; d=[Q(1)]
    for j in range(N+1):
        if j: d.append(acc[j]/(4*j))
        if j==N: break
        m=N-j+4
        coeff=[0]*(m+1)
        for r in (1,2,3):
            L=[a+j*b for a,b in zip(base[r][:m+1],linear[r][:m+1])]
            E=[1]
            for k in range(1,m+1):
                E.append(sum(choose[k-1][i-1]*L[i]*E[k-i] for i in range(1,k+1)))
            for k in range(5,m+1):
                coeff[k]+=sum(w*E[k-q]*sc[k]//sc[k-q] for q,w in weights[r].items() if q<=k)
        for k in range(5,m+1): acc[j+k-4]+=d[j]*Q(coeff[k],sc[k])
    return d

if __name__=='__main__':
    N=int(sys.argv[1]) if len(sys.argv)>1 else 100
    d=calc(N)
    expected=[Q(1),Q(55,24),Q(2581,1152),Q(255893,414720),Q(-4579919,7962624),Q(-1075995491,6688604160),Q(2191041054791,4815794995200),Q(4477083471901,23115815976960),Q(454915545479411,3170169048268800),Q(222043505853843095117,263631258054033408000),Q(2554196391424824098947,1610547321930095001600)]
    assert d[:min(N+1,len(expected))]==expected[:min(N+1,len(expected))]
    nonpositive=[j for j in range(6,N+1) if d[j]<=0]
    print('N =',N,'; exact OEIS matches =',min(N+1,len(expected)))
    print('Nonpositive indices in [6,N]:',nonpositive)
    for j in (6,10,20,50,100,150,200,300):
        if j<=N: print(j, 'normalized d_j:',float(d[j]*Q(4**j,factorial(j-1))))
    print('All rational signs checked without rounding.')

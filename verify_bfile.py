"""Compare both independent implementations with the OEIS b-files of A395976/A395977.

Usage: python3 verify_bfile.py [b395976.txt b395977.txt]
Without arguments the b-files are downloaded from oeis.org.
"""
import sys
import urllib.request
from fractions import Fraction as Q

import check_A395976 as saddle
import check_A395976_recurrence as recurrence


def read_bfile(source):
    if source.startswith('http'):
        req = urllib.request.Request(source, headers={'User-Agent': 'Mozilla/5.0'})
        text = urllib.request.urlopen(req).read().decode()
    else:
        text = open(source).read()
    rows = (line.split() for line in text.splitlines() if line.strip() and not line.startswith('#'))
    return {int(n): int(v) for n, v in rows}


if __name__ == '__main__':
    if len(sys.argv) == 3:
        num, den = read_bfile(sys.argv[1]), read_bfile(sys.argv[2])
    else:
        num = read_bfile('https://oeis.org/A395976/b395976.txt')
        den = read_bfile('https://oeis.org/A395977/b395977.txt')
    N = min(max(num), max(den))
    oeis = [Q(num[i], den[i]) for i in range(N + 1)]
    assert all(oeis[i].numerator == num[i] and oeis[i].denominator == den[i] for i in range(N + 1))
    print('b-file terms: 0 ..', N)
    print('recurrence implementation matches:', recurrence.calc(N) == oeis)
    print('saddle-point formula matches:     ', saddle.calc(N) == oeis)
    print('nonpositive numerators with n > 5:', [i for i in range(6, N + 1) if num[i] <= 0])

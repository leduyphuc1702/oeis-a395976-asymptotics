# Large-Order Asymptotics and Positivity of the Coefficients of OEIS A395976/A395977

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24521402.svg)](https://doi.org/10.5281/zenodo.24521402)

**Preprint (Zenodo DOI):** [https://doi.org/10.5281/zenodo.24521402](https://doi.org/10.5281/zenodo.24521402)

This repository contains the LaTeX manuscript and exact verification code for two results on the rational coefficients $d_j = \mathrm{A395976}(j)/\mathrm{A395977}(j)$ in the asymptotic expansion

$$a_m \sim \frac{m^m e^{4\sqrt m - m - 2}}{\sqrt{8\pi}} \sum_{j\ge0} d_j\, m^{-j/2}, \qquad a_m = \mathrm{A336293}(m) = \sum_{k=0}^m \binom mk^2 \binom{2k}k (m-k)!.$$

1. **Large-order asymptotics** (Vaclav Kotesovec's conjecture of May 21, 2026, in [A395976](https://oeis.org/A395976)):
   $$d_j \sim \frac{2e^2}{\pi}\frac{\Gamma(j)}{4^j} \sim \frac{j^{\,j-1/2}}{\sqrt\pi\, 2^{2j-3/2}\, e^{\,j-2}}.$$
2. **Positivity** (the question "For n > 5, are all a(n) positive?" in [A395976](https://oeis.org/A395976)): $d_j > 0$ for every $j \ge 6$. The proof is computer-assisted, with exact rational certificates and a uniform analytic tail bound.

## Contents

- `paper.tex`, `paper.pdf`: the manuscript.
- `check_A395976.py`: exact evaluation of the saddle-point formula for $d_j$ (SymPy, mpmath).
- `check_A395976_recurrence.py`: independent exact evaluation of the linear recurrence for $d_j$.
- `verify_bfile.py`: checks both implementations against all published terms of the OEIS b-files.
- `certify_A395976_positivity.py`: the finite certificates of the positivity proof (about 80 s).

## Reproducing

```bash
pip install sympy mpmath
python3 verify_bfile.py                 # both implementations vs OEIS b-files, terms 0..99
python3 certify_A395976_positivity.py   # prints ALL FINITE CERTIFICATES PASSED.
python3 check_A395976.py 35             # normalized coefficients approaching 1
```

## Status

Preprint, not yet peer-reviewed. Parts of the derivations and code were developed with AI assistance; every computational claim is backed by the exact scripts above.

## License

Manuscript text: CC BY 4.0. Code: MIT License.

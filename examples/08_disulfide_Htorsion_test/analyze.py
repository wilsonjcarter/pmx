#!/usr/bin/env python
"""Jarzynski (reverse) free-energy estimate from a set of reverse switching
works W_rev(B->A).

  <exp(-W_rev/kT)>_B = exp(-dG(B->A)/kT)      =>   dG(A->B) = -dG(B->A)

Also reports mean/std of the work (the <W> upper bound is exact only in the
limit of infinitely slow switching)."""
import sys, numpy as np

_trapz = getattr(np, "trapezoid", np.trapz)   # numpy>=2 renamed trapz
W = np.loadtxt(sys.argv[1] if len(sys.argv) > 1 else "works.txt", ndmin=1)
kT = 8.314462618e-3 * 300.0                       # kJ/mol at 300 K

dG_BA = -kT * np.log(np.mean(np.exp(-W / kT)))    # Jarzynski, B->A
dG_AB = -dG_BA

print(f"n switches           : {len(W)}")
print(f"<W_rev(B->A)>        : {W.mean():8.2f} +/- {W.std(ddof=1) if len(W)>1 else 0:.2f} kJ/mol")
print(f"Jarzynski dG(B->A)   : {dG_BA:8.2f} kJ/mol")
print(f"  => dG(A->B)        : {dG_AB:8.2f} kJ/mol = {dG_AB/4.184:.2f} kcal/mol")
print("\nNote: single-direction Jarzynski is biased for few/fast switches; for a")
print("converged dG run both directions and use BAR/Crooks (needs a well-defined")
print("state-A ensemble -> the split-atom/dummy hybrid, cf. the disulfide demo).")

"""Write the field at one nail texel (default: the centre, as in contact_spike.py) for dynamics.c.
Output (text): B0 (3, tesla), B90 (3), surface normal (3), kappa0 (1/(T^2 s)).
Usage: python center_field.py ../spike/fields.json field_center.txt [eta_PaS] [texel]"""
import json
import sys
import numpy as np

D = json.load(open(sys.argv[1]))
eta = float(sys.argv[3]) if len(sys.argv) > 3 else 300.0
cnt = D["count"]
t = int(sys.argv[4]) if len(sys.argv) > 4 else cnt // 2
B0 = np.array(D["B0"]).reshape(cnt, 3)[t]
B90 = np.array(D["B90"]).reshape(cnt, 3)[t]
N = np.array(D["N"]).reshape(cnt, 3)[t]
mu0, a, dchi = 4e-7 * np.pi, 6e-6, 3.0
kappa0 = 3 * np.pi * (a / 10) * dchi / (32 * mu0 * eta * a)
with open(sys.argv[2], "w") as f:
    f.write(" ".join("%.10e" % v for v in (*B0, *B90, *N, kappa0)) + "\n")
print("texel %d: |B0| = %.1f mT, |B90| = %.1f mT, kappa0 = %.0f" % (t, 1e3 * np.linalg.norm(B0), 1e3 * np.linalg.norm(B90), kappa0))

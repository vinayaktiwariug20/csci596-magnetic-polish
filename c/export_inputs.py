"""Write the inputs of platelet_spike.py (fields, per-flake kappa, initial normals) as one binary file for platelet.c.
The random draws repeat the spike's order (kappa first, then x0), so the C result can be compared texel by texel
with platelet_spike_maps.npz.  Usage: python export_inputs.py ../spike/fields.json inputs.bin [eta_PaS]"""
import json
import sys
import numpy as np

D = json.load(open(sys.argv[1]))
eta = float(sys.argv[3]) if len(sys.argv) > 3 else 300.0
cnt, NF = D["count"], 8
mu0, a, dchi = 4e-7 * np.pi, 6e-6, 3.0
h = a / 10.0
kappa0 = 3 * np.pi * h * dchi / (32 * mu0 * eta * a)
rng = np.random.default_rng(1)
kap = kappa0 * np.exp(0.55 * rng.standard_normal((cnt, NF)))
x0 = rng.standard_normal((cnt, NF, 3))
x0 /= np.linalg.norm(x0, axis=-1, keepdims=True)
with open(sys.argv[2], "wb") as f:
    np.array([cnt, NF], dtype=np.int32).tobytes()
    f.write(np.array([cnt, NF], dtype=np.int32).tobytes())
    for arr in (D["N"], D["B0"], D["B90"], kap.ravel(), x0.ravel()):
        f.write(np.asarray(arr, dtype=np.float64).tobytes())
print("wrote", sys.argv[2], "cnt =", cnt, "NF =", NF, "kappa0 = %.0f" % kappa0)

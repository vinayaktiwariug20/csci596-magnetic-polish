"""Reference for contact.c: runs overlap_counts() taken from spike/contact_spike.py (extracted with ast, so the
spike's own code is used) on random flakes, and writes the same flakes for contact.c.
Usage: python contact_ref.py N phi orient in.bin ref.bin     orient = random | aligned"""
import ast
import sys
import numpy as np

src = open("../spike/contact_spike.py").read()
fn = [n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef) and n.name == "overlap_counts"]
g = {"np": np, "a": 1.0, "h": 0.1}
exec(compile(ast.Module(body=fn, type_ignores=[]), "contact_spike.py", "exec"), g)

N, phi, orient = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
rng = np.random.default_rng(7)
L = (N * np.pi * 0.1 / phi) ** (1 / 3)
pos = rng.uniform(0, L, (N, 3))
if orient == "random":
    n = rng.standard_normal((N, 3))
else:   # mostly aligned: small spread about z, exercises the near-parallel slab branch
    n = np.array([0, 0, 1.0]) + 0.05 * rng.standard_normal((N, 3))
n /= np.linalg.norm(n, axis=1, keepdims=True)
with open(sys.argv[4], "wb") as f:
    f.write(np.int32(N).tobytes()); f.write(np.float64(L).tobytes())
    f.write(pos.astype(np.float64).tobytes()); f.write(n.astype(np.float64).tobytes())
ref = g["overlap_counts"](pos, n, L)
ref.astype(np.int32).tofile(sys.argv[5])
print("python: N = %d  phi = %.2f  %s  mean = %.4f   (random-orientation theory phi*pi*a/h = %.4f)" % (N, phi, orient, ref.mean(), phi * np.pi * 10))

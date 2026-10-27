"""Run soft.c at several nail texels and collect end flatness (no contact, walls only, and contact at given loadings).
Usage: python scan_line.py out.csv texel [texel ...]
Each texel: center_field.py writes its field, then soft.c runs ideal (no contact), free at 0.1% (walls only) and
free at each loading in LOADINGS, all in a 12a film with N flakes."""
import csv
import re
import subprocess
import sys

import numpy as np

N, FILM, SEED = 5000, 12, 1
LOADINGS = [0.03, 0.10]
EXE = "./soft"


def run(field, phi, mode, film):
    args = [EXE, field, str(phi), str(N), mode, str(SEED)] + ([str(film)] if film else [])
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
    line = [l for l in out.splitlines() if l.startswith(mode)][-1]
    get = lambda key: float(re.search(key + r" ([-0-9.]+)", line).group(1))
    prof = [l for l in out.splitlines() if "at end" in l]
    bands = [float(v) for v in re.findall(r"a ([0-9.]+) \(", prof[0])] if prof else []
    return get("hold1"), get("end"), bands, get("wall")


rows = []
for t in sys.argv[2:]:
    field = "field_%s.txt" % t
    subprocess.run([sys.executable, "center_field.py", "../spike/fields.json", field, "300", t],
                          capture_output=True, text=True, check=True)
    B = [float(v) for v in open(field).read().split()]
    b0, b90, nn = np.array(B[0:3]), np.array(B[3:6]), np.array(B[6:9])
    l = np.cross(b0, b90); l /= np.linalg.norm(l)
    lim = float(np.dot(l, nn / np.linalg.norm(nn)) ** 2)   # geometric limit flatness
    cases = [("none", 0.03, "ideal", 0), ("walls", 0.001, "free", FILM)] + [("%g%%" % (100 * p), p, "free", FILM) for p in LOADINGS]
    for name, phi, mode, film in cases:
        h1, end, bands, wall = run(field, phi, mode, film)
        rows.append(dict(texel=t, limit=round(lim, 4), case=name, hold1=h1, end=end,
                         nail_band=bands[0] if bands else "", mid_band=bands[len(bands) // 2] if bands else "", wall_s=wall))
        print("texel %s (limit %.3f) %-6s hold1 %.3f end %.3f  nail band %s  mid band %s  (%.0f s)" %
              (t, lim, name, h1, end, bands[0] if bands else "-", bands[len(bands) // 2] if bands else "-", wall), flush=True)

with open(sys.argv[1], "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)

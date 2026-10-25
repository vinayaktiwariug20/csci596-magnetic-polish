"""Compare out.bin from platelet.c with the maps the Python spike saved. Usage: python compare.py out.bin ../spike/platelet_spike_maps.npz"""
import sys
import numpy as np

out = np.fromfile(sys.argv[1], dtype=np.float64)
cnt = len(out) // 2
ref = np.load(sys.argv[2])
for name, c in (("platelet_2", out[:cnt]), ("platelet_3", out[cnt:])):
    r = ref[name]
    print("%-11s max|C-py| = %.2e   mean C = %.4f  mean py = %.4f" % (name, np.abs(c - r).max(), c.mean(), r.mean()))

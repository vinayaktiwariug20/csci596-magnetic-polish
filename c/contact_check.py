"""Per-flake comparison of contact.c output with the Python reference. Usage: python contact_check.py out.bin ref.bin"""
import sys
import numpy as np

a = np.fromfile(sys.argv[1], np.int32)
b = np.fromfile(sys.argv[2], np.int32)
print("   flakes with different counts: %d of %d" % (int((a != b).sum()), len(a)))

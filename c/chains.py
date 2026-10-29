"""Chains in a flake state written by soft.c: flakes linked when their centres are closer than RLINK and the line
between them is within ANGLE of the applied field B0. Prints chain-length statistics for each file.
Usage: python chains.py flakes.csv [more.csv ...]"""
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from read_flakes import load

RLINK = 2.2          # centre distance in flake radii a (touching edge to edge is 2)
ANGLE = 37.0         # max angle between the pair line and B0, deg (|cos| > 0.8)


def chains(head, d):
    box = head["box"]
    b = head["B0_T"] / np.linalg.norm(head["B0_T"])
    x = d[:, :3]
    # periodic in x and y; z is bounded by walls, so pad it to avoid wrapping
    tree = cKDTree(np.mod(x, [box[0], box[1], box[2] + 2 * RLINK]), boxsize=[box[0], box[1], box[2] + 2 * RLINK])
    pairs = tree.query_pairs(RLINK, output_type="ndarray")
    dv = x[pairs[:, 1]] - x[pairs[:, 0]]
    dv[:, :2] -= box[:2] * np.round(dv[:, :2] / box[:2])
    cosb = np.abs(dv @ b) / np.linalg.norm(dv, axis=1)
    p = pairs[cosb > np.cos(np.radians(ANGLE))]
    n = len(x)
    g = coo_matrix((np.ones(len(p)), (p[:, 0], p[:, 1])), shape=(n, n))
    _, lab = connected_components(g, directed=False)
    size = np.bincount(lab)
    per_flake = size[lab]
    return len(p) / n, per_flake.mean(), np.mean(per_flake >= 3), size.max()


if __name__ == "__main__":
    print("file                  links/flake  mean chain (per flake)  in chains >=3  longest")
    for path in sys.argv[1:]:
        head, d = load(path)
        lk, mean, f3, mx = chains(head, d)
        print("%-22s %8.3f %14.2f %16.3f %10d" % (path, lk, mean, f3, mx))

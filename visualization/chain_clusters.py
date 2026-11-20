"""Chains as connected clusters: static 45 s vs 0/90 switching 60 s.

Uses the same link rule as c/chains.py and c/chains_field.py (centres closer
than 2.2a, pair line within 37 deg of the field) and colours every connected
cluster of 3 or more flakes on its own; smaller ones are grey.  Links that
cross the periodic box edges are counted but not drawn.  Right panel: how
many flakes sit in clusters of each size.

Usage:
    python visualization/chain_clusters.py static.csv switching.csv --out chain_clusters.png
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "c"))
from read_flakes import load  # noqa: E402

RLINK, ANGLE = 2.2, 37.0


def clusters(head, d, field_key):
    box = head["box"]
    b = head[field_key] / np.linalg.norm(head[field_key])
    x = d[:, :3]
    span = [box[0], box[1], box[2] + 2 * RLINK]
    pairs = cKDTree(np.mod(x, span), boxsize=span).query_pairs(RLINK, output_type="ndarray")
    dv = x[pairs[:, 1]] - x[pairs[:, 0]]
    dv[:, :2] -= box[:2] * np.round(dv[:, :2] / box[:2])
    keep = np.abs(dv @ b) / np.linalg.norm(dv, axis=1) > np.cos(np.radians(ANGLE))
    linked = pairs[keep]
    n = len(x)
    g = coo_matrix((np.ones(len(linked)), (linked[:, 0], linked[:, 1])), shape=(n, n))
    _, labels = connected_components(g, directed=False)
    return labels, linked


def draw(ax, head, d, field_key, title, rng):
    labels, linked = clusters(head, d, field_key)
    sizes = np.bincount(labels)
    x = d[:, :2]
    box = head["box"][:2]
    big = np.flatnonzero(sizes >= 3)
    colour = np.full((len(x), 3), 0.82)
    tab = plt.get_cmap("tab20")
    palette = [tab(i)[:3] for i in range(20) if i not in (6, 7, 14, 15)]   # no reds, no greys
    order = big[np.argsort(-sizes[big])]
    size = np.full(len(x), 9.0)
    for rank, lab in enumerate(order):
        colour[labels == lab] = palette[rank % len(palette)] if rank else (0.8, 0.0, 0.0)
    size[labels == order[0]] = 22
    seg = np.stack([x[linked[:, 0]], x[linked[:, 1]]], 1)
    short = np.all(np.abs(seg[:, 1] - seg[:, 0]) < box / 2, axis=1)
    for (p, q) in seg[short]:
        ax.plot([p[0], q[0]], [p[1], q[1]], color="0.35", lw=0.5, zorder=1)
    ax.scatter(x[:, 0], x[:, 1], c=colour, s=size, edgecolors="none", zorder=2)
    biggest = sizes[order[0]]
    ax.set_title(f"{title}\nlargest cluster {biggest} flakes (large red dots), "
                 f"{len(big)} clusters of 3 or more", fontsize=10)
    ax.set_xlim(0, box[0]); ax.set_ylim(0, box[1])
    ax.set_aspect("equal")
    ax.set_xlabel("x / a"); ax.set_ylabel("y / a")
    fdir = head[field_key][:2] / np.linalg.norm(head[field_key][:2])
    c = np.array([4.5, 4.5])
    ax.annotate("", xy=c + 3 * fdir, xytext=c - 3 * fdir,
                arrowprops=dict(arrowstyle="<->", lw=1.8), zorder=3)
    ax.text(*(c + [1.2, 1.2]), field_key.replace("_T", ""), fontsize=10, weight="bold", zorder=3)
    return sizes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("static")
    ap.add_argument("switching")
    ap.add_argument("--out", default="chain_clusters.png")
    args = ap.parse_args()

    rng = np.random.default_rng(0)
    fig = plt.figure(figsize=(15, 5.4), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.8])
    hs, ds = load(args.static)
    hw, dw = load(args.switching)
    s_sizes = draw(fig.add_subplot(gs[0]), hs, ds, "B0_T", "Static field, 45 s (links along B0)", rng)
    w_sizes = draw(fig.add_subplot(gs[1]), hw, dw, "B90_T", "0°/90° switching, 60 s (links along B90)", rng)

    ax = fig.add_subplot(gs[2])
    bins = np.unique(np.round(np.logspace(0, np.log10(400), 14)).astype(int))
    for sizes, lab, c in [(s_sizes, "static 45 s", "C3"), (w_sizes, "switching 60 s", "C0")]:
        ax.hist(sizes, bins=bins, weights=sizes / sizes.sum(), histtype="step", lw=2, color=c, label=lab)
    ax.set_xscale("log")
    ax.set_xlabel("cluster size (flakes)")
    ax.set_ylabel("fraction of flakes in clusters of that size")
    ax.legend()
    ax.set_title("Where the flakes are", fontsize=10)
    fig.suptitle("Flake chains as connected clusters (3% loading, 2000 flakes, one seed; the box is periodic, so clusters wrap across its edges)", fontsize=13)
    fig.savefig(args.out, dpi=160)
    print(f"wrote {args.out}; largest static {s_sizes.max()}, switching {w_sizes.max()}")


if __name__ == "__main__":
    main()

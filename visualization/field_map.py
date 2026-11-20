"""Magnet field over the nail, and where it puts the cat-eye bead.

Panels 1-2: the two magnet positions (B0, B90) from spike/fields.json, seen
from above.  Colour is |B| in mT; streamlines follow the field component in
the nail surface; dashed contours give how steeply the field leaves the
surface (deg).
Panel 3: under 0/90 switching a flake can only satisfy both fields with its
normal along B0 x B90 (its plane contains both).  The colour is the angle of
that preferred normal from straight up; the bead (flakes facing a viewer and
light overhead) sits where it is small.

Usage: python visualization/field_map.py spike/fields.json --out field_map.png
"""

import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay


def unit(v):
    return v / np.maximum(np.linalg.norm(v, axis=-1, keepdims=True), 1e-30)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fields")
    ap.add_argument("--out", default="field_map.png")
    args = ap.parse_args()

    d = json.load(open(args.fields))
    n = d["count"]
    P = np.array(d["P"]).reshape(n, 3)
    N = unit(np.array(d["N"]).reshape(n, 3))
    B = {"B0": np.array(d["B0"]).reshape(n, 3), "B90": np.array(d["B90"]).reshape(n, 3)}

    tri = Delaunay(P[:, :2])
    gx = np.linspace(P[:, 0].min(), P[:, 0].max(), 160)
    gy = np.linspace(P[:, 1].min(), P[:, 1].max(), 220)
    GX, GY = np.meshgrid(gx, gy)
    pix = np.stack([GX.ravel(), GY.ravel()], 1)
    outside = (tri.find_simplex(pix) < 0).reshape(GY.shape)

    def grid(values):
        g = LinearNDInterpolator(tri, values)(pix).reshape(GY.shape + np.shape(values)[1:])
        g[outside] = np.nan
        return g

    fig, axes = plt.subplots(1, 3, figsize=(14, 5.6), constrained_layout=True)
    bmax = max(np.linalg.norm(b, axis=1).max() for b in B.values()) * 1e3

    for ax, (name, b) in zip(axes[:2], B.items()):
        mag = np.linalg.norm(b, axis=1) * 1e3
        normal = np.einsum("tk,tk->t", b, N)
        tangent = b - normal[:, None] * N                         # field in the nail surface
        out_deg = np.degrees(np.arcsin(np.clip(np.abs(normal) / np.linalg.norm(b, axis=1), 0, 1)))
        m = ax.pcolormesh(GX, GY, grid(mag), cmap="viridis", vmin=0, vmax=bmax, shading="gouraud")
        t = grid(tangent)
        ax.streamplot(gx, gy, t[..., 0], t[..., 1], color="white", density=0.9, linewidth=0.8,
                      arrowsize=0.9)
        cs = ax.contour(GX, GY, grid(out_deg), levels=[10, 20, 30, 40], colors="#ff6b4a",
                        linewidths=1.2, linestyles="--")
        ax.clabel(cs, fmt="%d°", fontsize=7)
        ax.set_title(f"{name} magnet position\n|B| (colour), surface field lines, "
                     "out-of-surface angle (dashed)", fontsize=10)
        fig.colorbar(m, ax=ax, shrink=0.7, label="|B| (mT)")

    pref = unit(np.cross(B["B0"], B["B90"]))
    ang = np.degrees(np.arccos(np.abs(pref[:, 2])))
    ax = axes[2]
    m = ax.pcolormesh(GX, GY, grid(ang), cmap="magma_r", vmin=0, vmax=90, shading="gouraud")
    cs = ax.contour(GX, GY, grid(ang), levels=[5, 10, 20, 45], colors="w", linewidths=0.8)
    ax.clabel(cs, fmt="%d°", fontsize=7)
    i = np.argmin(ang)
    ax.plot(*P[i, :2], "c*", ms=12, mec="k", label=f"bead centre (min {ang[i]:.1f}°)")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title("Switching: preferred flake normal B0×B90\nangle from straight up "
                 "(bead where small)", fontsize=10)
    fig.colorbar(m, ax=ax, shrink=0.7, label="angle (deg)")

    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
    fig.suptitle("Magnet field over the nail", fontsize=13)
    fig.text(0.5, -0.02, "Left two: colour |B|, white lines = field along the nail surface, red dashed = angle at which "
             "the field leaves the surface.  Right: the bead (flakes facing straight up) is where the tilt is small.",
             ha="center", fontsize=9, color="0.3")
    fig.savefig(args.out, dpi=170, bbox_inches="tight")
    print(f"wrote {args.out}; bead centre at ({P[i,0]:.2f}, {P[i,1]:.2f}) mm, "
          f"|B| range {np.linalg.norm(B['B0'],axis=1).min()*1e3:.0f}-{bmax:.0f} mT")


if __name__ == "__main__":
    main()

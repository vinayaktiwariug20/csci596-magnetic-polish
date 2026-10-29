"""Read a flake state written by soft.c (last argument out.csv) and make a quick-look figure.
Usage: python read_flakes.py flakes.csv [out.png]
load(path) returns (header dict, array with columns x, y, z, nx, ny, nz, kappa); lengths are in flake radii a.
Figure: left, flakes seen from above as their projected ellipses, coloured by tilt from the nail normal;
right, where the normals point (equal-area projection of the upper hemisphere, nail normal at the centre)."""
import sys

import numpy as np


def load(path):
    head, rows = {}, []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                words = line[1:].split()
                if words and words[0] in ("box", "B0_T", "B90_T", "nail_normal"):
                    head[words[0]] = np.array([float(w) for w in words[1:4]])
            elif line[0].isalpha():
                continue                      # column names
            else:
                rows.append([float(v) for v in line.split(",")])
    return head, np.array(rows)


def summary(head, d):
    nn = head.get("nail_normal", np.array([0, 0, 1.0]))
    c = d[:, 3:6] @ nn
    print("%d flakes, box %s a, flatness <(n.N)^2> %.3f, median tilt %.1f deg"
          % (len(d), np.round(head.get("box", [0, 0, 0]), 2), np.mean(c ** 2), np.degrees(np.median(np.arccos(np.minimum(1, np.abs(c)))))))


def figure(head, d, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import EllipseCollection

    nn = head.get("nail_normal", np.array([0, 0, 1.0]))
    n = d[:, 3:6] * np.sign(d[:, 3:6] @ nn)[:, None]          # normals on the upper hemisphere
    tilt = np.degrees(np.arccos(np.clip(n @ nn, -1, 1)))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5))
    # a disc of radius 1 with normal n projects onto the xy plane as an ellipse with axes 2 and 2|nz|,
    # the short axis along the in-plane part of n
    ang = np.degrees(np.arctan2(n[:, 1], n[:, 0]))
    ec = EllipseCollection(2 * np.abs(n[:, 2]) + 1e-3, 2 * np.ones(len(n)), ang, units="xy", offsets=d[:, :2],
                           offset_transform=a1.transData, array=tilt, cmap="viridis", clim=(0, 90), alpha=0.6)
    a1.add_collection(ec)
    box = head.get("box", d[:, :3].max(0))
    a1.set_xlim(0, box[0]); a1.set_ylim(0, box[1]); a1.set_aspect("equal")
    a1.set_xlabel("x (a)"); a1.set_ylabel("y (a)"); a1.set_title("flakes from above (colour: tilt from nail normal)")
    fig.colorbar(ec, ax=a1, label="tilt (deg)")
    # equal-area (Lambert) projection about the nail normal
    e1 = np.cross(nn, [1, 0, 0] if abs(nn[0]) < 0.9 else [0, 1, 0]); e1 /= np.linalg.norm(e1); e2 = np.cross(nn, e1)
    r = np.sqrt(2 * (1 - np.clip(n @ nn, -1, 1)))
    phi = np.arctan2(n @ e2, n @ e1)
    a2.scatter(r * np.cos(phi), r * np.sin(phi), s=2, alpha=0.3)
    a2.add_patch(plt.Circle((0, 0), np.sqrt(2), fill=False, color="grey"))
    a2.set_xlim(-1.5, 1.5); a2.set_ylim(-1.5, 1.5); a2.set_aspect("equal"); a2.axis("off")
    a2.set_title("flake normals (centre = along nail normal, rim = edge-on)")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    print("wrote", out)


if __name__ == "__main__":
    head, d = load(sys.argv[1])
    summary(head, d)
    figure(head, d, sys.argv[2] if len(sys.argv) > 2 else sys.argv[1].rsplit(".", 1)[0] + ".png")

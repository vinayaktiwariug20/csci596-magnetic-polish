"""Cat-eye nail under a moving light: static field vs 0/90 switching.

Extends render_brdf.py. The flake normals come from the same independent-
platelet dynamics (a = 6 um, h/a = 0.1, dchi = 3, eta = 300,
seed 1, log-normal rate spread 0.55, RK4 dt = 0.05 s, 32 flakes per point;
brushed start (10 deg rms) by default, --start random for the isotropic one;
static 45 s, switching 6 x (5 s B0 + 5 s B90)), so no flake-flake forces.

What is new: instead of a fixed half-vector (0,0,1), each frame places a
point light and the eye at finite distance above the curved nail, and every
nail point gets its own half-vector H = unit(unit(L-P) + unit(E-P)).  The
flakes reflect where their normals point along H, so as the light moves the
bright bead slides across the nail, as a real cat-eye does when you tilt
your hand.  The nail also gets a dark diffuse base coat and a glossy top-coat
highlight from its surface normal, for reference: the top-coat highlight
moves like a mirror, the bead moves according to the flake orientations.

Usage:
    python visualization/light_sweep.py spike/fields.json --out light_sweep.mp4
Needs numpy, scipy, matplotlib and ffmpeg on PATH.
"""

import argparse
import os
import shutil
import subprocess
import tempfile

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay


def unit(v):
    return v / np.maximum(np.linalg.norm(v, axis=-1, keepdims=True), 1e-30)


def evolve_normals(b0, b90, flakes, nail_normal=None, start_deg=10.0, eta=300.0):
    """Same dynamics and settings as render_brdf.py.

    nail_normal given: brushed start (normal along it, rms tilt start_deg),
    as in render_brdf.py and soft.c; None: random start.
    """
    count = len(b0)
    mu0 = 4e-7 * np.pi
    a = 6e-6
    h = a / 10.0
    dchi = 3.0
    kappa0 = 3 * np.pi * h * dchi / (32 * mu0 * eta * a)
    rng = np.random.default_rng(1)
    kappa = kappa0 * np.exp(0.55 * rng.standard_normal((count, flakes)))

    def rhs(x, field):
        bx = np.einsum("tfk,tk->tf", x, field)
        return -kappa[..., None] * bx[..., None] * (field[:, None, :] - bx[..., None] * x)

    def evolve(x, field, duration, dt=0.05):
        for _ in range(int(round(duration / dt))):
            k1 = rhs(x, field)
            k2 = rhs(unit(x + 0.5 * dt * k1), field)
            k3 = rhs(unit(x + 0.5 * dt * k2), field)
            k4 = rhs(unit(x + dt * k3), field)
            x = unit(x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4))
        return x

    noise = rng.standard_normal((count, flakes, 3))
    if nail_normal is None:
        x0 = unit(noise)
    else:
        x0 = unit(nail_normal[:, None, :] + np.radians(start_deg) / np.sqrt(2) * noise)
    static = evolve(x0.copy(), b0, 45.0)
    switching = x0.copy()
    for _ in range(6):
        switching = evolve(switching, b0, 5.0)
        switching = evolve(switching, b90, 5.0)
    return static, switching


def light_path(frames):
    """Tilt the light along the nail's long axis (y), then across it (x)."""
    half = frames // 2
    t = np.linspace(0, 2 * np.pi, half, endpoint=False)
    along = np.stack([np.zeros(half), 32 * np.sin(t)], 1)
    t = np.linspace(0, 2 * np.pi, frames - half, endpoint=False)
    across = np.stack([32 * np.sin(t), np.zeros(frames - half)], 1)
    return np.concatenate([along, across])          # (frames, 2) tilt in deg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fields")
    ap.add_argument("--flakes", type=int, default=32)
    ap.add_argument("--width-deg", type=float, default=8.0)
    ap.add_argument("--start", choices=["brushed", "random"], default="brushed")
    ap.add_argument("--frames", type=int, default=240)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--distance", type=float, default=250.0, help="light and eye distance, mm")
    ap.add_argument("--cache", default="light_sweep_normals.npz")
    ap.add_argument("--out", default="light_sweep.mp4")
    ap.add_argument("--still", default="light_sweep_still.png")
    args = ap.parse_args()

    d = json.load(open(args.fields))
    n = d["count"]
    P = np.array(d["P"]).reshape(n, 3)              # nail surface points, mm
    N = unit(np.array(d["N"]).reshape(n, 3))        # nail surface normals
    b0 = np.array(d["B0"]).reshape(n, 3)
    b90 = np.array(d["B90"]).reshape(n, 3)

    if os.path.exists(args.cache):
        c = np.load(args.cache)
        static, switching = c["static"], c["switching"]
    else:
        print("evolving flake normals (a few minutes, cached afterwards)...")
        static, switching = evolve_normals(b0, b90, args.flakes,
                                           N if args.start == "brushed" else None)
        np.savez_compressed(args.cache, static=static, switching=switching)

    # Raster (top-down view of the nail).  Interpolate per-point colours.
    tri = Delaunay(P[:, :2])
    pad = 0.6
    nx = 360
    x = np.linspace(P[:, 0].min() - pad, P[:, 0].max() + pad, nx)
    ny = int(nx * (np.ptp(P[:, 1]) + 2 * pad) / (np.ptp(P[:, 0]) + 2 * pad))
    y = np.linspace(P[:, 1].min() - pad, P[:, 1].max() + pad, ny)
    X, Y = np.meshgrid(x, y)
    pix = np.stack([X.ravel(), Y.ravel()], 1)
    inside = (tri.find_simplex(pix) >= 0).reshape(ny, nx)

    width = np.radians(args.width_deg)
    eye = np.array([0.0, 0.0, args.distance])
    base_col = np.array([0.20, 0.035, 0.11])        # dark plum polish
    flake_col = np.array([0.95, 0.92, 1.00])        # silver flakes
    gain = 2.2

    def shade(normals, tilt):
        ax_, ay_ = np.radians(tilt)
        ldir = unit(np.array([np.sin(ax_), np.sin(ay_), 1.0]))
        L = args.distance * ldir
        to_l = unit(L - P)
        to_e = unit(eye - P)
        H = unit(to_l + to_e)                       # per-point half-vector
        cosang = np.abs(np.einsum("tfk,tk->tf", normals, H))
        flake = np.exp(-np.arccos(np.clip(cosang, 0, 1)) / width).mean(1)
        diffuse = 0.25 + 0.75 * np.clip(np.einsum("tk,tk->t", N, to_l), 0, 1)
        coat = 0.9 * np.clip(np.einsum("tk,tk->t", N, H), 0, 1) ** 1500
        rgb = base_col[None] * diffuse[:, None] +gain * flake[:, None] * flake_col[None] + coat[:, None]
        rgb = 1 - np.exp(-1.6 * rgb)                # soft tone map
        img = LinearNDInterpolator(tri, rgb)(pix).reshape(ny, nx, 3)
        out = np.zeros((ny, nx, 4))
        out[..., :3] = np.nan_to_num(img)
        out[..., 3] = inside
        return out, flake

    tilts = light_path(args.frames)
    bg = "#0d0b10"
    fig, axes = plt.subplots(1, 2, figsize=(8, 6.4), facecolor=bg)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.86, bottom=0.12, wspace=0.05)
    ims = []
    for ax, title in zip(axes, ["Static field, 45 s", "0°/90° switching, 60 s"]):
        ax.set_facecolor(bg)
        ims.append(ax.imshow(np.zeros((ny, nx, 4)), origin="lower",
                             extent=[x[0], x[-1], y[0], y[-1]], interpolation="bilinear"))
        ax.set_title(title, color="white", fontsize=12)
        ax.set_axis_off()
    fig.suptitle("Cat-eye nail under a moving light", color="white", fontsize=14)
    fig.text(0.5, 0.025,
             f"Reflection proxy: independent flakes, {args.start} start, no flake-flake forces, {args.width_deg:g}° lobe; "
             "light and eye 25 cm above the nail",
             color="#a59fb0", ha="center", fontsize=8)
    label = fig.text(0.5, 0.075, "", color="white", ha="center", fontsize=10)

    tmp = tempfile.mkdtemp(prefix="sweep_")
    for i, tilt in enumerate(tilts):
        for im, normals in zip(ims, [static, switching]):
            im.set_data(shade(normals, tilt)[0])
        axis = "along nail" if tilt[0] == 0 else "across nail"
        label.set_text(f"light tilt {tilt[0] + tilt[1]:+5.1f}° {axis}")
        fig.savefig(os.path.join(tmp, f"f{i:04d}.png"), dpi=150, facecolor=bg)
        if i == args.frames // 8:                   # light tilted ~ +32 deg along the nail
            fig.savefig(args.still, dpi=150, facecolor=bg)
        if i % 40 == 0:
            print(f"frame {i}/{args.frames}")

    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(args.fps),
                    "-i", os.path.join(tmp, "f%04d.png"),
                    "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-c:v", "libx264",
                    "-pix_fmt", "yuv420p", "-crf", "18", args.out], check=True)
    shutil.rmtree(tmp)
    print(f"wrote {args.out} and {args.still}")


if __name__ == "__main__":
    main()

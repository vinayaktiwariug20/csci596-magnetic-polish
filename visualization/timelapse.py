"""Time-lapse of the bead forming: static field vs 0/90 switching, 0-60 s.

Same independent-platelet dynamics as render_brdf.py (a = 6 um, h/a = 0.1,
dchi = 3, eta = 300, seed 1, log-normal rate spread 0.55,
RK4 dt = 0.05 s), no flake-flake forces.  Both protocols run for 60 s here so
the clocks match (static: B0 throughout; switching: 5 s B0, 5 s B90, ...).
Flakes start brushed (10 deg rms about the nail normal) unless --start random.

Top row: where the flake normals point at the nail centre (2000 flakes,
equal-area projection, nail normal at the centre).  Bottom row: the whole
nail lit from overhead (32 flakes per point, light and eye 25 cm above).

Usage: python visualization/timelapse.py spike/fields.json --out timelapse.mp4
"""

import argparse
import json
import os
import shutil
import subprocess
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay


def unit(v):
    return v / np.maximum(np.linalg.norm(v, axis=-1, keepdims=True), 1e-30)


def kappa0(eta=300.0):
    mu0, a = 4e-7 * np.pi, 6e-6
    return 3 * np.pi * (a / 10) * 3.0 / (32 * mu0 * eta * a)


def run(field_at, x, kappa, total, snap_every, dt=0.05):
    """Integrate dn/dt with RK4; field_at(t) gives the field per point."""
    def rhs(x, f):
        bx = np.einsum("tfk,tk->tf", x, f)
        return -kappa[..., None] * bx[..., None] * (f[:, None, :] - bx[..., None] * x)

    snaps = [x.astype(np.float16)]
    steps, every = int(round(total / dt)), int(round(snap_every / dt))
    for s in range(steps):
        f = field_at(s * dt)
        k1 = rhs(x, f)
        k2 = rhs(unit(x + 0.5 * dt * k1), f)
        k3 = rhs(unit(x + 0.5 * dt * k2), f)
        k4 = rhs(unit(x + dt * k3), f)
        x = unit(x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4))
        if (s + 1) % every == 0:
            snaps.append(x.astype(np.float16))
    return np.array(snaps)


def lambert(n):
    n = n * np.sign(n[..., 2:3] + 1e-12)               # upper hemisphere
    r = np.sqrt(2 / (1 + n[..., 2]))
    return n[..., 0] * r, n[..., 1] * r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fields")
    ap.add_argument("--total", type=float, default=60.0)
    ap.add_argument("--snap", type=float, default=0.25)
    ap.add_argument("--start", choices=["brushed", "random"], default="brushed")
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--cache", default="timelapse_snaps.npz")
    ap.add_argument("--out", default="timelapse.mp4")
    args = ap.parse_args()

    d = json.load(open(args.fields))
    n = d["count"]
    P = np.array(d["P"]).reshape(n, 3)
    N = unit(np.array(d["N"]).reshape(n, 3))
    b0 = np.array(d["B0"]).reshape(n, 3)
    b90 = np.array(d["B90"]).reshape(n, 3)
    centre = np.argmin(np.linalg.norm(P[:, :2], axis=1))

    static_f = lambda t: b0
    switch_f = lambda t: b0 if int(t // 5.0) % 2 == 0 else b90

    if os.path.exists(args.cache):
        c = np.load(args.cache)
        nail = {"static": c["nail_static"], "switching": c["nail_switching"]}
        pt = {"static": c["pt_static"], "switching": c["pt_switching"]}
    else:
        rng = np.random.default_rng(1)
        k0 = kappa0()
        kn = k0 * np.exp(0.55 * rng.standard_normal((n, 32)))
        x0 = rng.standard_normal((n, 32, 3))
        kp = k0 * np.exp(0.55 * rng.standard_normal((1, 2000)))
        p0 = rng.standard_normal((1, 2000, 3))
        if args.start == "random":
            x0, p0 = unit(x0), unit(p0)
        else:                       # brushed, 10 deg rms, as in soft.c
            s = np.radians(10.0) / np.sqrt(2)
            x0 = unit(N[:, None, :] + s * x0)
            p0 = unit(N[centre][None, None, :] + s * p0)
        nail, pt = {}, {}
        for name, f in [("static", static_f), ("switching", switch_f)]:
            print(f"evolving {name} ...")
            nail[name] = run(f, x0.copy(), kn, args.total, args.snap)
            pt[name] = run(lambda t: f(t)[centre:centre + 1], p0.copy(), kp, args.total, args.snap)[:, 0]
        np.savez_compressed(args.cache, nail_static=nail["static"], nail_switching=nail["switching"],
                            pt_static=pt["static"], pt_switching=pt["switching"])

    # Nail shading (overhead light and eye), as in light_sweep.py at zero tilt.
    dist = 250.0
    eye = np.array([0, 0, dist])
    H = unit(unit(eye - P) * 2)
    tri = Delaunay(P[:, :2])
    nx = 260
    gx = np.linspace(P[:, 0].min() - 0.5, P[:, 0].max() + 0.5, nx)
    ny = int(nx * np.ptp(P[:, 1] + 0) / np.ptp(P[:, 0]))
    gy = np.linspace(P[:, 1].min() - 0.5, P[:, 1].max() + 0.5, ny)
    GX, GY = np.meshgrid(gx, gy)
    pix = np.stack([GX.ravel(), GY.ravel()], 1)
    inside = (tri.find_simplex(pix) >= 0).reshape(ny, nx)
    base = np.array([0.20, 0.035, 0.11])
    diffuse = 0.25 + 0.75 * np.clip(np.einsum("tk,tk->t", N, unit(eye - P)), 0, 1)
    width = np.radians(8.0)

    def nail_image(normals):
        cosang = np.abs(np.einsum("tfk,tk->tf", normals, H))
        flake = np.exp(-np.arccos(np.clip(cosang, 0, 1)) / width).mean(1)
        rgb = base[None] * diffuse[:, None] + 2.2 * flake[:, None] * np.array([0.95, 0.92, 1.0])
        rgb = 1 - np.exp(-1.6 * rgb)
        img = np.zeros((ny, nx, 4))
        img[..., :3] = np.nan_to_num(LinearNDInterpolator(tri, rgb)(pix).reshape(ny, nx, 3))
        img[..., 3] = inside
        return img

    bg = "#0d0b10"
    fig = plt.figure(figsize=(8, 9.4), facecolor=bg)
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.35], hspace=0.12, wspace=0.05,
                          left=0.03, right=0.97, top=0.88, bottom=0.07)
    titles = {"static": "Static field (B0)", "switching": "0°/90° switching"}
    dots, ims, phase = {}, {}, {}
    for j, name in enumerate(["static", "switching"]):
        ax = fig.add_subplot(gs[0, j], facecolor=bg)
        th = np.linspace(0, 2 * np.pi, 200)
        ax.plot(np.sqrt(2) * np.cos(th), np.sqrt(2) * np.sin(th), color="0.5", lw=1)
        for deg in (10, 30):
            r = 2 * np.sin(np.radians(deg) / 2)
            ax.plot(r * np.cos(th), r * np.sin(th), color="0.35", lw=0.6, ls=":")
            ax.text(r * 0.72, r * 0.72, f"{deg}°", color="0.55", fontsize=7)
        dots[name] = ax.scatter([], [], s=2, color="#e8dcff", alpha=0.6)
        ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal"); ax.set_axis_off()
        ax.set_title(titles[name], color="white", fontsize=12)
        phase[name] = ax.text(0, -1.62, "", color="#ffcf5c", ha="center", fontsize=10)
        ax2 = fig.add_subplot(gs[1, j], facecolor=bg)
        ims[name] = ax2.imshow(np.zeros((ny, nx, 4)), origin="lower",
                               extent=[gx[0], gx[-1], gy[0], gy[-1]], interpolation="bilinear")
        ax2.set_axis_off()
    fig.text(0.5, 0.955, "The bead forming", color="white", ha="center", fontsize=15)
    clock = fig.text(0.5, 0.915, "", color="white", ha="center", fontsize=12)
    fig.text(0.5, 0.505, "top: flake normals at the nail centre (equal-area, nail normal at centre)   "
             "bottom: the whole nail lit from overhead", color="#a59fb0", ha="center", fontsize=8)
    fig.text(0.5, 0.02, f"Independent flakes, {args.start} start, no flake-flake forces, 8° reflection lobe",
             color="#a59fb0", ha="center", fontsize=8)

    tmp = tempfile.mkdtemp(prefix="tl_")
    frames = len(pt["static"])
    for i in range(frames):
        t = i * args.snap
        clock.set_text(f"t = {t:4.1f} s")
        for name in ["static", "switching"]:
            u, v = lambert(pt[name][i].astype(np.float32))
            dots[name].set_offsets(np.c_[u, v])
            ims[name].set_data(nail_image(nail[name][i].astype(np.float32)))
            on = "B0" if name == "static" or int(min(t, args.total - 1e-9) // 5) % 2 == 0 else "B90"
            phase[name].set_text(f"field: {on}")
        fig.savefig(os.path.join(tmp, f"f{i:04d}.png"), dpi=110, facecolor=bg)
        if i % 40 == 0:
            print(f"frame {i}/{frames}")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(args.fps),
                    "-i", os.path.join(tmp, "f%04d.png"), "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", args.out], check=True)
    shutil.rmtree(tmp)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()

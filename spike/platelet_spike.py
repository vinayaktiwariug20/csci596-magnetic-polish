"""
Feasibility spike: does a PLATELET flake (easy plane, two axes) make a compact bright patch under
alternating perpendicular cat-eye holds, where the current ROD model can only print a bar?

Inputs: fields.json, exported from the sandbox (src/core field solver): the field B (tesla) at every nail
texel for the cat-eye bar at 0 deg and rotated 90 deg, plus the surface normals.

Model (all constants are order-of-magnitude estimates, not calibrated):
  platelet with normal n, easy plane, anisotropic susceptibility dchi:
      U = (V dchi / 2 mu0) (B.n)^2   -> torque makes n perpendicular to B (plane contains B)
      viscous rotation of a disc about a diameter: zeta = (32/3) eta a^3
      dn/dt = -kappa (B.n) (B - (B.n) n),   kappa = 3 pi h dchi / (32 mu0 eta a)
  rod with axis t (the current sandbox model): dt/dt = +kappa (B.t)(B - (B.t) t)
Each flake has its own kappa (log-normal spread, as in the sandbox's kSpread).

Usage: python platelet_spike.py [fields.json] [eta_PaS]
"""
import json
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

path = sys.argv[1] if len(sys.argv) > 1 else "fields.json"
eta = float(sys.argv[2]) if len(sys.argv) > 2 else 300.0       # Pa.s (gel, from the sandbox gel test)

D = json.load(open(path))
nu, nv, cnt = D["nu"], D["nv"], D["count"]
N = np.array(D["N"]).reshape(cnt, 3)
B0 = np.array(D["B0"]).reshape(cnt, 3)
B90 = np.array(D["B90"]).reshape(cnt, 3)

# --- constants (estimates) ---
mu0 = 4e-7 * np.pi
a = 6e-6            # flake radius (sandbox default flakeRadius)
h = a / 10.0        # thickness (aspect ratio 10: assumption)
dchi = 3.0          # sandbox 'chi'
kappa0 = 3 * np.pi * h * dchi / (32 * mu0 * eta * a)      # 1/(T^2 s)
spread = 0.55       # sandbox kSpread
NF = 8              # flakes per texel (raise for smoother maps; cost is linear)
rng = np.random.default_rng(1)
kap = kappa0 * np.exp(spread * rng.standard_normal((cnt, NF)))


def unit(v):
    return v / np.maximum(np.linalg.norm(v, axis=-1, keepdims=True), 1e-30)


def rhs(x, B, sign):
    bx = np.einsum("tfk,tk->tf", x, B)                     # B . x
    return sign * kap[..., None] * bx[..., None] * (B[:, None, :] - bx[..., None] * x)


def evolve(x, B, T, sign, dt=0.05):
    for _ in range(int(round(T / dt))):
        k1 = rhs(x, B, sign)
        k2 = rhs(unit(x + 0.5 * dt * k1), B, sign)
        k3 = rhs(unit(x + 0.5 * dt * k2), B, sign)
        k4 = rhs(unit(x + dt * k3), B, sign)
        x = unit(x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4))
    return x


def protocol(x0, sign, cycles, hold=5.0):
    """cycles of (hold at B0, hold at B90). Returns state after first hold and after all cycles."""
    x = evolve(x0.copy(), B0, hold, sign)
    first = x.copy()
    x = evolve(x, B90, hold, sign)
    for _ in range(cycles - 1):
        x = evolve(x, B0, hold, sign)
        x = evolve(x, B90, hold, sign)
    return first, x


x0 = unit(rng.standard_normal((cnt, NF, 3)))
CYC = 6
plate_first, plate_end = protocol(x0, -1.0, CYC)     # sign -1: n goes perpendicular to B
rod_first, rod_end = protocol(x0, +1.0, CYC)         # sign +1: axis goes parallel to B


def flat_platelet(x):          # 1 = face parallel to the nail surface (mirrors)
    return np.mean(np.einsum("tfk,tk->tf", x, N) ** 2, axis=1)


def flat_rod(x):               # 1 = rod lying along the surface
    return np.mean(1 - np.einsum("tfk,tk->tf", x, N) ** 2, axis=1)


# geometric limit: alternating projections converge to b0 x b90
lim = unit(np.cross(unit(B0), unit(B90)))
flat_limit = np.einsum("tk,tk->t", lim, N) ** 2

maps = {
    "rod, after hold 1 (single cat eye)": flat_rod(rod_first),
    "rod, after 6 switching cycles": flat_rod(rod_end),
    "platelet, after hold 1": flat_platelet(plate_first),
    "platelet, after 6 switching cycles": flat_platelet(plate_end),
    "geometric limit  n ~ B0 x B90": flat_limit,
}


def shape(m, thr=0.9):
    """area fraction above thr and aspect ratio of that region (nail-local coordinates)."""
    P = np.array(D["P"]).reshape(cnt, 3)[:, :2]
    s = m > thr
    if s.sum() < 5:
        return float(s.mean()), None
    c = np.cov((P[s] - P[s].mean(0)).T)
    w = np.sort(np.linalg.eigvalsh(c))
    return float(s.mean()), float(np.sqrt(w[1] / max(w[0], 1e-12)))


print("eta = %.0f Pa.s, kappa0 = %.0f 1/(T^2 s)  (k at 50 mT = %.2f /s)" % (eta, kappa0, kappa0 * 0.05**2))
print("%-40s %10s %8s %8s" % ("map", "area>0.9", "aspect", "mean"))
for k, m in maps.items():
    af, asp = shape(m)
    print("%-40s %10.3f %8s %8.3f" % (k, af, "-" if asp is None else "%.2f" % asp, m.mean()))

P2 = np.array(D["P"]).reshape(cnt, 3)[:, :2]
np.savez("platelet_spike_maps.npz", P=P2, **{k.split(",")[0].replace(" ", "_") + "_" + str(i): m for i, (k, m) in enumerate(maps.items())})
fig, ax = plt.subplots(1, 5, figsize=(17, 4.6))
for a_, (k, m) in zip(ax, maps.items()):
    sc = a_.scatter(P2[:, 0], P2[:, 1], c=m, s=14, vmin=0, vmax=1, cmap="magma", marker="s")
    a_.set_title(k, fontsize=9)
    a_.set_aspect("equal"); a_.set_xticks([]); a_.set_yticks([])
fig.colorbar(sc, ax=ax, shrink=0.8, label="mirror-like fraction (1 = flat)")
fig.suptitle("Bar held 5 s at 0 deg, then 5 s at 90 deg, x6 (gel %.0f Pa.s): rod vs platelet" % eta)
fig.savefig("platelet_spike.png", dpi=150, bbox_inches="tight")
print("saved platelet_spike.png")

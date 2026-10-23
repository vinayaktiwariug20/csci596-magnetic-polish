"""
Contact spike: does flake-flake CONTACT matter, and how hard is it for oriented discs?

Flakes are thin discs (radius a, thickness h = a/10, as in platelet_spike.py). Positions are uniform random
(ideal-gas placement, i.e. no flake-flake forces yet), so the count below is the expected number of OVERLAPPING
NEIGHBOURS per flake: a number near or above 1 means contact cannot be ignored.

Orientations come from the platelet model: n(t) under the cat-eye field at the nail centre, with the sandbox's
alignment-rate spread, starting random, held on B0 then switched to B90. We count overlaps at rest (random),
after the first hold (aligned, but rates differ) and DURING the switch (flakes rotating at different rates
sweep through each other).

Disc-disc overlap: thin discs overlap iff their chords on the line where the two planes meet overlap.
Near-parallel discs fall back to a slab test (|dz| < h and in-plane distance < 2a).

Usage: python contact_spike.py [fields.json] [eta_PaS]
"""
import json
import sys
import numpy as np

path = sys.argv[1] if len(sys.argv) > 1 else "fields.json"
eta = float(sys.argv[2]) if len(sys.argv) > 2 else 300.0
D = json.load(open(path))
cnt = D["count"]
B0 = np.array(D["B0"]).reshape(cnt, 3)[cnt // 2]
B90 = np.array(D["B90"]).reshape(cnt, 3)[cnt // 2]

mu0 = 4e-7 * np.pi
a = 1.0                       # lengths in units of the flake radius
h = 0.1                       # thickness / radius
kappa0 = 3 * np.pi * (0.1 * 6e-6) * 3.0 / (32 * mu0 * eta * 6e-6)
rng = np.random.default_rng(2)


def unit(v):
    return v / np.maximum(np.linalg.norm(v, axis=-1, keepdims=True), 1e-30)


def rhs(x, B, kap):
    bx = x @ B
    return -kap[:, None] * bx[:, None] * (B[None, :] - bx[:, None] * x)


def evolve(x, B, kap, T, dt=0.05):
    for _ in range(int(round(T / dt))):
        k1 = rhs(x, B, kap)
        k2 = rhs(unit(x + 0.5 * dt * k1), B, kap)
        k3 = rhs(unit(x + 0.5 * dt * k2), B, kap)
        k4 = rhs(unit(x + dt * k3), B, kap)
        x = unit(x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4))
    return x


def overlap_counts(pos, n, L):
    """Number of overlapping neighbours per flake (periodic box of side L)."""
    N = len(pos)
    cnts = np.zeros(N)
    for i in range(N - 1):
        j = np.arange(i + 1, N)
        d = pos[j] - pos[i]
        d -= L * np.round(d / L)
        near = np.einsum("jk,jk->j", d, d) < (2 * a) ** 2
        if not near.any():
            continue
        j, d = j[near], d[near]
        ni, nj = n[i], n[j]
        c = np.cross(ni, nj)
        cn = np.linalg.norm(c, axis=1)
        hit = np.zeros(len(j), bool)
        # general case: chords on the line of intersection of the two planes
        g = cn > 0.05
        if g.any():
            dd = c[g] / cn[g, None]
            M = np.stack([np.broadcast_to(ni, (g.sum(), 3)), nj[g], dd], axis=1)   # rows: n_i, n_j, dir
            rhs_ = np.stack([np.zeros(g.sum()), np.einsum("jk,jk->j", nj[g], d[g]), np.zeros(g.sum())], axis=1)
            p0 = np.linalg.solve(M, rhs_[..., None])[..., 0]       # point on the line, relative to disc i centre
            ti = -np.einsum("jk,jk->j", p0, dd)                     # centre of i projected on the line
            di2 = np.einsum("jk,jk->j", p0 + ti[:, None] * dd, p0 + ti[:, None] * dd)
            tj = np.einsum("jk,jk->j", d[g] - p0, dd)               # centre of j projected on the line
            pj = p0 + tj[:, None] * dd - d[g]
            dj2 = np.einsum("jk,jk->j", pj, pj)
            ok = (di2 < a * a) & (dj2 < a * a)
            hi = np.sqrt(np.maximum(a * a - di2, 0))
            hj = np.sqrt(np.maximum(a * a - dj2, 0))
            ok &= np.abs(ti - tj) < (hi + hj)
            hit[g] = ok
        # near-parallel: slab test
        p = ~g
        if p.any():
            dz = np.abs(np.einsum("jk,k->j", d[p], ni))
            inpl = np.linalg.norm(d[p] - np.einsum("jk,k->j", d[p], ni)[:, None] * ni, axis=1)
            hit[p] = (dz < h) & (inpl < 2 * a)
        cnts[i] += hit.sum()
        cnts[j[hit]] += 1
    return cnts


# stages of the protocol
N = 2500
kap = kappa0 * np.exp(0.55 * rng.standard_normal(N))
x_rand = unit(rng.standard_normal((N, 3)))
print("kappa0 = %.0f 1/(T^2 s); eta = %.0f Pa.s" % (kappa0, eta))
stages = {"random (at rest)": x_rand}
x = x_rand.copy()
x = evolve(x, B0, kap, 2.0);  stages["mid-rotation, hold 1 (2 s)"] = x.copy()
x = evolve(x, B0, kap, 8.0);  stages["after hold 1 (10 s, aligned)"] = x.copy()
x = evolve(x, B90, kap, 1.5); stages["mid-rotation, switch to 90 deg (1.5 s)"] = x.copy()
x = evolve(x, B90, kap, 8.5); stages["after hold 2 (10 s, aligned)"] = x.copy()

print("\nmean overlapping neighbours per flake (ideal-gas positions), aspect ratio 10")
phis = [0.01, 0.03, 0.05, 0.10]
print("%-42s" % "stage" + "".join("phi=%-6.2f" % p for p in phis))
res = {}
for name, xs in stages.items():
    row = []
    for phi in phis:
        # N flakes of volume pi a^2 h in a box: phi = N pi a^2 h / L^3
        L = (N * np.pi * a * a * h / phi) ** (1 / 3)
        pos = rng.uniform(0, L, (N, 3))
        row.append(overlap_counts(pos, xs, L).mean())
    res[name] = row
    print("%-42s" % name + "".join("%-10.2f" % v for v in row))

# how much of the aligned state's contact is the in-plane/parallel kind
print("\nspread of normals after hold 1: mean angle to mean normal = %.1f deg" %
      np.degrees(np.arccos(np.clip(np.abs(stages["after hold 1 (10 s, aligned)"] @ unit(stages["after hold 1 (10 s, aligned)"].mean(0))), 0, 1)).mean()))
json.dump({"phis": phis, "res": res}, open("contact_spike_result.json", "w"), indent=1)

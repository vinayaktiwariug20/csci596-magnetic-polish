"""How much do the flakes' own magnetic fields matter? Two estimates, no dynamics.

1. Film (mean field). The magnetised film is a thin sheet, so its field inside cuts the component of B along the
   nail normal N: B_in = B_t + B_N / (1 + chi_NN) N, with chi_NN = phi [chi_in (1 - (n.N)^2) + chi_perp (n.N)^2]
   (linear flakes; n solved self-consistently as the geometric limit of the reduced fields). Gives the shifted
   end flatness at each scan texel.
2. Neighbours (local field). Flakes at random positions in a cube, all at the limit orientation, each a disc of
   radius a a ring of magnetic charge on its rim. The field of the other flakes (cutoff RC) is averaged over each disc; its
   component along n tilts the flake by about b_n / |B| (energy (V dchi / 2 mu0) (B.n)^2). Reports rms b_n / |B|.
   The mean over a cutoff sphere is ~0 by symmetry; the film part is item 1.
Usage: python dipole_estimate.py
"""
import numpy as np

A, H = 1.0, 0.1                      # disc radius and thickness (units of a)
V = np.pi * A * A * H
TEXELS = [3152, 3148, 3144, 3140, 3136, 3132, 3128]
PHIS = [0.01, 0.03, 0.10]
CHIS = [(4.0, 1.0), (8.0, 1.0)]      # (chi_in, chi_perp): model value dchi = 3, and a stronger bracket


def unit(v):
    return v / np.linalg.norm(v)


def limit(b0, b90, nn, phi, chi_in, chi_perp):
    n = unit(np.cross(b0, b90))
    for _ in range(50):
        c2 = np.dot(n, nn) ** 2
        s = phi * (chi_in * (1 - c2) + chi_perp * c2)
        red = lambda b: b - (s / (1 + s)) * np.dot(b, nn) * nn
        n = unit(np.cross(red(b0), red(b90)))
    return np.dot(n, nn) ** 2


print("1. Film field: end flatness (n.N)^2 at each texel")
for chi_in, chi_perp in CHIS:
    print("   chi_in %.0f, chi_perp %.0f" % (chi_in, chi_perp))
    print("   texel  " + "  ".join("phi %4.0f%%" % (100 * p) for p in [0.0] + PHIS))
    for t in TEXELS:
        f = np.loadtxt("field_%d.txt" % t)
        b0, b90, nn = f[0:3], f[3:6], unit(f[6:9])
        row = [limit(b0, b90, nn, p, chi_in, chi_perp) for p in [0.0] + PHIS]
        print("   %5d  " % t + "  ".join("%9.3f" % r for r in row))

# 2. local field, centre texel: B along x (hold at 0 deg), n = z. Non-overlapping parallel discs (rejection).
# A uniformly in-plane magnetised disc acts as a magnetic charge sigma = M cos(theta) on its rim (height H),
# so each source disc is P rim charges; the field is averaged over target points spread over the target disc.
rng = np.random.default_rng(1)
P, KT, RC, NF, NE = 48, 61, 6.0, 4000, 600
th = 2 * np.pi * (np.arange(P) + 0.5) / P
rim = np.stack([A * np.cos(th), A * np.sin(th), np.zeros(P)], 1)
qrim = np.cos(th) * H * A * (2 * np.pi / P)                  # charge per rim element, per unit M
u = rng.random((KT, 2)); tgt = np.stack([A * np.sqrt(u[:, 0]) * np.cos(2 * np.pi * u[:, 1]),
                                          A * np.sqrt(u[:, 0]) * np.sin(2 * np.pi * u[:, 1]), np.zeros(KT)], 1)


def place(phi):
    L = (NF * V / phi) ** (1 / 3)
    pos = np.zeros((NF, 3)); k = 0; tries = 0
    while k < NF:
        p = rng.random(3) * L; tries += 1
        d = pos[:k] - p; d -= L * np.round(d / L)
        if not np.any((np.hypot(d[:, 0], d[:, 1]) < 2 * A) & (np.abs(d[:, 2]) < H)):
            pos[k] = p; k += 1
    return pos, L


print("\n2. Neighbour field at the centre texel (n = z, B along x), %d flakes sampled, cutoff %.0fa" % (NE, RC))
print("   phi   chi_in   b_n/|B| median   90th pct   rms    flatness loss (mean sin^2 tilt)")
for phi in PHIS:
    pos, L = place(phi)
    for chi_in, _ in CHIS:
        M = chi_in                           # units: B_ext = 1, mu0 = 1, so M = chi_in
        bn = np.zeros(NE)
        for i in range(NE):
            d = pos - pos[i]; d -= L * np.round(d / L)
            near = np.where((np.linalg.norm(d, axis=1) < RC) & (np.arange(NF) != i))[0]
            src = (d[near][:, None, :] + rim[None]).reshape(-1, 3)
            q = np.tile(qrim, len(near)) * M
            r = tgt[:, None, :] - src[None]
            r2 = np.maximum(np.sum(r * r, -1), (0.5 * H) ** 2)
            bn[i] = np.mean(np.sum(q[None] * r[..., 2] / (4 * np.pi * r2 ** 1.5), 1))
        a = np.abs(bn)
        print("   %3.0f%%   %4.0f      %8.4f     %8.4f  %8.4f   %.4f" % (100 * phi, chi_in, np.median(a),
              np.percentile(a, 90), np.sqrt(np.mean(bn ** 2)), np.mean(np.sin(np.arctan(bn)) ** 2)))

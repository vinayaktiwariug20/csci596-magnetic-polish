/* Thin-disc overlap test shared by contact.c and dynamics.c (port of overlap_counts() in spike/contact_spike.py).
 * Discs have radius A and thickness H (lengths in units of the flake radius).
 * Overlap: chords on the line where the two planes meet must overlap; near-parallel pairs
 * (|n_i x n_j| <= 0.05) use a slab test. */
#ifndef OVERLAP_H
#define OVERLAP_H
#include <math.h>

#define A 1.0
#define H 0.1

/* d = centre of j minus centre of i (minimum image) */
static int overlap(const double *ni, const double *nj, const double *d) {
    double c[3] = {ni[1] * nj[2] - ni[2] * nj[1], ni[2] * nj[0] - ni[0] * nj[2], ni[0] * nj[1] - ni[1] * nj[0]};
    double cn = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    if (cn > 0.05) {
        double dd[3] = {c[0] / cn, c[1] / cn, c[2] / cn};
        /* point p0 on both planes with dd.p0 = 0: rows ni, nj, dd; right side 0, nj.d, 0 (Cramer's rule) */
        double r = nj[0] * d[0] + nj[1] * d[1] + nj[2] * d[2];
        double det = ni[0] * (nj[1] * dd[2] - nj[2] * dd[1]) - ni[1] * (nj[0] * dd[2] - nj[2] * dd[0]) + ni[2] * (nj[0] * dd[1] - nj[1] * dd[0]);
        double p0[3];
        p0[0] = (-ni[1] * (r * dd[2]) + ni[2] * (r * dd[1])) / det;
        p0[1] = (ni[0] * (r * dd[2]) - ni[2] * (r * dd[0])) / det;
        p0[2] = (-ni[0] * (r * dd[1]) + ni[1] * (r * dd[0])) / det;
        double ti = -(p0[0] * dd[0] + p0[1] * dd[1] + p0[2] * dd[2]);
        double di2 = 0, dj2 = 0, tj = 0;
        for (int k = 0; k < 3; k++) { double q = p0[k] + ti * dd[k]; di2 += q * q; tj += (d[k] - p0[k]) * dd[k]; }
        for (int k = 0; k < 3; k++) { double q = p0[k] + tj * dd[k] - d[k]; dj2 += q * q; }
        if (di2 >= A * A || dj2 >= A * A) return 0;
        return fabs(ti - tj) < sqrt(A * A - di2) + sqrt(A * A - dj2);
    }
    double dz = d[0] * ni[0] + d[1] * ni[1] + d[2] * ni[2], ip2 = 0;
    for (int k = 0; k < 3; k++) { double q = d[k] - dz * ni[k]; ip2 += q * q; }
    return fabs(dz) < H && ip2 < 4 * A * A;
}

/* The slab test is not symmetric in (i, j); always test from the lower index, as the Python spike does. */
static int overlap_ij(int i, int j, const double *ni, const double *nj, const double *d) {
    if (i < j) return overlap(ni, nj, d);
    double m[3] = {-d[0], -d[1], -d[2]};
    return overlap(nj, ni, m);
}

#endif

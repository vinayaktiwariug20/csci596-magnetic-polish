/* Platelet rotation under the switching cat-eye protocol, with or without hard flake-flake contact.
 *
 * N flakes (thin discs, overlap.h) sit at fixed random positions in a periodic cube at volume fraction phi,
 * all in the uniform field of one nail texel.  Each flake rotates as in platelet.c:
 *   dn/dt = -kappa (B.n) (B - (B.n) n),  RK4, dt = 0.05 s, log-normal kappa spread 0.55
 * Protocol: 5 s on B0, 5 s on B90, repeated 6 times (platelet_spike.py).
 *
 * mode ideal: flakes ignore each other (the platelet model).
 * mode hard:  a rotation step is accepted only if it creates no overlap; a blocked flake tries a half and a
 *             quarter step, else stays put.  Positions never change, so blocked flakes cannot slide away:
 *             this is an upper bound on the effect of contact (ideal is the lower bound).
 * Updates are Jacobi-style so the result does not depend on loop order or thread count:
 *   pass 1: candidate for i, checked against every neighbour's OLD orientation;
 *   pass 2: i is reverted if a lower-index neighbour that also moved overlaps i's candidate.
 * Every candidate was checked against all old orientations, so reverting is always safe and the
 * no-overlap invariant holds after every step.
 * Both modes start from the same overlap-free state (random sequential addition, random orientations).
 *
 * Build: gcc -O2 -fopenmp -o dynamics dynamics.c -lm
 * Run:   ./dynamics field_center.txt phi N ideal|hard seed
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "overlap.h"

#ifndef DT
#define DT 0.05
#endif
#define HOLD 5.0
#define CYCLES 6
#define SPREAD 0.55

static unsigned long long rs;
static double urand(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (rs >> 11) * (1.0 / 9007199254740992.0); }
static double grand(void) { double u = urand(), v = urand(); return sqrt(-2 * log(u + 1e-300)) * cos(2 * M_PI * v); }

static void unit(double *v) {
    double s = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (s < 1e-30) s = 1e-30;
    v[0] /= s; v[1] /= s; v[2] /= s;
}

static void rhs(const double *x, const double *B, double kap, double *d) {
    double bx = B[0] * x[0] + B[1] * x[1] + B[2] * x[2];
    for (int k = 0; k < 3; k++) d[k] = -kap * bx * (B[k] - bx * x[k]);
}

/* one RK4 step of length dt from x into y */
static void rk4(const double *x, const double *B, double kap, double dt, double *y) {
    double k1[3], k2[3], k3[3], k4[3], t[3];
    rhs(x, B, kap, k1);
    for (int k = 0; k < 3; k++) t[k] = x[k] + 0.5 * dt * k1[k];
    unit(t); rhs(t, B, kap, k2);
    for (int k = 0; k < 3; k++) t[k] = x[k] + 0.5 * dt * k2[k];
    unit(t); rhs(t, B, kap, k3);
    for (int k = 0; k < 3; k++) t[k] = x[k] + dt * k3[k];
    unit(t); rhs(t, B, kap, k4);
    for (int k = 0; k < 3; k++) y[k] = x[k] + dt / 6.0 * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]);
    unit(y);
}

/* fixed geometry: neighbours within 2A in CSR form, with the separation vector of each pair */
static int N, *nbstart, *nbj;
static double L, *pos, *nbd;
static double *x, *cand, *kap;
static int *moved;

static double min_image(double v) { return v - L * round(v / L); }

static int hits_old(int i, const double *ni) {
    for (int p = nbstart[i]; p < nbstart[i + 1]; p++)
        if (overlap_ij(i, nbj[p], ni, x + 3 * nbj[p], nbd + 3 * p)) return 1;
    return 0;
}

static long long blocked;     /* flake-steps whose full step was rejected */

static void step(const double *B, int hard) {
    long long blk = 0;
    #pragma omp parallel for schedule(static) reduction(+:blk)
    for (int i = 0; i < N; i++) {
        double *c = cand + 3 * i;
        moved[i] = 0;
        for (int tr = 0; tr < (hard ? 3 : 1); tr++) {
            rk4(x + 3 * i, B, kap[i], DT / (1 << tr), c);
            if (!hard || !hits_old(i, c)) { moved[i] = 1; break; }
            if (tr == 0) blk++;
        }
        if (!moved[i]) memcpy(c, x + 3 * i, sizeof(double) * 3);
    }
    if (hard) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; i++) {
            if (!moved[i]) continue;
            for (int p = nbstart[i]; p < nbstart[i + 1]; p++) {
                int j = nbj[p];
                if (j < i && moved[j] && overlap_ij(i, j, cand + 3 * i, cand + 3 * j, nbd + 3 * p)) { moved[i] = -1; break; }
            }
        }
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++)
        if (moved[i] == 1) memcpy(x + 3 * i, cand + 3 * i, sizeof(double) * 3);
    blocked += blk;
}

static void hold(const double *B, int hard) {
    int n = (int)floor(HOLD / DT + 0.5);
    for (int s = 0; s < n; s++) step(B, hard);
}

static double mean_flat(const double *Nn) {
    double s = 0;
    for (int i = 0; i < N; i++) { double c = x[3 * i] * Nn[0] + x[3 * i + 1] * Nn[1] + x[3 * i + 2] * Nn[2]; s += c * c; }
    return s / N;
}

static long long count_overlapping_pairs(void) {
    long long c = 0;
    for (int i = 0; i < N; i++)
        for (int p = nbstart[i]; p < nbstart[i + 1]; p++)
            if (nbj[p] > i && overlap_ij(i, nbj[p], x + 3 * i, x + 3 * nbj[p], nbd + 3 * p)) c++;
    return c;
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: %s field_center.txt phi N ideal|hard seed\n", argv[0]); return 1; }
    double B0[3], B90[3], Nn[3], kappa0;
    FILE *f = fopen(argv[1], "r");
    if (!f || fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", B0, B0 + 1, B0 + 2, B90, B90 + 1, B90 + 2, Nn, Nn + 1, Nn + 2, &kappa0) != 10) {
        fprintf(stderr, "bad field file\n"); return 1;
    }
    fclose(f);
    double phi = atof(argv[2]);
    N = atoi(argv[3]);
    int hard = !strcmp(argv[4], "hard");
    rs = 88172645463325252ULL ^ (unsigned long long)atoll(argv[5]);
    L = cbrt(N * M_PI * A * A * H / phi);

    /* random sequential addition into a cell grid (cell side >= 2A) */
    int nc = (int)floor(L / (2 * A));
    if (nc < 3) { fprintf(stderr, "box too small (L = %g), raise N\n", L); return 1; }
    double cs = L / nc;
    int ncell = nc * nc * nc, *head = malloc(sizeof(int) * ncell), *next = malloc(sizeof(int) * N), *cellof = malloc(sizeof(int) * N);
    for (int c = 0; c < ncell; c++) head[c] = -1;
    pos = malloc(sizeof(double) * 3 * N); x = malloc(sizeof(double) * 3 * N);
    long long tries = 0;
    for (int i = 0; i < N;) {
        if (++tries > 2000LL * N) { fprintf(stderr, "insertion jammed at %d of %d flakes\n", i, N); return 1; }
        double p[3], n[3];
        for (int k = 0; k < 3; k++) { p[k] = urand() * L; n[k] = grand(); }
        unit(n);
        int cx = (int)(p[0] / cs) % nc, cy = (int)(p[1] / cs) % nc, cz = (int)(p[2] / cs) % nc, ok = 1;
        for (int ox = -1; ox <= 1 && ok; ox++) for (int oy = -1; oy <= 1 && ok; oy++) for (int oz = -1; oz <= 1 && ok; oz++) {
            int nb = (((cx + ox + nc) % nc) * nc + (cy + oy + nc) % nc) * nc + (cz + oz + nc) % nc;
            for (int j = head[nb]; j >= 0 && ok; j = next[j]) {
                double d[3], r2 = 0;
                for (int k = 0; k < 3; k++) { d[k] = min_image(pos[3 * j + k] - p[k]); r2 += d[k] * d[k]; }
                if (r2 < 4 * A * A && overlap_ij(i, j, n, x + 3 * j, d)) ok = 0;
            }
        }
        if (!ok) continue;
        memcpy(pos + 3 * i, p, sizeof p); memcpy(x + 3 * i, n, sizeof n);
        cellof[i] = (cx * nc + cy) * nc + cz;
        next[i] = head[cellof[i]]; head[cellof[i]] = i;
        i++;
    }

    /* neighbour lists (fixed, since positions never change) */
    nbstart = malloc(sizeof(int) * (N + 1));
    int cap = 16 * N, np = 0;
    nbj = malloc(sizeof(int) * cap); nbd = malloc(sizeof(double) * 3 * cap);
    for (int i = 0; i < N; i++) {
        nbstart[i] = np;
        int c = cellof[i], cx = c / (nc * nc), cy = (c / nc) % nc, cz = c % nc;
        for (int ox = -1; ox <= 1; ox++) for (int oy = -1; oy <= 1; oy++) for (int oz = -1; oz <= 1; oz++) {
            int nb = (((cx + ox + nc) % nc) * nc + (cy + oy + nc) % nc) * nc + (cz + oz + nc) % nc;
            for (int j = head[nb]; j >= 0; j = next[j]) {
                if (j == i) continue;
                double d[3], r2 = 0;
                for (int k = 0; k < 3; k++) { d[k] = min_image(pos[3 * j + k] - pos[3 * i + k]); r2 += d[k] * d[k]; }
                if (r2 >= 4 * A * A) continue;
                if (np == cap) { cap *= 2; nbj = realloc(nbj, sizeof(int) * cap); nbd = realloc(nbd, sizeof(double) * 3 * cap); }
                nbj[np] = j; memcpy(nbd + 3 * np, d, sizeof d); np++;
            }
        }
    }
    nbstart[N] = np;

    kap = malloc(sizeof(double) * N); cand = malloc(sizeof(double) * 3 * N); moved = malloc(sizeof(int) * N);
    for (int i = 0; i < N; i++) kap[i] = kappa0 * exp(SPREAD * grand());

    double lim[3] = {B0[1] * B90[2] - B0[2] * B90[1], B0[2] * B90[0] - B0[0] * B90[2], B0[0] * B90[1] - B0[1] * B90[0]};
    unit(lim);
    double flim = pow(lim[0] * Nn[0] + lim[1] * Nn[1] + lim[2] * Nn[2], 2);

    double t0 = 0, t1 = 0, f0 = mean_flat(Nn), f1 = 0;
#ifdef _OPENMP
    t0 = omp_get_wtime();
#endif
    for (int cyc = 0; cyc < CYCLES; cyc++) {
        hold(B0, hard);
        if (cyc == 0) f1 = mean_flat(Nn);
        hold(B90, hard);
    }
#ifdef _OPENMP
    t1 = omp_get_wtime();
#endif
    long long steps = (long long)N * CYCLES * 2 * (int)floor(HOLD / DT + 0.5);
    printf("%-5s phi=%.2f N=%d neighbours/flake=%.2f  flat: start %.3f  hold1 %.3f  end %.3f  (limit %.3f)  "
           "blocked %.1f%%  overlaps at end %lld  wall %.1f s\n",
           argv[4], phi, N, (double)np / N, f0, f1, mean_flat(Nn), flim, 100.0 * blocked / steps,
           count_overlapping_pairs(), t1 - t0);
    return 0;
}

/* Platelet rotation AND translation under the switching cat-eye protocol, with hard flake-flake contact
 * resolved by projection (overdamped hard-body dynamics).
 *
 * Contact shape: a "spherodisc", the set of points within H/2 of a zero-thickness disc of radius A
 * (the platelet analogue of a spherocylinder).  Two flakes overlap when the distance D between their
 * zero-thickness discs is below H; the penetration is pen = H - D.
 * Disc-disc distance: for non-crossing discs the closest pair always includes a rim point, so
 *   D = min( min_theta dist(rim_i(theta), disc_j), min_phi dist(rim_j(phi), disc_i) ),
 * each found by sampling the rim and refining with golden-section search.
 *
 * Each time step (lengths in units of A, time in s):
 *   1. free magnetic rotation, RK4 as in platelet.c: dn/dt = -kappa (B.n) (B - (B.n) n);
 *   2. projection sweeps: every overlapping pair is pushed apart along the closest-point direction e.
 *      A push g on flake i at contact point r moves it by the disc mobilities (translation 16 eta A face-on,
 *      (32/3) eta A edge-on, rotation (32/3) eta A^3):  dpos = g_n n + 1.5 g_t,  dtheta = 1.5 r x g.
 *      g is sized so that the pair closes its penetration (both flakes share it by mobility), with
 *      under-relaxation because a flake can have several contacts.  Sweeps are Jacobi (each flake sums its
 *      own pushes; no Newton's third law, no write races) and stop when max pen < TOL * H.
 *   3. tunnelling backstop: a fast flake can pass through a neighbour's edge within one step, and crossed
 *      discs have no meaningful push direction.  Crossings (non-parallel discs whose chords on the line
 *      where their planes meet overlap: the H-independent branch of overlap.h) get no push; after the
 *      sweeps every flake still crossing a neighbour returns to its state at the start of the step.
 *      Reverting only moves towards the previous crossing-free state, so this loop terminates.
 * In the limit dt -> 0 this is infinitely stiff contact with viscous drag, without a stiffness parameter.
 *
 * modes: ideal (no contact), fixed (positions frozen, contacts resolved by rotation only),
 *        free (contacts resolved by rotation and translation).
 *
 * Build: gcc -O2 -fopenmp -o soft soft.c -lm
 * Run:   ./soft field_center.txt phi N ideal|fixed|free seed
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
#define DT 0.01
#endif
#define HOLD 5.0
#define CYCLES 6
#define SPREAD 0.55
#define NSAMP 16
#define NGOLD 12
#define MAXSWEEP 20
#define TOL 0.02
#define RELAX 0.7

static unsigned long long rs;
static double urand(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (rs >> 11) * (1.0 / 9007199254740992.0); }
static double grand(void) { double u = urand(), v = urand(); return sqrt(-2 * log(u + 1e-300)) * cos(2 * M_PI * v); }

static void unit(double *v) {
    double s = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (s < 1e-30) s = 1e-30;
    v[0] /= s; v[1] /= s; v[2] /= s;
}
static double dot(const double *a, const double *b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
static void cross(const double *a, const double *b, double *c) {
    c[0] = a[1] * b[2] - a[2] * b[1]; c[1] = a[2] * b[0] - a[0] * b[2]; c[2] = a[0] * b[1] - a[1] * b[0];
}

static void rhs(const double *x, const double *B, double kap, double *d) {
    double bx = dot(B, x);
    for (int k = 0; k < 3; k++) d[k] = -kap * bx * (B[k] - bx * x[k]);
}

/* one RK4 step of length dt, in place (same arithmetic as dynamics.c) */
static void rk4(double *x, const double *B, double kap, double dt) {
    double k1[3], k2[3], k3[3], k4[3], t[3];
    rhs(x, B, kap, k1);
    for (int k = 0; k < 3; k++) t[k] = x[k] + 0.5 * dt * k1[k];
    unit(t); rhs(t, B, kap, k2);
    for (int k = 0; k < 3; k++) t[k] = x[k] + 0.5 * dt * k2[k];
    unit(t); rhs(t, B, kap, k3);
    for (int k = 0; k < 3; k++) t[k] = x[k] + dt * k3[k];
    unit(t); rhs(t, B, kap, k4);
    for (int k = 0; k < 3; k++) x[k] += dt / 6.0 * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]);
    unit(x);
}

/* orthonormal u, v spanning the plane of normal n */
static void basis(const double *n, double *u, double *v) {
    double t[3] = {1, 0, 0};
    if (fabs(n[0]) > 0.9) { t[0] = 0; t[1] = 1; }
    cross(n, t, u); unit(u); cross(n, u, v);
}

/* closest point q on the disc (centre c, normal n, radius A) to point p */
static void point_disc(const double *p, const double *c, const double *n, double *q) {
    double w[3] = {p[0] - c[0], p[1] - c[1], p[2] - c[2]}, z = dot(w, n);
    for (int k = 0; k < 3; k++) w[k] -= z * n[k];
    double r = sqrt(dot(w, w)), s = r > A ? A / r : 1.0;
    for (int k = 0; k < 3; k++) q[k] = c[k] + s * w[k];
}

static double cs_tab[NSAMP], sn_tab[NSAMP];

/* squared distance from the rim point (cos, sin) of disc (c1, u1, v1) to disc (c2, n2); p on rim, q on disc 2 */
static double rim_d2(double cs, double sn, const double *c1, const double *u1, const double *v1, const double *c2, const double *n2, double *p, double *q) {
    for (int k = 0; k < 3; k++) p[k] = c1[k] + A * (cs * u1[k] + sn * v1[k]);
    point_disc(p, c2, n2, q);
    double d[3] = {q[0] - p[0], q[1] - p[1], q[2] - p[2]};
    return dot(d, d);
}

/* min over the rim of disc 1 of the distance to disc 2 (sample, then golden-section refine) */
static double rim_min(const double *c1, const double *n1, const double *c2, const double *n2, double *pb, double *qb) {
    double u[3], v[3], p[3], q[3], best = 1e300;
    int sb = 0;
    basis(n1, u, v);
    for (int s = 0; s < NSAMP; s++) {
        double d2 = rim_d2(cs_tab[s], sn_tab[s], c1, u, v, c2, n2, p, q);
        if (d2 < best) { best = d2; sb = s; }
    }
    const double g = 0.6180339887498949, w = 2 * M_PI / NSAMP, tb = sb * w;
    double lo = tb - w, hi = tb + w, x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
    double f1 = rim_d2(cos(x1), sin(x1), c1, u, v, c2, n2, p, q), f2 = rim_d2(cos(x2), sin(x2), c1, u, v, c2, n2, p, q);
    for (int it = 0; it < NGOLD; it++) {
        if (f1 < f2) { hi = x2; x2 = x1; f2 = f1; x1 = hi - g * (hi - lo); f1 = rim_d2(cos(x1), sin(x1), c1, u, v, c2, n2, p, q); }
        else { lo = x1; x1 = x2; f1 = f2; x2 = lo + g * (hi - lo); f2 = rim_d2(cos(x2), sin(x2), c1, u, v, c2, n2, p, q); }
    }
    double th = 0.5 * (lo + hi), d2 = rim_d2(cos(th), sin(th), c1, u, v, c2, n2, pb, qb);
    if (best < d2) d2 = rim_d2(cs_tab[sb], sn_tab[sb], c1, u, v, c2, n2, pb, qb);
    return sqrt(d2);
}

/* distance between disc i (centre at origin) and disc j (centre d); p on disc i, q on disc j */
static double disc_dist(const double *ni, const double *d, const double *nj, double *p, double *q) {
    static const double o[3] = {0, 0, 0};
    double p1[3], q1[3], p2[3], q2[3];
    double d1 = rim_min(o, ni, d, nj, p1, q1);
    double d2 = rim_min(d, nj, o, ni, p2, q2);
    if (d1 <= d2) { memcpy(p, p1, sizeof p1); memcpy(q, q1, sizeof q1); return d1; }
    memcpy(p, q2, sizeof q2); memcpy(q, p2, sizeof p2); return d2;
}

/* cheap necessary condition for D < H: disc j must reach within H of the plane of i, and vice versa */
static int may_touch(const double *ni, const double *d, const double *nj) {
    double c[3];
    cross(ni, nj, c);
    double s = sqrt(dot(c, c)) * A + H;
    return fabs(dot(d, ni)) < s && fabs(dot(d, nj)) < s;
}

static int N, mode, nc, ncell, *head, *next, *flag;
static double L, cs, *pos, *x, *dpos, *drot, *kap, *xold, *posold;

/* zero-thickness discs cut through each other (near-parallel pairs never count: their distance is reliable) */
static int crossing(const double *ni, const double *nj, const double *d) {
    double c[3];
    cross(ni, nj, c);
    return dot(c, c) > 0.05 * 0.05 && overlap(ni, nj, d);
}

static double min_image(double v) { return v - L * round(v / L); }

static void build_cells(void) {
    for (int c = 0; c < ncell; c++) head[c] = -1;
    for (int i = 0; i < N; i++) {
        int c = (((int)(pos[3 * i] / cs) % nc) * nc + (int)(pos[3 * i + 1] / cs) % nc) * nc + (int)(pos[3 * i + 2] / cs) % nc;
        next[i] = head[c]; head[c] = i;
    }
}

/* how far a unit push along e at point r on a flake with normal n closes the gap (mobility along e) */
static double mobility(const double *n, const double *r, const double *e, int translate) {
    double m = 0, t[3];
    if (translate) { double en = dot(e, n); m += en * en + 1.5 * (1 - en * en); }
    cross(r, e, t);
    double tn = dot(t, n);
    m += 1.5 * (dot(t, t) - tn * tn);     /* spin about n does not move a disc */
    return m;
}

/* one Jacobi projection sweep; returns the largest penetration found (before the corrections) */
static double sweep(long long *ncont, long long *ncross) {
    double mp = 0;
    long long nct = 0, ncr = 0;
    #pragma omp parallel for schedule(dynamic, 64) reduction(max:mp) reduction(+:nct, ncr)
    for (int i = 0; i < N; i++) {
        double DP[3] = {0, 0, 0}, DR[3] = {0, 0, 0}, *ni = x + 3 * i;
        int cx = (int)(pos[3 * i] / cs) % nc, cy = (int)(pos[3 * i + 1] / cs) % nc, cz = (int)(pos[3 * i + 2] / cs) % nc;
        for (int ox = -1; ox <= 1; ox++) for (int oy = -1; oy <= 1; oy++) for (int oz = -1; oz <= 1; oz++) {
            int nb = (((cx + ox + nc) % nc) * nc + (cy + oy + nc) % nc) * nc + (cz + oz + nc) % nc;
            for (int j = head[nb]; j >= 0; j = next[j]) {
                if (j == i) continue;
                double d[3], r2 = 0, *nj = x + 3 * j;
                for (int k = 0; k < 3; k++) { d[k] = min_image(pos[3 * j + k] - pos[3 * i + k]); r2 += d[k] * d[k]; }
                if (r2 >= (2 * A + H) * (2 * A + H) || !may_touch(ni, d, nj)) continue;
                if (crossing(ni, nj, d)) { ncr++; continue; }        /* left to the backstop */
                double p[3], q[3], D = disc_dist(ni, d, nj, p, q);
                if (D >= H) continue;
                double e[3], qj[3];
                for (int k = 0; k < 3; k++) e[k] = (q[k] - p[k]) / D;
                for (int k = 0; k < 3; k++) qj[k] = q[k] - d[k];      /* contact point relative to centre j */
                double pen = H - D;
                double w = mobility(ni, p, e, mode == 2) + mobility(nj, qj, e, mode == 2);
                double gm = RELAX * pen / (w > 1e-12 ? w : 1e-12), g[3], t[3];
                for (int k = 0; k < 3; k++) g[k] = -gm * e[k];         /* push on i, away from j */
                if (mode == 2) { double gn = dot(g, ni); for (int k = 0; k < 3; k++) DP[k] += gn * ni[k] + 1.5 * (g[k] - gn * ni[k]); }
                cross(p, g, t);
                for (int k = 0; k < 3; k++) DR[k] += 1.5 * t[k];
                nct++; if (pen > mp) mp = pen;
            }
        }
        memcpy(dpos + 3 * i, DP, sizeof DP); memcpy(drot + 3 * i, DR, sizeof DR);
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        double *n = x + 3 * i, wn[3];
        cross(drot + 3 * i, n, wn);
        for (int k = 0; k < 3; k++) n[k] += wn[k];
        unit(n);
        if (mode == 2) for (int k = 0; k < 3; k++) { double y = pos[3 * i + k] + dpos[3 * i + k]; pos[3 * i + k] = y - L * floor(y / L); }
    }
    *ncont += nct; *ncross += ncr;
    return mp;
}

/* send every flake that crosses a neighbour back to its start-of-step state; returns the number reverted */
static long long revert_crossings(void) {
    long long total = 0;
    for (int pass = 0; pass < 50; pass++) {
        build_cells();
        long long nf = 0;
        #pragma omp parallel for schedule(dynamic, 64) reduction(+:nf)
        for (int i = 0; i < N; i++) {
            flag[i] = 0;
            int cx = (int)(pos[3 * i] / cs) % nc, cy = (int)(pos[3 * i + 1] / cs) % nc, cz = (int)(pos[3 * i + 2] / cs) % nc;
            for (int ox = -1; ox <= 1 && !flag[i]; ox++) for (int oy = -1; oy <= 1 && !flag[i]; oy++) for (int oz = -1; oz <= 1 && !flag[i]; oz++) {
                int nb = (((cx + ox + nc) % nc) * nc + (cy + oy + nc) % nc) * nc + (cz + oz + nc) % nc;
                for (int j = head[nb]; j >= 0; j = next[j]) {
                    if (j == i) continue;
                    double d[3], r2 = 0;
                    for (int k = 0; k < 3; k++) { d[k] = min_image(pos[3 * j + k] - pos[3 * i + k]); r2 += d[k] * d[k]; }
                    if (r2 < 4 * A * A && crossing(x + 3 * i, x + 3 * j, d)) { flag[i] = 1; break; }
                }
            }
            nf += flag[i];
        }
        if (!nf) return total;
        for (int i = 0; i < N; i++)
            if (flag[i]) { memcpy(x + 3 * i, xold + 3 * i, sizeof(double) * 3); memcpy(pos + 3 * i, posold + 3 * i, sizeof(double) * 3); }
        total += nf;
    }
    fprintf(stderr, "crossings not cleared after 50 revert passes\n");
    return total;
}

static double mean_flat(const double *Nn) {
    double s = 0;
    for (int i = 0; i < N; i++) { double c = dot(x + 3 * i, Nn); s += c * c; }
    return s / N;
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: %s field_center.txt phi N ideal|fixed|free seed\n", argv[0]); return 1; }
    double B0[3], B90[3], Nn[3], kappa0;
    FILE *f = fopen(argv[1], "r");
    if (!f || fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", B0, B0 + 1, B0 + 2, B90, B90 + 1, B90 + 2, Nn, Nn + 1, Nn + 2, &kappa0) != 10) {
        fprintf(stderr, "bad field file\n"); return 1;
    }
    fclose(f);
    double phi = atof(argv[2]);
    N = atoi(argv[3]);
    mode = !strcmp(argv[4], "free") ? 2 : !strcmp(argv[4], "fixed") ? 1 : 0;
    rs = 88172645463325252ULL ^ (unsigned long long)atoll(argv[5]);
    L = cbrt(N * M_PI * A * A * H / phi);
    for (int s = 0; s < NSAMP; s++) { cs_tab[s] = cos(2 * M_PI * s / NSAMP); sn_tab[s] = sin(2 * M_PI * s / NSAMP); }

    /* random sequential addition: same overlap test, cell size and RNG order as dynamics.c */
    nc = (int)floor(L / (2 * A));
    if (nc < 3) { fprintf(stderr, "box too small (L = %g), raise N\n", L); return 1; }
    cs = L / nc; ncell = nc * nc * nc;
    head = malloc(sizeof(int) * ncell); next = malloc(sizeof(int) * N);
    for (int c = 0; c < ncell; c++) head[c] = -1;
    pos = malloc(sizeof(double) * 3 * N); x = malloc(sizeof(double) * 3 * N);
    dpos = malloc(sizeof(double) * 3 * N); drot = malloc(sizeof(double) * 3 * N); kap = malloc(sizeof(double) * N);
    xold = malloc(sizeof(double) * 3 * N); posold = malloc(sizeof(double) * 3 * N); flag = malloc(sizeof(int) * N);
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
        int c = (cx * nc + cy) * nc + cz;
        next[i] = head[c]; head[c] = i;
        i++;
    }
    for (int i = 0; i < N; i++) kap[i] = kappa0 * exp(SPREAD * grand());

    /* from here on cells must hold spherodisc contacts: side >= 2A + H */
    nc = (int)floor(L / (2 * A + H));
    if (nc < 3) { fprintf(stderr, "box too small (L = %g), raise N\n", L); return 1; }
    cs = L / nc; ncell = nc * nc * nc;
    free(head); head = malloc(sizeof(int) * ncell);
    double *pos0 = malloc(sizeof(double) * 3 * N);
    memcpy(pos0, pos, sizeof(double) * 3 * N);

    double t0 = 0, t1 = 0, fs = mean_flat(Nn), f1 = 0, worst = 0;
    long long ncont = 0, ncross = 0, nsweeps = 0, unconverged = 0, nsteps = 0, nrevert = 0;
#ifdef _OPENMP
    t0 = omp_get_wtime();
#endif
    int nhold = (int)floor(HOLD / DT + 0.5);
    for (int cyc = 0; cyc < CYCLES; cyc++) {
        for (int ph = 0; ph < 2; ph++) {
            const double *B = ph ? B90 : B0;
            for (int s = 0; s < nhold; s++, nsteps++) {
                if (mode) { memcpy(xold, x, sizeof(double) * 3 * N); memcpy(posold, pos, sizeof(double) * 3 * N); }
                #pragma omp parallel for schedule(static)
                for (int i = 0; i < N; i++) rk4(x + 3 * i, B, kap[i], DT);
                if (!mode) continue;
                int it;
                double mp = 0;
                for (it = 0; it < MAXSWEEP; it++) {
                    if (mode == 2 || it == 0) build_cells();
                    mp = sweep(&ncont, &ncross);
                    nsweeps++;
                    if (mp < TOL * H) break;
                }
                if (it == MAXSWEEP) { unconverged++; if (mp > worst) worst = mp; }
                nrevert += revert_crossings();
            }
            if (cyc == 0 && ph == 0) f1 = mean_flat(Nn);
        }
    }
#ifdef _OPENMP
    t1 = omp_get_wtime();
#endif
    double disp = 0;
    for (int i = 0; i < N; i++) {
        double r2 = 0;
        for (int k = 0; k < 3; k++) { double d = min_image(pos[3 * i + k] - pos0[3 * i + k]); r2 += d * d; }
        disp += sqrt(r2);
    }
    printf("%-5s phi=%.2f N=%d dt=%g  flat: start %.3f  hold1 %.3f  end %.3f  sweeps/step %.1f  unconverged steps %lld (worst pen %.2fH)  "
           "crossings seen %lld  reverted flake-steps %.3f%%  mean shift %.2fa  wall %.1f s\n",
           argv[4], phi, N, DT, fs, f1, mean_flat(Nn), (double)nsweeps / nsteps, unconverged, worst / H, ncross, 100.0 * nrevert / ((double)N * nsteps), disp / N, t1 - t0);
    return 0;
}

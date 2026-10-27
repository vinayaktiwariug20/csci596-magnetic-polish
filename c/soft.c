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
 * Box: a periodic cube by default.  With a film thickness Lz (in units of A) the box is periodic in x and y
 * and bounded by two flat walls, the nail at z = 0 and the film top at z = Lz; the fields are first rotated
 * so that the local nail normal is z.
 * A flake keeps its lowest (highest) rim point H/2 from a wall; walls take part in the projection sweeps
 * but do not move, so the flake takes the whole correction.  Flatness is then also reported by height.
 *
 * Build: gcc -O2 -fopenmp -o soft soft.c -lm
 * Run:   ./soft field_center.txt phi N ideal|fixed|free seed [film_thickness_in_a]
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

static int N, mode, wall, ncv[3], ncell, *head, *next, *flag;
static double Lb[3], csv[3], *pos, *x, *dpos, *drot, *kap, *xold, *posold;

/* zero-thickness discs cut through each other (near-parallel pairs never count: their distance is reliable) */
static int crossing(const double *ni, const double *nj, const double *d) {
    double c[3];
    cross(ni, nj, c);
    return dot(c, c) > 0.05 * 0.05 && overlap(ni, nj, d);
}

/* minimum image along axis k; z is not periodic when there are walls */
static double mi(double v, int k) { return (wall && k == 2) ? v : v - Lb[k] * round(v / Lb[k]); }

static void cell_of(const double *p, int *c) {
    for (int k = 0; k < 3; k++) {
        c[k] = (int)(p[k] / csv[k]);
        if (c[k] >= ncv[k]) c[k] = ncv[k] - 1;
        if (c[k] < 0) c[k] = 0;
    }
}

/* neighbouring cell (ox, oy, oz) of cell c, or -1 beyond a wall */
static int nb_cell(const int *c, int ox, int oy, int oz) {
    int z = c[2] + oz;
    if (wall) { if (z < 0 || z >= ncv[2]) return -1; }
    else z = (z + ncv[2]) % ncv[2];
    return (((c[0] + ox + ncv[0]) % ncv[0]) * ncv[1] + (c[1] + oy + ncv[1]) % ncv[1]) * ncv[2] + z;
}

static void set_grid(double side) {
    ncell = 1;
    for (int k = 0; k < 3; k++) {
        ncv[k] = (int)floor(Lb[k] / side);
        if (ncv[k] < 1) ncv[k] = 1;
        if (ncv[k] < 3 && !(wall && k == 2)) { fprintf(stderr, "box too small (L = %g), raise N\n", Lb[k]); exit(1); }
        csv[k] = Lb[k] / ncv[k]; ncell *= ncv[k];
    }
}

static void build_cells(void) {
    for (int c = 0; c < ncell; c++) head[c] = -1;
    for (int i = 0; i < N; i++) {
        int c[3];
        cell_of(pos + 3 * i, c);
        int id = (c[0] * ncv[1] + c[1]) * ncv[2] + c[2];
        next[i] = head[id]; head[id] = i;
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

/* push g on a flake (normal n) at contact point r: accumulate displacement DP and rotation DR */
static void apply_push(const double *n, const double *r, const double *g, double *DP, double *DR) {
    double t[3];
    if (mode == 2) { double gn = dot(g, n); for (int k = 0; k < 3; k++) DP[k] += gn * n[k] + 1.5 * (g[k] - gn * n[k]); }
    cross(r, g, t);
    for (int k = 0; k < 3; k++) DR[k] += 1.5 * t[k];
}

/* penetration of a flake (centre height z, normal n) into the wall on side sd (-1: nail at z = 0,
 * +1: film top at z = Lz); r = lowest (highest) rim point relative to the centre */
static double wall_pen(double z, const double *n, int sd, double *r) {
    double s = sqrt(fmax(0.0, 1 - n[2] * n[2]));
    if (s > 1e-9) { for (int k = 0; k < 3; k++) r[k] = A * (sd * (k == 2) - sd * n[2] * n[k]) / s; }
    else r[0] = r[1] = r[2] = 0;
    return sd < 0 ? H / 2 - (z - A * s) : H / 2 - (Lb[2] - z - A * s);
}

/* one Jacobi projection sweep; returns the largest penetration found (before the corrections) */
static double sweep(long long *ncont, long long *ncross) {
    double mp = 0;
    long long nct = 0, ncr = 0;
    #pragma omp parallel for schedule(dynamic, 64) reduction(max:mp) reduction(+:nct, ncr)
    for (int i = 0; i < N; i++) {
        double DP[3] = {0, 0, 0}, DR[3] = {0, 0, 0}, *ni = x + 3 * i;
        int c[3];
        cell_of(pos + 3 * i, c);
        for (int ox = -1; ox <= 1; ox++) for (int oy = -1; oy <= 1; oy++) for (int oz = -1; oz <= 1; oz++) {
            int nb = nb_cell(c, ox, oy, oz);
            if (nb < 0) continue;
            for (int j = head[nb]; j >= 0; j = next[j]) {
                if (j == i) continue;
                double d[3], r2 = 0, *nj = x + 3 * j;
                for (int k = 0; k < 3; k++) { d[k] = mi(pos[3 * j + k] - pos[3 * i + k], k); r2 += d[k] * d[k]; }
                if (r2 >= (2 * A + H) * (2 * A + H) || !may_touch(ni, d, nj)) continue;
                if (crossing(ni, nj, d)) { ncr++; continue; }        /* left to the backstop */
                double p[3], q[3], D = disc_dist(ni, d, nj, p, q);
                if (D >= H) continue;
                double e[3], qj[3];
                for (int k = 0; k < 3; k++) e[k] = (q[k] - p[k]) / D;
                for (int k = 0; k < 3; k++) qj[k] = q[k] - d[k];      /* contact point relative to centre j */
                double pen = H - D;
                double w = mobility(ni, p, e, mode == 2) + mobility(nj, qj, e, mode == 2);
                double gm = RELAX * pen / (w > 1e-12 ? w : 1e-12), g[3];
                for (int k = 0; k < 3; k++) g[k] = -gm * e[k];         /* push on i, away from j */
                apply_push(ni, p, g, DP, DR);
                nct++; if (pen > mp) mp = pen;
            }
        }
        /* walls do not move, so the flake takes the whole correction */
        for (int sd = -1; wall && sd <= 1; sd += 2) {
            double r[3], pen = wall_pen(pos[3 * i + 2], ni, sd, r);
            if (pen <= 0) continue;
            double e[3] = {0, 0, sd}, w = mobility(ni, r, e, mode == 2);
            if (w < 1e-12) continue;      /* frozen flat flake on the wall: nothing can move */
            double g[3] = {0, 0, -sd * RELAX * pen / w};
            apply_push(ni, r, g, DP, DR);
            nct++; if (pen > mp) mp = pen;
        }
        memcpy(dpos + 3 * i, DP, sizeof DP); memcpy(drot + 3 * i, DR, sizeof DR);
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        double *n = x + 3 * i, wn[3];
        cross(drot + 3 * i, n, wn);
        for (int k = 0; k < 3; k++) n[k] += wn[k];
        unit(n);
        if (mode == 2) for (int k = 0; k < 3; k++) {
            double y = pos[3 * i + k] + dpos[3 * i + k];
            pos[3 * i + k] = (wall && k == 2) ? y : y - Lb[k] * floor(y / Lb[k]);
        }
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
            int c[3];
            cell_of(pos + 3 * i, c);
            for (int ox = -1; ox <= 1 && !flag[i]; ox++) for (int oy = -1; oy <= 1 && !flag[i]; oy++) for (int oz = -1; oz <= 1 && !flag[i]; oz++) {
                int nb = nb_cell(c, ox, oy, oz);
                if (nb < 0) continue;
                for (int j = head[nb]; j >= 0; j = next[j]) {
                    if (j == i) continue;
                    double d[3], r2 = 0;
                    for (int k = 0; k < 3; k++) { d[k] = mi(pos[3 * j + k] - pos[3 * i + k], k); r2 += d[k] * d[k]; }
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

/* flatness in height bins of width 2A (one flake diameter) above the nail; wall runs only */
static void flat_profile(const double *Nn, const char *label) {
    int nb = (int)ceil(Lb[2] / (2 * A));
    printf("  %s flatness by height above the nail (count):", label);
    for (int b = 0; b < nb; b++) {
        double s = 0;
        int n = 0;
        for (int i = 0; i < N; i++) {
            if ((int)(pos[3 * i + 2] / (2 * A)) != b) continue;
            double c = dot(x + 3 * i, Nn); s += c * c; n++;
        }
        printf("  %g-%ga %.3f (%d)", 2.0 * b, fmin(2.0 * b + 2, Lb[2]), n ? s / n : 0.0, n);
    }
    printf("\n");
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: %s field_center.txt phi N ideal|fixed|free seed [film_thickness_in_a]\n", argv[0]); return 1; }
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
    double Lz = argc > 6 ? atof(argv[6]) : 0;
    wall = Lz > 0;
    if (wall) {
        /* rotate the fields so that the local nail normal is z (Rodrigues: R u = u + v x u + v x (v x u) / (1 + c)) */
        unit(Nn);
        double z[3] = {0, 0, 1}, v[3], c = dot(Nn, z), *vec[3] = {B0, B90, Nn};
        if (c < -0.999) { fprintf(stderr, "nail normal points along -z\n"); return 1; }
        cross(Nn, z, v);
        for (int m = 0; m < 3; m++) {
            double *u = vec[m], vu[3], vvu[3];
            cross(v, u, vu); cross(v, vu, vvu);
            for (int k = 0; k < 3; k++) u[k] += vu[k] + vvu[k] / (1 + c);
        }
        Lb[0] = Lb[1] = sqrt(N * M_PI * A * A * H / (phi * Lz)); Lb[2] = Lz;
    } else Lb[0] = Lb[1] = Lb[2] = cbrt(N * M_PI * A * A * H / phi);
    for (int s = 0; s < NSAMP; s++) { cs_tab[s] = cos(2 * M_PI * s / NSAMP); sn_tab[s] = sin(2 * M_PI * s / NSAMP); }

    /* random sequential addition: same overlap test, cell size and RNG order as dynamics.c */
    set_grid(2 * A);
    head = malloc(sizeof(int) * ncell); next = malloc(sizeof(int) * N);
    for (int c = 0; c < ncell; c++) head[c] = -1;
    pos = malloc(sizeof(double) * 3 * N); x = malloc(sizeof(double) * 3 * N);
    dpos = malloc(sizeof(double) * 3 * N); drot = malloc(sizeof(double) * 3 * N); kap = malloc(sizeof(double) * N);
    xold = malloc(sizeof(double) * 3 * N); posold = malloc(sizeof(double) * 3 * N); flag = malloc(sizeof(int) * N);
    long long tries = 0;
    for (int i = 0; i < N;) {
        if (++tries > 2000LL * N) { fprintf(stderr, "insertion jammed at %d of %d flakes\n", i, N); return 1; }
        double p[3], n[3], r[3];
        for (int k = 0; k < 3; k++) { p[k] = urand() * Lb[k]; n[k] = grand(); }
        unit(n);
        if (wall && (wall_pen(p[2], n, -1, r) > 0 || wall_pen(p[2], n, 1, r) > 0)) continue;
        int c[3], ok = 1;
        cell_of(p, c);
        for (int ox = -1; ox <= 1 && ok; ox++) for (int oy = -1; oy <= 1 && ok; oy++) for (int oz = -1; oz <= 1 && ok; oz++) {
            int nb = nb_cell(c, ox, oy, oz);
            if (nb < 0) continue;
            for (int j = head[nb]; j >= 0 && ok; j = next[j]) {
                double d[3], r2 = 0;
                for (int k = 0; k < 3; k++) { d[k] = mi(pos[3 * j + k] - p[k], k); r2 += d[k] * d[k]; }
                if (r2 < 4 * A * A && overlap_ij(i, j, n, x + 3 * j, d)) ok = 0;
            }
        }
        if (!ok) continue;
        memcpy(pos + 3 * i, p, sizeof p); memcpy(x + 3 * i, n, sizeof n);
        int id = (c[0] * ncv[1] + c[1]) * ncv[2] + c[2];
        next[i] = head[id]; head[id] = i;
        i++;
    }
    for (int i = 0; i < N; i++) kap[i] = kappa0 * exp(SPREAD * grand());

    /* from here on cells must hold spherodisc contacts: side >= 2A + H */
    set_grid(2 * A + H);
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
            if (cyc == 0 && ph == 0) { f1 = mean_flat(Nn); if (wall) flat_profile(Nn, "after hold 1"); }
        }
    }
#ifdef _OPENMP
    t1 = omp_get_wtime();
#endif
    double disp = 0;
    for (int i = 0; i < N; i++) {
        double r2 = 0;
        for (int k = 0; k < 3; k++) { double d = mi(pos[3 * i + k] - pos0[3 * i + k], k); r2 += d * d; }
        disp += sqrt(r2);
    }
    if (wall) flat_profile(Nn, "at end      ");
    printf("%-5s phi=%.2f N=%d dt=%g film=%ga  flat: start %.3f  hold1 %.3f  end %.3f  sweeps/step %.1f  unconverged steps %lld (worst pen %.2fH)  "
           "crossings seen %lld  reverted flake-steps %.3f%%  mean shift %.2fa  wall %.1f s\n",
           argv[4], phi, N, DT, Lz, fs, f1, mean_flat(Nn), (double)nsweeps / nsteps, unconverged, worst / H, ncross, 100.0 * nrevert / ((double)N * nsteps), disp / N, t1 - t0);
    return 0;
}

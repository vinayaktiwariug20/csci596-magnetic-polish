/* Platelet orientation model, C port of spike/platelet_spike.py (platelet branch only).
 *   dn/dt = -kappa (B.n) (B - (B.n) n),  RK4 with renormalisation, dt = 0.05 s
 * Protocol: 5 s on B0, 5 s on B90, repeated 6 times.  Texels are independent, so the outer loop is the
 * parallel loop (OpenMP if compiled with -fopenmp).
 * Build: gcc -O2 -fopenmp -o platelet platelet.c -lm      Run: ./platelet inputs.bin out.bin
 * out.bin = cnt doubles (mean (n.N)^2 after hold 1), cnt doubles (same, after all cycles). */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define DT 0.05
#define HOLD 5.0
#define CYCLES 6

static void unit(double *v) {
    double s = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (s < 1e-30) s = 1e-30;
    v[0] /= s; v[1] /= s; v[2] /= s;
}

static void rhs(const double *x, const double *B, double kap, double *d) {
    double bx = B[0] * x[0] + B[1] * x[1] + B[2] * x[2];
    for (int k = 0; k < 3; k++) d[k] = -kap * bx * (B[k] - bx * x[k]);
}

static void step(double *x, const double *B, double kap) {
    double k1[3], k2[3], k3[3], k4[3], t[3];
    rhs(x, B, kap, k1);
    for (int k = 0; k < 3; k++) t[k] = x[k] + 0.5 * DT * k1[k];
    unit(t); rhs(t, B, kap, k2);
    for (int k = 0; k < 3; k++) t[k] = x[k] + 0.5 * DT * k2[k];
    unit(t); rhs(t, B, kap, k3);
    for (int k = 0; k < 3; k++) t[k] = x[k] + DT * k3[k];
    unit(t); rhs(t, B, kap, k4);
    for (int k = 0; k < 3; k++) x[k] += DT / 6.0 * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]);
    unit(x);
}

static void hold(double *x, const double *B, double kap) {
    int n = (int)floor(HOLD / DT + 0.5);
    for (int s = 0; s < n; s++) step(x, B, kap);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s inputs.bin out.bin\n", argv[0]); return 1; }
    FILE *f = fopen(argv[1], "rb");
    if (!f) { perror(argv[1]); return 1; }
    int hdr[2];
    if (fread(hdr, sizeof(int), 2, f) != 2) return 1;
    int cnt = hdr[0], NF = hdr[1];
    double *N = malloc(sizeof(double) * 3 * cnt), *B0 = malloc(sizeof(double) * 3 * cnt), *B90 = malloc(sizeof(double) * 3 * cnt);
    double *kap = malloc(sizeof(double) * cnt * NF), *x = malloc(sizeof(double) * 3 * cnt * NF);
    if (fread(N, 8, 3 * cnt, f) != 3 * cnt || fread(B0, 8, 3 * cnt, f) != 3 * cnt || fread(B90, 8, 3 * cnt, f) != 3 * cnt ||
        fread(kap, 8, cnt * NF, f) != (size_t)cnt * NF || fread(x, 8, 3 * (size_t)cnt * NF, f) != 3 * (size_t)cnt * NF) {
        fprintf(stderr, "short read\n"); return 1;
    }
    fclose(f);
    double *first = malloc(sizeof(double) * cnt), *end = malloc(sizeof(double) * cnt);

    double t0 = 0;
#ifdef _OPENMP
    t0 = omp_get_wtime();
#endif
    #pragma omp parallel for schedule(static)
    for (int t = 0; t < cnt; t++) {
        double sf = 0, se = 0;
        for (int i = 0; i < NF; i++) {
            double *xi = x + 3 * ((size_t)t * NF + i), k = kap[(size_t)t * NF + i];
            hold(xi, B0 + 3 * t, k);
            double c = xi[0] * N[3 * t] + xi[1] * N[3 * t + 1] + xi[2] * N[3 * t + 2];
            sf += c * c;
            hold(xi, B90 + 3 * t, k);
            for (int cyc = 1; cyc < CYCLES; cyc++) { hold(xi, B0 + 3 * t, k); hold(xi, B90 + 3 * t, k); }
            c = xi[0] * N[3 * t] + xi[1] * N[3 * t + 1] + xi[2] * N[3 * t + 2];
            se += c * c;
        }
        first[t] = sf / NF; end[t] = se / NF;
    }
#ifdef _OPENMP
    printf("threads = %d, wall = %.3f s\n", omp_get_max_threads(), omp_get_wtime() - t0);
#endif
    FILE *o = fopen(argv[2], "wb");
    fwrite(first, 8, cnt, o); fwrite(end, 8, cnt, o); fclose(o);
    return 0;
}

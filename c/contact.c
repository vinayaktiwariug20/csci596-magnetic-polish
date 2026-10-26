/* Flake-flake contact count, C port of overlap_counts() in spike/contact_spike.py.
 * Flakes are thin discs (radius A = 1, thickness H = 0.1) in a periodic cube of side L.
 * Overlap test in overlap.h.  Neighbours come from a cell list with cell side >= 2A.
 * Each flake counts all its partners itself (no Newton's third law), so the OpenMP loop has no races.
 *
 * Build: gcc -O2 -fopenmp -o contact contact.c -lm
 * Run:   ./contact file in.bin out.bin       in.bin = int32 N, double L, N*3 pos, N*3 normals; out = N int32 counts
 *        ./contact gen N phi seed            random positions and random normals, prints the mean count
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "overlap.h"

/* counts[i] = number of flakes overlapping flake i */
static void count_overlaps(int N, double L, const double *pos, const double *nrm, int *counts) {
    int nc = (int)floor(L / (2 * A));
    if (nc < 3) { fprintf(stderr, "box too small for cell list (L = %g)\n", L); exit(1); }
    double cs = L / nc;
    int ncell = nc * nc * nc;
    int *head = malloc(sizeof(int) * ncell), *next = malloc(sizeof(int) * N), *cell = malloc(sizeof(int) * N);
    for (int c = 0; c < ncell; c++) head[c] = -1;
    for (int i = 0; i < N; i++) {
        int cx = (int)(pos[3 * i] / cs) % nc, cy = (int)(pos[3 * i + 1] / cs) % nc, cz = (int)(pos[3 * i + 2] / cs) % nc;
        cell[i] = (cx * nc + cy) * nc + cz;
        next[i] = head[cell[i]]; head[cell[i]] = i;
    }
    #pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < N; i++) {
        int c = cell[i], cx = c / (nc * nc), cy = (c / nc) % nc, cz = c % nc, cnt = 0;
        for (int ox = -1; ox <= 1; ox++) for (int oy = -1; oy <= 1; oy++) for (int oz = -1; oz <= 1; oz++) {
            int nb = (((cx + ox + nc) % nc) * nc + (cy + oy + nc) % nc) * nc + (cz + oz + nc) % nc;
            for (int j = head[nb]; j >= 0; j = next[j]) {
                if (j == i) continue;
                double d[3], r2 = 0;
                for (int k = 0; k < 3; k++) { d[k] = pos[3 * j + k] - pos[3 * i + k]; d[k] -= L * round(d[k] / L); r2 += d[k] * d[k]; }
                if (r2 < 4 * A * A) cnt += overlap_ij(i, j, nrm + 3 * i, nrm + 3 * j, d);
            }
        }
        counts[i] = cnt;
    }
    free(head); free(next); free(cell);
}

static unsigned long long rs;
static double urand(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (rs >> 11) * (1.0 / 9007199254740992.0); }

int main(int argc, char **argv) {
    int N; double L, *pos, *nrm;
    if (argc >= 4 && !strcmp(argv[1], "file")) {
        FILE *f = fopen(argv[2], "rb");
        if (!f || fread(&N, 4, 1, f) != 1 || fread(&L, 8, 1, f) != 1) { fprintf(stderr, "bad input\n"); return 1; }
        pos = malloc(sizeof(double) * 3 * N); nrm = malloc(sizeof(double) * 3 * N);
        if (fread(pos, 8, 3 * N, f) != 3 * (size_t)N || fread(nrm, 8, 3 * N, f) != 3 * (size_t)N) { fprintf(stderr, "short read\n"); return 1; }
        fclose(f);
    } else if (argc >= 5 && !strcmp(argv[1], "gen")) {
        N = atoi(argv[2]); double phi = atof(argv[3]); rs = 88172645463325252ULL ^ (unsigned long long)atoll(argv[4]);
        L = cbrt(N * M_PI * A * A * H / phi);
        pos = malloc(sizeof(double) * 3 * N); nrm = malloc(sizeof(double) * 3 * N);
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < 3; k++) pos[3 * i + k] = urand() * L;
            double z = 2 * urand() - 1, ph = 2 * M_PI * urand(), s = sqrt(1 - z * z);
            nrm[3 * i] = s * cos(ph); nrm[3 * i + 1] = s * sin(ph); nrm[3 * i + 2] = z;
        }
    } else { fprintf(stderr, "usage: %s file in.bin out.bin | gen N phi seed\n", argv[0]); return 1; }

    int *counts = malloc(sizeof(int) * N), nth = 1;
    double t0 = 0, t1 = 0;
#ifdef _OPENMP
    nth = omp_get_max_threads();
    t0 = omp_get_wtime();
#endif
    count_overlaps(N, L, pos, nrm, counts);
#ifdef _OPENMP
    t1 = omp_get_wtime();
#endif
    long long tot = 0;
    for (int i = 0; i < N; i++) tot += counts[i];
    printf("N = %d  L = %.2f  mean overlaps per flake = %.4f  threads = %d  wall = %.3f s\n", N, L, (double)tot / N, nth, t1 - t0);
    if (!strcmp(argv[1], "file")) { FILE *o = fopen(argv[3], "wb"); fwrite(counts, 4, N, o); fclose(o); }
    return 0;
}

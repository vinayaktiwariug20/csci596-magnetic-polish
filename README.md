# CSCI 596 final project: platelet flakes in magnetic nail polish

Team: Vinayak Tiwari, Shreeya Halwasia, Ishita Joshi (USC CSCI 596, Fall 2026).

## What this is
A particle-level simulation of magnetic platelet flakes in viscous polish, driven by a magnetic field solver and
parallelized on USC Discovery. The question: does giving each flake a full orientation (a platelet, not a rod)
reproduce the compact bright bead that real users get by turning the magnet 90 degrees, which the existing
rod model cannot? A second question is whether flake-flake contact changes the result, which would also give
the parallel code real communication to do.

## What existed before the course, and what is new
- **Before the course (Aug 2026): the baseline simulator.** Magnetic field solver, rod-based flake model with
  viscosity thickening, 3D viewer, tests. Not part of this repository.
  Live: https://vinayaktiwariug20.github.io/magnetic-lacquer/ and code: https://github.com/vinayaktiwariug20/magnetic-lacquer
- **New for the course (this repository, from Oct 2026):** platelet orientation model, contact detection,
  parallelization (OpenMP, MPI), validation against the baseline and nail photos, scaling study.

## Contents
- `PLAN.md`: project plan, roles, timeline, risks, and the results of the first feasibility tests.
- `spike/platelet_spike.py`: Python feasibility test, rod vs platelet under 0 then 90 degree magnet switching.
  Result: `spike/platelet_spike.png` (rod gives a band, platelet gives a compact central patch).
- `spike/contact_spike.py`: overlap counts for randomly placed oriented thin discs at 1 to 10 percent loading.
  Result: `spike/contact_spike_result.json`.
- `spike/fields.json`: magnetic field at every nail point, exported from the baseline simulator.

## Run the spikes
    cd spike
    python platelet_spike.py fields.json 300    # needs numpy, matplotlib
    python contact_spike.py fields.json 300     # needs numpy

## Status
Feasibility only. All constants are order-of-magnitude estimates, and both spikes are toy models. See `PLAN.md`
for what each test did and did not show. The C implementation and parallel versions are the next step.

## C port (Oct 24)

`c/platelet.c` is the platelet orientation model in C (RK4, same protocol as `spike/platelet_spike.py`), with an OpenMP loop over nail texels.

    cd c
    python export_inputs.py ../spike/fields.json inputs.bin
    gcc -O2 -fopenmp -o platelet_omp platelet.c -lm
    OMP_NUM_THREADS=8 ./platelet_omp inputs.bin out.bin
    python compare.py out.bin ../spike/platelet_spike_maps.npz

Check: matches the Python spike to 1e-14 (6305 texels x 8 flakes). Laptop timing, 8 cores: 6.24 s (1 thread), 0.88 s (8 threads), 7.1x. This is the embarrassingly parallel floor; flake-flake contact is the next step.

## Contact count in C (Oct 26)

`c/contact.c` is the disc-disc overlap test from `spike/contact_spike.py` with a cell-list neighbour search and an OpenMP loop over flakes (each flake counts its own partners, so no write races).

    cd c
    gcc -O2 -fopenmp -o contact contact.c -lm
    python contact_ref.py 4000 0.05 random in.bin ref.bin    # Python reference using the spike's own function
    ./contact file in.bin out.bin && python contact_check.py out.bin ref.bin
    OMP_NUM_THREADS=8 ./contact gen 2000000 0.05 1           # large random run

Checks:
- Identical per-flake counts to the Python spike for random and nearly aligned flakes at 1, 5 and 10% volume fraction (4000 flakes each).
- 2 million random flakes at 5%: 1.5715 overlaps per flake vs excluded-volume theory phi*pi*a/h = 1.5708.

Laptop timing, 2 million flakes: 16.1 s (1 thread), 8.6 s (2), 4.8 s (4), 3.1 s (8), 5.3x. Flakes are not sorted by cell yet, so memory access is scattered; sorting is the obvious next optimisation.

## Contact dynamics: upper bound (Oct 26)

`c/dynamics.c`: flakes at fixed random positions (overlap-free start), all in the field of the nail-centre texel, run through the switching protocol (5 s at 0 deg, 5 s at 90 deg, x6). Mode `ideal` ignores contact. Mode `hard` accepts a rotation step only if it creates no overlap (tries full, half, quarter step). Positions never move, so blocked flakes cannot slide away; `hard` is an upper bound on the effect of contact and `ideal` the lower bound.

    cd c
    python center_field.py ../spike/fields.json field_center.txt
    gcc -O2 -fopenmp -o dynamics dynamics.c -lm
    ./dynamics field_center.txt 0.03 20000 hard 1

Mean flatness (n.N)^2 at the nail centre, 20000 flakes (1 = all mirror-like; geometric limit here is 1.000):

| volume fraction | ideal, end | hard, end | full steps blocked |
|---|---|---|---|
| 1% | 1.000 | 0.962 | 6.9% |
| 3% | 1.000 | 0.873 | 22% |
| 10% | 1.000 | 0.593 | 66% |

Checks: identical output on 1 and 8 threads; no overlaps at the end of any hard run; halving dt twice moves 3% from 0.873 to 0.857 and 10% from 0.593 to 0.583; 80000 flakes or another seed change the result by at most 0.004.

Reading: if flakes could not move at all, contact would cost about 13% of the shine at 3% loading and 40% at 10%. The real effect is somewhere between this and zero; letting flakes translate (overdamped contact forces) is what narrows it. Bulk periodic box only, no nail surface or film walls yet.

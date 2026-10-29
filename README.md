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

## Start here (team, Oct 28)
Question now: do magnetic forces between flakes explain the loading effect measured in real polish (alignment
degrades above about 5 vol%, Pereira et al. 2017), which contact alone does not reproduce? Roles are in `PLAN.md`.
The dated sections below are the log of results so far; the last few matter most.

Build and run the main model (gcc with OpenMP; on Windows use Git Bash or MSYS2, on Discovery `module load gcc`):

    cd c
    gcc -O2 -fopenmp -Wall soft.c -o soft -lm
    python center_field.py ../spike/fields.json field_3152.txt 300 3152    # field at one nail point (centre)
    ./soft field_3152.txt 0.03 5000 free 1 12             # 3% loading, 5000 flakes, contact, 12a film
    ./soft field_3152.txt 0.10 5000 free 1 12 4 5 10      # + film field (chi_in 4), static 5 s hold, brushed start

Arguments: field file, volume fraction, number of flakes, `ideal` (no contact) / `fixed` / `free` (sliding
contact), seed, film thickness in flake radii (0 = periodic cube), chi_in (film field, 0 = off), static hold in s
(0 = 0/90 switching), start spread in degrees (-1 = random start). Lengths are in flake radii a = 6 um.

Where things live in `c/soft.c`:
- `rk4` / `rhs`: magnetic torque on one flake (dn/dt = -kappa (B.n)(B - (B.n) n)). A field from other flakes
  would be added to B here, per flake.
- `build_cells`, `nb_cell`, `cell_of`: cell lists (cell side 2a + H). A dipole cutoff of a few a needs larger cells
  or a second list; the same neighbour loop is what an MPI version splits across ranks (ghost layers, as in
  assignments 4 and 5).
- `sweep`, `apply_push`, `mobility`: contact resolution (projection, Jacobi sweeps). Flake mobilities in
  `mobility` also convert a magnetic force or torque into motion.
- Each flake updates only itself from neighbours' old state, so OpenMP needs no locks and results do not
  depend on thread count; keep that property.

For the visualisation: `soft` currently prints summaries only. Writing flake positions and normals to a file
(for rendering a BRDF lobe per nail point) is a small addition; ask Vinayak or open an issue.

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

## Contact with sliding flakes (Oct 27)

`c/soft.c`: same setup as `dynamics.c`, but flakes may also translate. Each flake is a "spherodisc" (points within H/2 of a zero-thickness disc). Each step: free magnetic rotation (RK4), then Jacobi projection sweeps that push overlapping pairs apart along the closest-point direction, split between the two flakes by viscous disc mobilities (16 eta a face-on, 32/3 eta a edge-on, 32/3 eta a^3 rotation). Flakes that tunnel through a neighbour within one step return to their start-of-step state. No stiffness parameter; dt = 0.01 s.

    cd c
    gcc -O2 -fopenmp -o soft soft.c -lm
    ./soft field_center.txt 0.03 5000 free 1      # modes: ideal | fixed | free

Mean flatness at the nail centre at the end of the protocol (5000 flakes; no contact gives 1.000):

| volume fraction | rejection, frozen (dynamics.c) | projection, frozen | projection, sliding |
|---|---|---|---|
| 3% | 0.880 | 0.913 | 0.998 |
| 10% | 0.593 | 0.637 | 0.993 |

Checks: `ideal` mode at dt = 0.05 reproduces dynamics.c exactly (0.499 after hold 1, 1.000 at end); halving dt moves the 3% sliding result from 0.998 to 1.000. Sliding runs converge (3%: every step; 10%: 49 of 6000 steps leave a residual penetration of at most 0.05 H); 0.5% (3%) and 1.6% (10%) of flake-steps are reverted. Mean displacement over the run: 0.13 a (3%), 0.44 a (10%). Frozen projection never fully converges (rotation alone cannot separate some pairs), so it is only a rough cross-check of the rejection bound.

Reading: in the bulk, contact barely changes the final alignment even at 10% loading, because flakes slide a fraction of a radius out of each other's way. The large loss in the frozen bound comes from pinning, not from contact as such. In overdamped dynamics viscosity sets the time scale of both magnetic and contact motion, so it does not change this comparison, only how far flakes get within a hold. Open: the nail surface (a wall that blocks sliding), larger boxes, and performance (10%: 600 s for 5000 flakes on 8 threads, dominated by the disc-distance search).

Box-size check (Oct 27): 10% sliding with 20000 flakes and another seed gives 0.478 after hold 1 and 0.993 at the end (5000 flakes: 0.479 and 0.993), mean shift 0.44 a, 1.5% reverted flake-steps; 1384 s on 8 threads.

## Nail surface and film top (Oct 27)

`soft.c` takes an optional film thickness (in flake radii): the box is then periodic sideways and bounded by a wall at the nail (z = 0) and one at the film top (z = Lz). A flake keeps its lowest (highest) rim point H/2 from a wall; walls join the projection sweeps but do not move. Without the argument the cube runs reproduce the previous results exactly (3%: 0.493 after hold 1, 0.998 at end, same crossing count).

    ./soft field_center.txt 0.03 5000 free 1 12      # 12a = 72 um film

Flatness by height above the nail, 12a film, 5000 sliding flakes, nail centre:

| loading | after hold 1: 0-2a / middle / 10-12a | end: all bands | end: mean |
|---|---|---|---|
| 0.1% (walls only) | 0.558 / 0.49-0.50 / 0.546 | 0.999-1.000 | 1.000 |
| 3% | 0.564 / 0.49 / 0.560 | 0.995-1.000 | 0.998 |
| 10% | 0.609 / 0.46-0.52 / 0.594 | 0.988-0.997 | 0.994 |

Reading: the walls do not stop alignment at the nail centre. Steeply tilted flakes cannot sit within one radius of a wall, so the wall layers hold fewer, flatter flakes (0.1%: 549 flakes in 0-2a vs about 950 per middle band) and start ahead; by the end every band is flat. This is expected here because the target orientation at the centre (n along B0 x B90) is parallel to the nail, which a wall never blocks. Off-centre texels, where the target is tilted, are the case where a wall could matter; that needs the wall tilted relative to the box (the code now refuses a nail normal that is not along z).

## Across the nail (Oct 27)

`soft.c` now rotates the fields so the local nail normal is z (Rodrigues; the stored normal is renormalised first), so wall runs work at any texel. `c/scan_line.py` runs, per texel: no contact, walls only (0.1%), and 3% and 10% with walls (12a film, 5000 sliding flakes); results in `c/scan_x.csv`. Seven texels from the centre (3152) to where the target orientation stands on edge (3128); the map is symmetric about the centre.

    python scan_line.py scan_x.csv 3152 3148 3144 3140 3136 3132 3128

| texel | geometric limit | end: none / walls / 3% / 10% | after hold 1: none / 10% |
|---|---|---|---|
| 3152 (centre) | 1.000 | 1.000 / 1.000 / 0.998 / 0.994 | 0.499 / 0.521 |
| 3148 | 0.919 | 0.919 / 0.919 / 0.917 / 0.916 | 0.457 / 0.484 |
| 3144 | 0.708 | 0.708 / 0.708 / 0.706 / 0.707 | 0.351 / 0.381 |
| 3140 | 0.442 | 0.442 / 0.442 / 0.441 / 0.442 | 0.219 / 0.244 |
| 3136 | 0.203 | 0.203 / 0.203 / 0.203 / 0.203 | 0.100 / 0.122 |
| 3132 | 0.048 | 0.048 / 0.048 / 0.049 / 0.049 | 0.024 / 0.045 |
| 3128 | 0.000 | 0.000 / 0.000 / 0.001 / 0.001 | 0.000 / 0.021 |

Reading: after six switching cycles every texel sits on its geometric limit (n along B0 x B90) to within 0.006, with or without contact and walls, nail band included. So in this model the converged bead shape is set by field geometry alone. (Corrected Oct 28: an earlier version said contact raises flatness after the first hold by 0.02 to 0.03. That compared runs with different starting states: inserting flakes between walls makes them start flatter at higher loading. With the same start, the first-hold difference vanishes as dt goes to 0; see "Static magnet" below.)

## Flakes' own magnetic field: estimate (Oct 27)

`c/dipole_estimate.py` (no dynamics). Linear flakes, chi_in 4 / chi_perp 1 (model dchi = 3) and a stronger bracket chi_in 8.

1. Film field. The magnetised layer is a thin sheet, so inside it the field component along the nail normal drops by 1/(1 + chi_NN). This tilts B0 and B90 toward the nail and moves the geometric limit. End flatness, chi_in 4 (8):

| texel | 0% | 1% | 3% | 10% |
|---|---|---|---|---|
| 3152 | 1.000 | 1.000 | 1.000 | 1.000 |
| 3148 | 0.919 | 0.921 | 0.924 (0.925) | 0.934 (0.937) |
| 3144 | 0.708 | 0.715 (0.720) | 0.729 (0.740) | 0.768 (0.790) |
| 3140 | 0.442 | 0.455 (0.465) | 0.479 (0.504) | 0.547 (0.601) |
| 3136 | 0.203 | 0.214 (0.224) | 0.235 (0.263) | 0.303 (0.375) |
| 3132 | 0.048 | 0.052 (0.055) | 0.059 (0.071) | 0.087 (0.129) |

2. Neighbour field. Non-overlapping flakes at the centre texel, each a rim of magnetic charge, cutoff 6a. Typical normal field from neighbours relative to |B| (median / 90th percentile): 1% 0.01-0.02 / 0.05-0.09, 3% 0.03-0.06 / 0.08-0.15, 10% 0.05-0.11 / 0.14-0.28. Flatness loss if each flake tilts by that much: at most 0.007 (chi_in 4) and 0.025 (chi_in 8) at 10%.

Reading: the largest effect is the film field, a first-order shift of the bead off-centre at 3-10% loading (0.20 to 0.30-0.38 at texel 3136). It needs no pair sums: it is a correction to the field. Neighbour fields tilt flakes by a few degrees, comparable to the contact transient. Not yet checked: dipole forces moving flakes into chains. Rough drift at 10% is about one radius over the full 60 s protocol, independent of viscosity (rotation and drift both scale as 1/eta).

## Film field in the dynamics (Oct 27)

`c/soft.c` takes an optional last argument `chi_in` (chi_perp = 1): at every step the field component along the nail normal is cut to B_z / (1 + chi_zz), with chi_zz = phi <chi_in (1 - nz^2) + nz^2> from the current flake orientations. Example: `./soft field_3136.txt 0.10 5000 free 1 12 4`.

- No contact (`ideal`, 2000 flakes, film 12a): end flatness at all 7 scan texels, phi 1/3/10%, chi_in 4 and 8, equals the self-consistent estimate of `dipole_estimate.py` to three decimals (e.g. texel 3136: 0.214/0.235/0.303 for chi_in 4, 0.224/0.263/0.375 for chi_in 8).
- With contact (`free`, 5000 flakes, texel 3136, chi_in 4): 3% end 0.235, 10% end 0.304. Contact again changes nothing at the end.

Reading: in this model, loading moves the bead through the film field alone; contact does not add to it. Literature notes are in PLAN.md.

## Static magnet and a time-step artifact (Oct 28)

`soft` arguments 8 and 9: one static hold in B0 (`static_hold_s`) and a start near the nail plane (`start_spread_deg`, -1 = random). Output: rms angle by which flake planes miss B0. Centre texel, film 12a, 5000 flakes, 5 s, contact (`free`), dt 0.01:

| loading | 1% | 3% | 5% | 8% | 10% | 13% |
|---|---|---|---|---|---|---|
| brushed start (10 deg), rms miss | 0.26 | 0.47 | 0.61 | 0.74 | 0.76 | 0.77 deg |
| random start, rms miss | 4.2 | 6.1 | 9.5 | 11.2 | 11.4 | 13.4 deg |

The random-start misalignment does not relax with longer holds (5, 15, 30 s identical), but it is a numerical artifact: about 7% of flakes per step rotate into a neighbour in one step, are sent back by the crossing backstop and never feel a contact push. At 10%: dt 0.01 / 0.005 / 0.0025 gives reverted flake-steps 7.1 / 3.0 / 1.5%, rms miss 11.4 / 7.6 / 5.7 deg, mean miss 2.78 / 1.24 / 0.71 deg (heading to 0). First-hold flatness with contact 0.521 / 0.530 / 0.533 against 0.533 without contact.

Reading: with a converged time step, sliding contact changes neither the first hold nor the end state. The switching results above were already dt-converged at the end (0.998 vs 1.000), so they stand. Contact therefore cannot explain the loading effect Pereira et al. measured, which leaves interactions between flakes (magnetic forces, chaining) as the candidate. For anything near contact, use dt 0.0025 or smaller, or add a per-step rotation limit (the rim of a flake should move less than about H/4 per step); this applies to the flake-force work too.

## Rotation limit with contact (Oct 28)

With contact, `soft` now splits each step into substeps so that no flake turns by more than MAXROT = 0.05 H/a (the rim moves at most H/20) per substep; the bound uses each flake's current magnetic rotation rate. Without contact nothing changes. Override at compile time with `-DMAXROT=(0.1*H/A)`.

Convergence at 10%, static 5 s, random start (rms / mean angle by which flake planes miss B0): no limit 11.4 / 2.78 deg; H/4 4.7 / 0.51; H/10 3.0 / 0.25; H/20 1.95 / 0.14. Both fall roughly in proportion to the limit, so what remains is still numerical and heads to 0.

Static sweep with the limit (centre, film 12a, 5000 flakes, 5 s, random start), rms miss: 1% 0.6, 3% 0.03, 5% 1.5, 8% 1.3, 10% 1.9, 13% 1.9 deg (mean at most 0.19 deg). First-hold flatness 0.51 to 0.54, as without contact.

Switching protocol check (3%, centre, film 12a): end 1.000, no crossings, 60 s on 8 threads (41 s before).

Reading: with sliding contact and a converged step, contact does not keep flakes from aligning in a static field at any loading up to 13%. The loading effect measured by Pereira et al. is not reproduced by contact.

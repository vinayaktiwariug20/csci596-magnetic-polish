# CSCI 596 final project plan (draft)

Team: Vinayak Tiwari, Shreeya Halwasia, Ishita Joshi (roles below are drafts; edit freely).

## One-line project
Particle-level simulation of magnetic platelet flakes in viscous polish, driven by the existing field solver,
parallelized and benchmarked on USC Discovery, to test whether a full-orientation (platelet) flake model
reproduces effects the current independent-rod model cannot.

## Background: what exists vs. what is new
- **Existing (built before the course began, Aug 7-20, 2026):** exact magnetic field solver, rod-ensemble finish
  model with a clock (viscosity thickening), shader, 3D viewer, tests, README.
- **New for CSCI 596 (this folder, from Oct 2026):** per-flake platelet orientation simulation, parallelization,
  validation against the existing model, scaling study.
- This split must be stated in the final report and told to the professor.

## Research question
The current model treats every flake as a nematic rod with one axis, so it cannot distinguish mono-axial from
bi-axial orientation, and it does not reproduce the compact bright bead from alternating perpendicular holds
(it prints a bar, or a cross with chaining). Does giving each flake a full orientation frame (a platelet with an
easy plane) under a time-varying field reproduce that bead?

Secondary tests (cheap, no new code):
1. Set gel `eta0` to a few hundred Pa.s and try the spinning-horseshoe bead. (Gel is currently defined with
   eta0 = 0.6 Pa.s, which makes lag impossible.)
2. Hybrid velvet / cat eye: rotate the bridging bar ~45 deg and flip its polarity; check the diagonal and the
   view-dependence reported by real users.

## Levels
- **Level 1 (floor):** single-flake platelet dynamics in C: torque from the field, viscous rotational drag,
  time-stepping; compare with the rod model for a static field.
- **Level 2:** ensemble of flakes per texel under the time-varying field (perpendicular switching); compare the
  predicted bright patch with the rod model and with photographs.
- **Level 3 (stretch):** gel rheology (yield stress), flake-flake interactions, 3D.

## Roles (draft)
| Role | Scope |
|---|---|
| A. Simulation core | Platelet orientation dynamics in C |
| B. Parallelization & benchmarks | OpenMP/MPI on Discovery, Slurm scripts, scaling plot |
| C. Validation & viewer | Compare with existing finish model and photos, feed results into viewer, README |

Each role must be written specifically in the final report (who did what).

## Timeline
| When | Goal |
|---|---|
| by Oct 21 | Decision made, roles written down, professor meeting scheduled |
| by Nov 3 | Level 1 working |
| Nov 4-17 | Parallel version running and validated (OpenMP then MPI) on Discovery |
| Nov 18-29 | Scaling plot, viewer integration, README, pitch; professor talk done by Nov 20 |
| Nov 30 / Dec 2 | One-minute pitch |
| Dec 11 | Repository due; every member submits the same URL on Brightspace |

## Draft one-minute pitch
Magnetic nail polish gets its bright line because flakes follow the magnet's field and together act as a mirror.
Current previews treat each flake as an independent rod, so they cannot reproduce effects that depend on the
flake's second axis, such as the bead made by alternating perpendicular holds. We simulate platelet flakes with
full orientation in viscous polish, parallelized with OpenMP and MPI on USC's Discovery cluster, and test
whether that reproduces the missing effect. This moves magnetic-finish previews from trial and error toward
physically grounded prediction.

## Risks
- Platelet model loses the closed-form stable integrator, so use small time steps and test stability early.
- Polish constants are estimates, not measured; state this in the report.
- Scope creep: Level 1 is the floor.

## References to read first
- Rosenberg, Kantorovich, Ivanov, Camp, Nanoscale 2024 (magnetic nanoplatelets, Brownian dynamics).
- Patents on orienting magnetic flakes (mono-axial vs bi-axial orientation), e.g. US 9027479.
- Course notes: FinalNote.pdf (story, Neno's ten questions, team roles).

## Spike result (Oct 14, 2026): platelet vs rod under 0/90 deg switching
Script: `spike/platelet_spike.py`, fields exported from the sandbox solver (`spike/fields.json`), output `spike/platelet_spike.png`.
Setup: cat-eye bar held 5 s at 0 deg then 5 s at 90 deg, 6 cycles, gel 300 Pa.s, 8 flakes/texel, log-normal rate spread.
- Rod (current model): single seam, then a broad horizontal band along the LAST direction. No central patch.
- Platelet (easy plane, two axes): after switching, a compact round bright patch at the centre where the two seams cross
  (flat area ~1.3% of the nail, aspect 1.16). It matches the geometric limit n ~ B0 x B90.
- Caveats: toy model (ideal easy plane, no Brownian motion, independent flakes, random initial orientation,
  estimated constants); the "flatness" proxy is not the true optical response; one thin bright line in the rod map
  is unexplained. The patch POSITION does not depend on viscosity; its formation TIME does.
- Next: vary cycles, hold time, viscosity, polydispersity; compare patch size with photographs of real nails.

## Contact spike (Oct 19, 2026): do platelet-platelet overlaps matter?
Script: `spike/contact_spike.py` (thin discs, aspect ratio 10, ideal-gas positions, orientations from the platelet model).
Overlap test checked against the excluded-volume formula for random thin discs: 1.60 vs 1.57 overlaps per flake at 5% loading.
Mean overlapping neighbours per flake, by volume fraction:
| stage | 1% | 3% | 5% | 10% |
|---|---|---|---|---|
| random (at rest) | 0.33 | 0.96 | 1.60 | 3.09 |
| after hold 1 (aligned) | 0.27 | 0.75 | 1.29 | 2.56 |
| mid-switch to 90 deg | 0.07 | 0.25 | 0.39 | 0.83 |
| after hold 2 (aligned) | 0.10 | 0.23 | 0.40 | 0.80 |
- Contact cannot be ignored above a few percent loading. Real loading is unknown (not in the repo), so this is a sweep.
- Overlap is NOT higher during the switch than at rest; it falls as flakes align. No evidence yet that contact changes the bead.
- Positions are uncorrelated with orientation (no forces, no exclusion), so these are would-be overlaps, not a dynamics result.

### Real loading: what the search found (Oct 19, 2026)
- No magnetic-polish-specific loading found. No SDS for a cat-eye polish turned up (only raw iron oxide powders).
- WO2006037907A1 (orienting magnetic bodies in cosmetics, nail polish named): magnetic substances ~0.05-50 wt% in the composition, preferred 1-30 wt%. Platelet = largest dimension / thickness >= 5. Size D50 0.5-300 um preferred. No film thickness.
- Kobo K-RAY Mag pigments (SEP-Magnetic-001, Mar 2026): mica/TiO2 or iron-oxide pearls with black iron oxide, particle size 10-60 um. No use level in the PDF.
- US5470562: pigment PASTE 5-40% of a nail polish; paste itself ~26-41 wt% pigment. Not magnetic specific; no particle or film thickness.
- Converted (assumes flake density ~3 g/cm3, fluid ~1; assumptions, not sourced): 1 wt% -> 0.34 vol%, 5 -> 1.7, 10 -> 3.6, 30 -> 12.5.
- Sandbox flakeRadius = 6 um (12 um diameter) is at the low end of the 10-60 um pigment size.
- Still needed: measure a real bottle (dry-down weight, magnetic separation, or microscope count).

### Loading plan (decided Oct 19, 2026): sweep and fit, no bottle measurement
No microscope or precision scale available, and phone photos cannot resolve 10-60 um flakes reliably.
- Treat flake volume fraction as an unknown parameter. Sweep 1%, 3%, 10% (plausible wet range ~2-4%, from the patent weight-percent range converted with assumed densities).
- Report results as conditional: contact matters if loading is above about 3%.
- Fit loading to nail photos in the validation step (role C); report it as a fitted value with an error bar.
- Optional: email a brand for pigment weight percent or an SDS; ask the professor about lab access.

## Professor meeting (Oct 19, 2026) and debris baseline research (Oct 20, 2026)
Meeting: he liked both the debris idea and the polish project, polish more. Asked for details, the simulator link, and wrote down the team.
He said we may build on an existing baseline rather than do orbital from scratch. He has not yet answered whether independent flakes are enough parallel content.

Open baseline pieces for debris (no full open Kessler cascade simulator found):
- Fragmentation: Kessler.rs (Rust, MIT) and KessPy (Python wrapper), NASA standard breakup model; early-stage, accuracy not guaranteed before 1.0. (An unrelated "Kessler" ML library exists at Surrey/Oxford.)
- Orbit motion: SGP4 implementations (python-sgp4 etc.), hapsira (fork of poliastro, archived Oct 2023), Orekit (Java). Analytic near-circular orbits may be enough.
- Collision detection with MPI/OpenMP: REBOUND (tree and plane-sweep; MPI essential tree). 2012 paper; current API not checked. Risk: its MPI machinery could make the parallel contribution someone else's.
- Full models: NASA LEGEND (Fortran 90), ESA MASTER (free portal access). Open-source status not confirmed.
To check if debris is chosen: KessPy output vs the published breakup formula; REBOUND docs for orbits around a central body; whether SGP4 is needed.

### Verification of the debris baseline pieces (Oct 20, 2026), corrects the entry above
- KessPy 0.2.0: NASA Standard Breakup Model, Python wrapper over Kessler.rs. Last release Dec 6, 2023 (stale). Accuracy not guaranteed before 1.0. MIT stated for Kessler; kesspy's own license not stated on its page.
- hapsira 0.18.0: MIT, last release Dec 24, 2023 (NOT actively maintained, earlier note said "maintained").
- sgp4 2.27: MIT, released Jul 2026, wraps Vallado's C++ code, vectorized (SatrecArray) when the C++ build is installed. Needs TLE/OMM input data. Active.
- REBOUND 5.2.2 (Oct 5, 2026): GPL-3.0 (page inconsistent on "only" vs "or later"). Collision detection routines and OpenMP yes; MPI "supported for some special use cases only". Earlier note on full MPI collision detection came from the 2012 paper and is outdated. GPL matters if we copy code into our repo.

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

Reading: after six switching cycles every texel sits on its geometric limit (n along B0 x B90) to within 0.006, with or without contact and walls, nail band included. So in this model the converged bead shape is set by field geometry alone. Contact shows up only in the transient: after the first hold, 10% loading raises flatness by 0.02 to 0.03 at every texel (crowding pushes flakes towards a common plane). If contact matters for real nails, it is in short or single holds, not after full switching.

## Literature (Oct 27)

- Pereira, Paes Leme, Marschner, Rusinkiewicz, "Printing Anisotropic Appearance with Magnetic Flakes", ACM Trans. Graph. 36(4) 123 (SIGGRAPH 2017). https://gfx.cs.princeton.edu/pubs/Pereira_2017_PAA/flakes.pdf
  Closest to us: magnetic nail polish and UV gel with commercial flakes (mica coated with iron oxide and TiO2, 5-25 and 10-60 um). Time-varying fields (circular, spiral, line) shape the orientation distribution; measured BRDF lobes. Loading test at 1, 6, 12, 18, 30 wt% under a bar magnet: above 12 wt% the lobe widens, "they begin to physically interfere and can no longer be oriented parallel". At a density of roughly 3 g/cm3 for flakes and 1.1 for gel, 12 wt% is about 5 vol% and 30 wt% about 13 vol%. A measured loading effect we can try to reproduce (static field, orientation spread vs loading).
- Le Ferrand and Arrieta, "Magnetically driven in-plane modulation of the 3D orientation of vertical ferromagnetic flakes", Soft Matter (2022), arXiv:2201.05097. Ni flakes 15.8 um x 0.58 um (aspect 27, close to our a = 6 um, h = 0.6 um), 1 vol% in epoxy (13 Pa s) or PDMS (3.5 Pa s), rotating magnets, COMSOL fields predict in-plane orientation. Fields kept below 3 mT to avoid magnetophoresis toward the magnet; even so, ferromagnetic flakes formed long chains along the field under rotating fields (not seen under static fields in the viscous, curing epoxy). Paramagnetic flakes are reported not to chain.
- Erb, Segmehl, Charilaou, Loeffler, Studart, "Non-linear alignment dynamics in suspensions of platelets under rotating magnetic fields", Soft Matter 8, 7604 (2012). Platelets decorated with magnetic nanoparticles: synchronous rotation below a critical frequency, biaxial alignment in the plane of rotation above it. Our 0/90 switching is a two-step version of this.
- Rosenberg, Kantorovich, Ivanov, Camp, Nanoscale (2024) and PCCP (2023): Brownian dynamics of dipolar nanoplatelets (permanent dipoles, clustering). Different regime (nanometre, thermal), shows platelet clustering differs from spheres.
- Not found: any physics paper or simulation of magnetic nail polish itself; no measured magnetisation of cosmetic flakes.

Consequences for the model: (1) cosmetic flakes are mica coated with iron oxide, so their susceptibility is lower than nickel's; chi_in 4 is a guess, keep the bracket. (2) Magnetophoresis (field gradients near the magnet pull flakes) is a further effect we leave out. (3) Chaining is observed for ferromagnetic flakes; whether iron-oxide-coated flakes chain at 160 mT is open.

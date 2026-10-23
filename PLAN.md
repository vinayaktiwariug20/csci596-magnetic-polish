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

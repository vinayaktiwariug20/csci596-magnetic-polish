"""Write a VMD script that draws a soft.c flake state as 3D discs.

Each flake is a disc of radius a (= 1) and thickness 0.1a at its position,
facing along its normal, coloured by its tilt from the nail normal
(blue = lying flat, red = 45 deg or more).  The nail (z = 0) is a dark plate
and the film top (z = Lz) is outlined.

Usage:
    python visualization/flakes_vmd.py state.csv flakes.tcl --image flakes.tga [--batch]
    vmd -e flakes.tcl                                   (interactive: rotate it yourself)
    vmd -dispdev text -e flakes.tcl                     (with --batch: render and quit)
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "c"))
from read_flakes import load  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("state")
    ap.add_argument("tcl")
    ap.add_argument("--image", default="flakes.tga")
    ap.add_argument("--title", default="")
    ap.add_argument("--batch", action="store_true")
    args = ap.parse_args()

    head, d = load(args.state)
    box = head["box"]
    x, nrm = d[:, :3], d[:, 3:6]
    nrm = nrm * np.sign(nrm[:, 2:3] + 1e-12)
    tilt = np.degrees(np.arccos(np.clip(nrm[:, 2], 0, 1)))
    half = 0.05 * nrm

    L = ["display projection Orthographic", "display depthcue off", "axes location Off",
         "color Display Background gray", "color scale method BWR",
         "display resize 1200 800", "mol new", "graphics top materials on",
         "graphics top material Diffuse",
         "display shadows off"]
    # nail plate and film-top outline
    L += ["graphics top color black",
          f"graphics top triangle {{0 0 -0.3}} {{{box[0]} 0 -0.3}} {{{box[0]} {box[1]} -0.3}}",
          f"graphics top triangle {{0 0 -0.3}} {{{box[0]} {box[1]} -0.3}} {{0 {box[1]} -0.3}}",
          "graphics top color black"]
    corners = [(0, 0), (box[0], 0), (box[0], box[1]), (0, box[1])]
    for z in (0, box[2]):
        for (a, b), (c, e) in zip(corners, corners[1:] + corners[:1]):
            L.append(f"graphics top line {{{a} {b} {z}}} {{{c} {e} {z}}} width 2 style dashed")
    for a, b in corners:
        L.append(f"graphics top line {{{a} {b} 0}} {{{a} {b} {box[2]}}} width 2 style dashed")
    for p, h, t in zip(x, half, tilt):
        cid = 33 + int(min(t, 45) / 45 * 1023)
        L.append(f"graphics top color {cid}")
        L.append("graphics top cylinder {%.3f %.3f %.3f} {%.3f %.3f %.3f} radius 1.0 resolution 14 filled yes"
                 % (*(p - h), *(p + h)))
    L += ["display resetview", "rotate z by -25", "rotate x by -58", "scale by 1.2",
          f"render TachyonInternal {args.image}"]
    if args.batch:
        L.append("quit")
    with open(args.tcl, "w") as f:
        f.write("\n".join(L) + "\n")
    print(f"wrote {args.tcl}: {len(x)} flakes, median tilt {np.median(tilt):.1f} deg")


if __name__ == "__main__":
    main()

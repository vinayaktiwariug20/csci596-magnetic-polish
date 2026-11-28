"""Render a reflection proxy from simulated platelet orientations.

Flakes start brushed (default): normal along the local nail normal with an
rms tilt of --start-deg (10 deg), as in the Discovery runs (soft.c, last
argument).  --start random gives the earlier isotropic start.
"""

import argparse
import json

import matplotlib.pyplot as plt
import numpy as np


def unit(v):
    return v / np.maximum(
        np.linalg.norm(v, axis=-1, keepdims=True),
        1e-30,
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("fields")
    parser.add_argument("--eta", type=float, default=300.0)
    parser.add_argument("--flakes", type=int, default=32)
    parser.add_argument("--width-deg", type=float, default=8.0)
    parser.add_argument("--start", choices=["brushed", "random"], default="brushed")
    parser.add_argument("--start-deg", type=float, default=10.0,
                        help="rms tilt of the brushed start (deg)")
    parser.add_argument("--output", default="brdf_render.png")

    args = parser.parse_args()

    with open(args.fields) as f:
        data = json.load(f)

    count = data["count"]

    p = np.array(data["P"]).reshape(count, 3)
    b0 = np.array(data["B0"]).reshape(count, 3)
    b90 = np.array(data["B90"]).reshape(count, 3)
    nail_normal = unit(np.array(data["N"]).reshape(count, 3))

    mu0 = 4e-7 * np.pi
    a = 6e-6
    h = a / 10.0
    dchi = 3.0

    kappa0 = (
        3 * np.pi * h * dchi
        / (32 * mu0 * args.eta * a)
    )

    rng = np.random.default_rng(1)

    spread = 0.55

    kappa = kappa0 * np.exp(
        spread
        * rng.standard_normal(
            (count, args.flakes)
        )
    )

    def rhs(x, field):
        bx = np.einsum(
            "tfk,tk->tf",
            x,
            field,
        )

        return (
            -kappa[..., None]
            * bx[..., None]
            * (
                field[:, None, :]
                - bx[..., None] * x
            )
        )

    def evolve(x, field, duration, dt=0.05):
        steps = int(round(duration / dt))

        for _ in range(steps):
            k1 = rhs(x, field)

            k2 = rhs(
                unit(x + 0.5 * dt * k1),
                field,
            )

            k3 = rhs(
                unit(x + 0.5 * dt * k2),
                field,
            )

            k4 = rhs(
                unit(x + dt * k3),
                field,
            )

            x = unit(
                x
                + dt
                / 6
                * (
                    k1
                    + 2 * k2
                    + 2 * k3
                    + k4
                )
            )

        return x

    noise = rng.standard_normal((count, args.flakes, 3))

    if args.start == "random":
        x0 = unit(noise)
    else:
        # same construction as soft.c: N + (spread / sqrt 2) * gaussian
        x0 = unit(
            nail_normal[:, None, :]
            + np.radians(args.start_deg) / np.sqrt(2) * noise
        )

    static = evolve(
        x0.copy(),
        b0,
        45.0,
    )

    switching = x0.copy()

    for _ in range(6):
        switching = evolve(
            switching,
            b0,
            5.0,
        )

        switching = evolve(
            switching,
            b90,
            5.0,
        )

    half_vector = np.array(
        [0.0, 0.0, 1.0]
    )

    width = np.radians(
        args.width_deg
    )

    def reflection(x):
        cos_angle = np.abs(
            np.einsum(
                "tfk,k->tf",
                x,
                half_vector,
            )
        )

        angle = np.arccos(
            np.clip(
                cos_angle,
                0.0,
                1.0,
            )
        )

        # Laplacian-like angular reflection lobe.
        contribution = np.exp(
            -angle / width
        )

        return contribution.mean(axis=1)

    static_intensity = reflection(static)
    switching_intensity = reflection(switching)

    scale = max(
        static_intensity.max(),
        switching_intensity.max(),
    )

    static_intensity /= scale
    switching_intensity /= scale

    print(
        "static:    max %.3f  mean %.3f"
        % (
            static_intensity.max(),
            static_intensity.mean(),
        )
    )

    print(
        "switching: max %.3f  mean %.3f"
        % (
            switching_intensity.max(),
            switching_intensity.mean(),
        )
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(9, 6),
        constrained_layout=True,
    )

    panels = [
        (
            axes[0],
            static_intensity,
            "Static field, 45 s",
        ),
        (
            axes[1],
            switching_intensity,
            "0°/90° switching, 60 s",
        ),
    ]

    scatter = None

    for ax, intensity, title in panels:
        scatter = ax.scatter(
            p[:, 0],
            p[:, 1],
            c=intensity,
            s=18,
            marker="s",
            cmap="magma",
            vmin=0,
            vmax=1,
        )

        ax.set_title(title)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])

    fig.colorbar(
        scatter,
        ax=axes,
        label="Normalized reflected intensity",
        shrink=0.75,
    )

    start_label = (
        "random start" if args.start == "random"
        else f"brushed start ({args.start_deg:g}° rms)"
    )

    fig.suptitle(
        f"Platelet reflection proxy\nIndependent flakes, {start_label}, no flake-flake forces, {args.width_deg:g}° width",
        fontsize=14,
    )

    fig.savefig(
        args.output,
        dpi=220,
        bbox_inches="tight",
    )

    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()

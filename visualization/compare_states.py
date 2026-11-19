"""Compare static and switching flake states."""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "c"))

from read_flakes import load


RLINK = 2.2
ANGLE = 37.0


def oriented_normals(d):
    n = d[:, 3:6].copy()
    n /= np.linalg.norm(n, axis=1)[:, None]

    n[n[:, 2] < 0] *= -1
    return n


def orientation_stats(d):
    n = oriented_normals(d)

    flatness = np.mean(n[:, 2] ** 2)
    tilt = np.degrees(
        np.arccos(np.clip(n[:, 2], -1.0, 1.0))
    )

    return n, flatness, np.median(tilt), tilt


def chain_data(head, d, field_key):
    box = head["box"]

    b = head[field_key]
    b = b / np.linalg.norm(b)

    x = d[:, :3]

    tree = cKDTree(
        np.mod(
            x,
            [box[0], box[1], box[2] + 2 * RLINK],
        ),
        boxsize=[
            box[0],
            box[1],
            box[2] + 2 * RLINK,
        ],
    )

    pairs = tree.query_pairs(
        RLINK,
        output_type="ndarray",
    )

    if len(pairs) == 0:
        return pairs, (0.0, 1.0, 0.0, 1)

    dv = x[pairs[:, 1]] - x[pairs[:, 0]]

    dv[:, :2] -= box[:2] * np.round(
        dv[:, :2] / box[:2]
    )

    cosb = np.abs(dv @ b) / np.linalg.norm(
        dv,
        axis=1,
    )

    linked = pairs[
        cosb > np.cos(np.radians(ANGLE))
    ]

    n = len(x)

    graph = coo_matrix(
        (
            np.ones(len(linked)),
            (linked[:, 0], linked[:, 1]),
        ),
        shape=(n, n),
    )

    _, labels = connected_components(
        graph,
        directed=False,
    )

    sizes = np.bincount(labels)
    per_flake = sizes[labels]

    stats = (
        len(linked) / n,
        per_flake.mean(),
        np.mean(per_flake >= 3),
        sizes.max(),
    )

    return linked, stats


def visible_segments(head, d, linked):
    box = head["box"]
    x = d[:, :3]

    segments = []

    for i, j in linked:
        raw = x[j, :2] - x[i, :2]

        minimum = raw - box[:2] * np.round(
            raw / box[:2]
        )

        if not np.allclose(raw, minimum):
            continue

        segments.append(
            [x[i, :2], x[j, :2]]
        )

    return segments


def lambert_projection(normals):
    nz = np.clip(normals[:, 2], -1.0, 1.0)

    scale = np.sqrt(
        2.0 / np.maximum(1.0 + nz, 1e-12)
    )

    return (
        scale * normals[:, 0],
        scale * normals[:, 1],
    )


def plot_state(ax_xy, ax_n, path, label, field_key):
    head, d = load(path)

    normals, flatness, median_tilt, tilt = orientation_stats(d)

    linked, stats = chain_data(
        head,
        d,
        field_key,
    )

    links_per_flake, mean_chain, _, _ = stats

    segments = visible_segments(
        head,
        d,
        linked,
    )

    if segments:
        lines = LineCollection(
            segments,
            linewidths=0.4,
            alpha=0.22,
        )
        ax_xy.add_collection(lines)

    scatter = ax_xy.scatter(
        d[:, 0],
        d[:, 1],
        c=tilt,
        s=7,
        vmin=0,
        vmax=45,
        cmap="viridis",
        rasterized=True,
    )

    box = head["box"]

    ax_xy.set_xlim(0, box[0])
    ax_xy.set_ylim(0, box[1])
    ax_xy.set_aspect("equal")

    ax_xy.set_xlabel("x / a")
    ax_xy.set_ylabel("y / a")

    ax_xy.set_title(
        f"{label}\n"
        f"flatness {flatness:.3f} | "
        f"median tilt {median_tilt:.1f}°\n"
        f"links per flake {links_per_flake:.3f} | "
        f"mean chain {mean_chain:.1f}"
    )

    b = head[field_key]
    bxy = b[:2]

    if np.linalg.norm(bxy) > 0:
        bxy = bxy / np.linalg.norm(bxy)

        origin = np.array(
            [0.15 * box[0], 0.15 * box[1]]
        )

        length = 0.12 * min(box[0], box[1])

        ax_xy.annotate(
            field_key.replace("_T", ""),
            xy=origin + length * bxy,
            xytext=origin,
            arrowprops=dict(
                arrowstyle="->",
                linewidth=1.3,
            ),
        )

    u, v = lambert_projection(normals)

    ax_n.scatter(
        u,
        v,
        c=tilt,
        s=7,
        vmin=0,
        vmax=45,
        cmap="viridis",
        rasterized=True,
    )

    theta = np.linspace(
        0,
        2 * np.pi,
        400,
    )

    radius = np.sqrt(2.0)

    ax_n.plot(
        radius * np.cos(theta),
        radius * np.sin(theta),
        linewidth=1.0,
    )

    for angle in (10, 20, 30):
        ring = 2 * np.sin(
            np.radians(angle) / 2
        )

        ax_n.plot(
            ring * np.cos(theta),
            ring * np.sin(theta),
            linewidth=0.6,
            linestyle=":",
            alpha=0.6,
        )

        ax_n.text(
            ring + 0.03,
            0.07,
            f"{angle}°",
            fontsize=8,
            alpha=0.7,
        )

    ax_n.axhline(
        0,
        linewidth=0.5,
        alpha=0.35,
    )

    ax_n.axvline(
        0,
        linewidth=0.5,
        alpha=0.35,
    )

    ax_n.set_xlim(-1.5, 1.5)
    ax_n.set_ylim(-1.5, 1.5)
    ax_n.set_aspect("equal")

    ax_n.set_xlabel("Lambert x")
    ax_n.set_ylabel("Lambert y")
    ax_n.set_title("Flake normals")

    return scatter


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "static_csv",
    )

    parser.add_argument(
        "switching_csv",
    )

    parser.add_argument(
        "--output",
        default="static_vs_switching.png",
    )

    args = parser.parse_args()

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(12.5, 11.5),
        constrained_layout=True,
    )

    scatter = plot_state(
        axes[0, 0],
        axes[1, 0],
        args.static_csv,
        "Static field, 45 s",
        "B0_T",
    )

    plot_state(
        axes[0, 1],
        axes[1, 1],
        args.switching_csv,
        "0°/90° switching, 60 s",
        "B90_T",
    )

    fig.colorbar(
        scatter,
        ax=axes,
        label="Tilt from nail normal (degrees)",
        shrink=0.78,
    )

    fig.suptitle(
        "Static field vs 0°/90° switching",
        fontsize=15,
    )

    fig.savefig(
        args.output,
        dpi=220,
        bbox_inches="tight",
    )

    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()

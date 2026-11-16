"""Chain statistics along B0 or B90."""

import argparse

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from read_flakes import load


RLINK = 2.2
ANGLE = 37.0


def chains(head, d, field_key):
    box = head["box"]

    b = head[field_key]
    b /= np.linalg.norm(b)

    x = d[:, :3]

    tree = cKDTree(
        np.mod(x, [box[0], box[1], box[2] + 2 * RLINK]),
        boxsize=[box[0], box[1], box[2] + 2 * RLINK],
    )

    pairs = tree.query_pairs(
        RLINK,
        output_type="ndarray",
    )

    if len(pairs) == 0:
        return 0.0, 1.0, 0.0, 1

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

    return (
        len(linked) / n,
        per_flake.mean(),
        np.mean(per_flake >= 3),
        sizes.max(),
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("file")

    parser.add_argument(
        "--field",
        choices=["B0", "B90"],
        default="B0",
    )

    args = parser.parse_args()

    head, d = load(args.file)

    links, mean_chain, fraction_three, longest = chains(
        head,
        d,
        f"{args.field}_T",
    )

    print(f"field: {args.field}")
    print(f"links/flake: {links:.3f}")
    print(f"mean chain (per flake): {mean_chain:.2f}")
    print(f"fraction in chains >=3: {fraction_three:.3f}")
    print(f"longest chain: {longest}")


if __name__ == "__main__":
    main()

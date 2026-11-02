import numpy as np
from numpy.typing import NDArray
from pymatgen.core import Structure, SymmOp

def pair_opposite_vectors(
    directions: NDArray[np.float64],
) -> list[tuple[int, int]]:
    """
    Pair six approximately unit vectors into three opposite pairs.

    The pairing minimizes dot products: opposite vectors have dot ≈ -1.
    """
    remaining = set(range(len(directions)))
    pairs: list[tuple[int, int]] = []

    while remaining:
        i = min(remaining)
        remaining.remove(i)

        if not remaining:
            raise ValueError("Odd number of ligand directions.")

        j = min(
            remaining,
            key=lambda k: np.dot(directions[i], directions[k]),
        )
        remaining.remove(j)
        pairs.append((i, j))

    print("\nOpposite ligand pairs:")
    for i, j in pairs:
        dot = np.dot(directions[i], directions[j])
        angle = np.degrees(np.arccos(np.clip(dot, -1.0, 1.0)))
        print(f"  {i} ↔ {j}: dot = {dot:.10f}, angle = {angle:.6f}°")

    return pairs


def verify_octahedron_orientation(
    structure: Structure,
    metal_index: int,
    oxygen: str = "O",
    neighbor_radius: float = 3.0,
) -> None:
    """Print the six nearest ligand directions after rotation."""
    center = structure[metal_index]

    neighbors = [
        nn
        for nn in structure.get_neighbors(center, neighbor_radius)
        if nn.specie.symbol == oxygen
    ]
    neighbors.sort(key=lambda nn: nn.nn_distance)
    neighbors = neighbors[:6]

    vectors = np.array(
        [nn.coords - center.coords for nn in neighbors],
        dtype=float,
    )
    directions = vectors / np.linalg.norm(vectors, axis=1)[:, None]

    print("\nLigand directions after rotation:")
    for direction in directions:
        print(
            "  "
            + np.array2string(
                direction,
                precision=8,
                suppress_small=True,
            )
        )

"""
tesserae — k-ary derivative structure enumeration via T ⋊ P semidirect product decomposition.

All enumeration functions return a list of colorings, where each coloring is a
list of integer species labels (0-indexed) of length equal to the number of sites
in the supercell.
"""

from typing import List

PointGroup = List[List[List[int]]]
"""3×3 integer rotation matrices forming the parent point group."""

Coloring = List[int]
"""Species labels for each site, 0-indexed."""

def enumerate_supercell(
    supercell_matrix: List[List[int]],
    composition: List[int],
    parent_point_group: PointGroup,
) -> List[Coloring]:
    """Enumerate all symmetry-inequivalent colorings of a supercell.

    Args:
        supercell_matrix: 3×3 integer matrix defining the supercell.
        composition: Number of sites per species, summing to |det(supercell_matrix)|.
        parent_point_group: List of 3×3 integer rotation matrices (the full group,
            not just generators).

    Returns:
        List of colorings, one per inequivalent structure.
    """
    ...

def enumerate_supercell_primitive(
    supercell_matrix: List[List[int]],
    composition: List[int],
    parent_point_group: PointGroup,
) -> List[Coloring]:
    """Like enumerate_supercell but excludes superperiodic structures.

    Superperiodic structures are those representable by a smaller supercell
    (invariant under a non-identity translation). Matches the enumlib convention.
    """
    ...

def enumerate_at_index(
    parent_point_group: PointGroup,
    index: int,
    composition: List[int],
) -> List[Coloring]:
    """Enumerate colorings across all inequivalent supercells at a given volume index.

    Sweeps all Hermite Normal Form matrices with determinant `index`, reduces
    them to symmetry-inequivalent representatives under `parent_point_group`,
    and enumerates colorings with `composition` for each.

    Args:
        parent_point_group: Full point group (list of 3×3 integer matrices).
        index: Volume index = determinant of the supercell matrix.
        composition: Number of sites per species, summing to `index`.
    """
    ...

def enumerate_at_index_primitive(
    parent_point_group: PointGroup,
    index: int,
    composition: List[int],
) -> List[Coloring]:
    """Like enumerate_at_index but excludes superperiodic structures."""
    ...

def inequivalent_hnfs(
    hnfs: List[List[List[int]]],
    parent_point_group: PointGroup,
) -> List[List[List[int]]]:
    """Filter a list of HNF matrices to symmetry-inequivalent representatives.

    Two HNFs are equivalent if one can be obtained from the other by a
    point-group rotation. Returns one representative per equivalence class.

    Args:
        hnfs: List of 3×3 (or 2×2) upper-triangular integer matrices.
        parent_point_group: List of rotation matrices.
    """
    ...

def enumerate_from_cif(
    cif_path: str,
    supercell_matrix: List[List[int]],
    composition: List[int],
    symprec: float = 1e-5,
) -> List[Coloring]:
    """Enumerate colorings of a fixed supercell, reading symmetry from a CIF file.

    Args:
        cif_path: Path to a CIF file with explicit atom sites.
        supercell_matrix: 3×3 integer matrix defining the supercell.
        composition: Number of sites per species, summing to n_sites × |det(supercell_matrix)|.
        symprec: Symmetry tolerance for spglib (default 1e-5).
    """
    ...

def enumerate_from_cif_primitive(
    cif_path: str,
    supercell_matrix: List[List[int]],
    composition: List[int],
    symprec: float = 1e-5,
) -> List[Coloring]:
    """Like enumerate_from_cif but excludes superperiodic structures."""
    ...

def enumerate_from_cif_at_index(
    cif_path: str,
    index: int,
    composition: List[int],
    symprec: float = 1e-5,
) -> List[Coloring]:
    """Enumerate colorings across all inequivalent supercells at a volume index.

    Reads the parent structure's symmetry from the CIF file, then sweeps all
    inequivalent supercell shapes at the given volume index.

    Args:
        cif_path: Path to a CIF file with explicit atom sites.
        index: Volume index = determinant of the supercell matrix.
        composition: Number of sites per species, summing to `index`.
        symprec: Symmetry tolerance for spglib (default 1e-5).
    """
    ...

def enumerate_from_cif_at_index_primitive(
    cif_path: str,
    index: int,
    composition: List[int],
    symprec: float = 1e-5,
) -> List[Coloring]:
    """Like enumerate_from_cif_at_index but excludes superperiodic structures."""
    ...

def enumerate_from_poscar(
    poscar_path: str,
    supercell_matrix: List[List[int]],
    composition: List[int],
    symprec: float = 1e-5,
) -> List[Coloring]:
    """Enumerate colorings of a fixed supercell, reading symmetry from a POSCAR file."""
    ...

def enumerate_from_poscar_primitive(
    poscar_path: str,
    supercell_matrix: List[List[int]],
    composition: List[int],
    symprec: float = 1e-5,
) -> List[Coloring]:
    """Like enumerate_from_poscar but excludes superperiodic structures."""
    ...

def enumerate_from_poscar_at_index(
    poscar_path: str,
    index: int,
    composition: List[int],
    symprec: float = 1e-5,
) -> List[Coloring]:
    """Enumerate colorings across all inequivalent supercells at a volume index.

    Reads the parent structure's symmetry from the POSCAR file, then sweeps all
    inequivalent supercell shapes at the given volume index.
    """
    ...

def enumerate_from_poscar_at_index_primitive(
    poscar_path: str,
    index: int,
    composition: List[int],
    symprec: float = 1e-5,
) -> List[Coloring]:
    """Like enumerate_from_poscar_at_index but excludes superperiodic structures."""
    ...

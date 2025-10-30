# dna_montecarlo.py
import random
import math
import csv
from typing import Tuple, List, Iterable
import numpy as np

# ------------------------
# Task 0: simulation box
# ------------------------
def define_box(xmin: float, xmax: float, ymin: float, ymax: float, zmin: float, zmax: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Return (min_corner, max_corner) as numpy arrays in the chosen units (ANGSTROM recommended).
    """
    min_corner = np.array([xmin, ymin, zmin], dtype=float)
    max_corner = np.array([xmax, ymax, zmax], dtype=float)
    assert np.all(max_corner > min_corner), "Box max must be larger than min on each axis"
    return min_corner, max_corner

# ------------------------
# Task 1: random point in box
# ------------------------
def random_point_in_box(min_corner: np.ndarray, max_corner: np.ndarray) -> np.ndarray:
    """Uniform random point inside axis aligned box."""
    r = np.random.random(3)
    return min_corner + r * (max_corner - min_corner)

# ------------------------
# Task 2: random sphere
# ------------------------
def random_sphere_in_box(min_corner: np.ndarray, max_corner: np.ndarray, min_radius: float, max_radius: float) -> Tuple[np.ndarray, float]:
    """
    Return (center, radius). Ensures center is placed so sphere lies fully inside the box.
    """
    assert max_radius > 0 and min_radius >= 0 and max_radius >= min_radius
    r = random.uniform(min_radius, max_radius)
    # choose center so the entire sphere fits inside box
    low = min_corner + r
    high = max_corner - r
    assert np.all(high >= low), "Box too small for sphere of chosen radius range"
    center = np.array([random.uniform(low[i], high[i]) for i in range(3)], dtype=float)
    return center, r

# ------------------------
# Task 3: point inside sphere
# ------------------------
def point_in_sphere(point: np.ndarray, center: np.ndarray, radius: float) -> bool:
    return np.sum((point - center)**2) <= radius**2

# Vectorized version for many points
def points_in_sphere(points: np.ndarray, center: np.ndarray, radius: float) -> np.ndarray:
    """Return boolean mask for points inside sphere. points: (N,3) array."""
    d2 = np.sum((points - center)**2, axis=1)
    return d2 <= radius**2

# ------------------------
# Monte Carlo helpers
# ------------------------
def estimate_fraction_inside(objects: Iterable[Tuple[np.ndarray, float]],
                             min_corner: np.ndarray,
                             max_corner: np.ndarray,
                             n_points: int,
                             rng_seed: int = None) -> float:
    """
    objects: iterable of (center, radius)
    Returns fraction of random points that are inside at least one object.
    """
    if rng_seed is not None:
        np.random.seed(rng_seed)
    pts = np.random.random((n_points, 3)) * (max_corner - min_corner) + min_corner
    inside_any = np.zeros(n_points, dtype=bool)
    for center, radius in objects:
        inside_any |= points_in_sphere(pts, center, radius)
    frac = inside_any.sum() / float(n_points)
    return frac

def box_volume(min_corner: np.ndarray, max_corner: np.ndarray) -> float:
    return float(np.prod(max_corner - min_corner))

# ------------------------
# Task 4 & 5: single sphere check and pi estimate
# ------------------------
def single_sphere_test(center: np.ndarray, radius: float, min_corner: np.ndarray, max_corner: np.ndarray, n_points: int=100000, rng_seed: int=None):
    """
    Monte Carlo estimate of sphere volume and comparison with analytical value.
    Also can estimate pi via slice method: use random (x,y) inside bounding square to estimate circle area if desired.
    Returns dict with montecarlo_volume, analytical_volume, relative_error.
    """
    frac = estimate_fraction_inside([(center, radius)], min_corner, max_corner, n_points, rng_seed)
    v_box = box_volume(min_corner, max_corner)
    v_est = frac * v_box
    v_analytical = 4/3 * math.pi * radius**3
    rel_err = abs(v_est - v_analytical) / v_analytical
    return {"v_est": v_est, "v_analytical": v_analytical, "rel_err": rel_err, "fraction": frac}

def estimate_pi_via_montecarlo_in_square(n_points: int=200000, rng_seed: int=None) -> float:
    """
    Classic 2D Monte Carlo: inside unit circle in square [-1,1]x[-1,1].
    """
    if rng_seed is not None:
        np.random.seed(rng_seed)
    pts = np.random.uniform(-1, 1, size=(n_points, 2))
    inside = np.sum(np.sum(pts**2, axis=1) <= 1.0)
    pi_est = 4.0 * inside / n_points
    return pi_est

# ------------------------
# Task 6 & 7: many random spheres
# ------------------------
def generate_random_spheres(min_corner: np.ndarray, max_corner: np.ndarray, n_spheres: int, min_radius: float, max_radius: float) -> List[Tuple[np.ndarray, float]]:
    return [random_sphere_in_box(min_corner, max_corner, min_radius, max_radius) for _ in range(n_spheres)]

# ------------------------
# Task 8: read DNA coordinates file
# ------------------------
def read_xyz_with_types(filename: str) -> List[Tuple[np.ndarray, str]]:
    """
    Expected file format: whitespace csv with columns x y z element (angstrom). 
    Skips empty or comment lines. Returns list of (np.array([x,y,z]), element_symbol).
    """
    rows = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            # accept lines with at least 4 columns and last is element
            if len(parts) < 4:
                continue
            x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
            elem = parts[3]
            rows.append((np.array([x, y, z], dtype=float), elem))
    return rows

# ------------------------
# Helper: element -> covalent/van-der-Waals radius (angstrom)
# Use a small table; extend as needed.
# ------------------------
ELEMENT_RADII = {
    'H': 1.20,   # van der Waals approximate
    'C': 1.70,
    'N': 1.55,
    'O': 1.52,
    'P': 1.80,
    'S': 1.80,
    # add others if needed
}

def element_radius(element: str) -> float:
    el = element.capitalize()
    if el in ELEMENT_RADII:
        return ELEMENT_RADII[el]
    raise KeyError(f"Element radius not available for {element}")

# ------------------------
# Task 9: build objects list from DNA file
# ------------------------
def atoms_to_spheres(atom_list: Iterable[Tuple[np.ndarray, str]]) -> List[Tuple[np.ndarray, float]]:
    return [(pos, element_radius(elem)) for pos, elem in atom_list]

def bounding_box_for_atoms(atom_list: Iterable[Tuple[np.ndarray, str]], padding: float = 2.0) -> Tuple[np.ndarray, np.ndarray]:
    positions = np.array([pos for pos, _ in atom_list])
    mins = positions.min(axis=0) - padding
    maxs = positions.max(axis=0) + padding
    return mins, maxs

# ------------------------
# Task 10: DNA volume by Monte Carlo
# ------------------------
def dna_volume_montecarlo(atom_file: str, n_points: int = 200_000, rng_seed: int = None) -> dict:
    atoms = read_xyz_with_types(atom_file)
    spheres = atoms_to_spheres(atoms)
    min_corner, max_corner = bounding_box_for_atoms(atoms)
    frac = estimate_fraction_inside(spheres, min_corner, max_corner, n_points, rng_seed)
    vol = frac * box_volume(min_corner, max_corner)
    return {
        "n_atoms": len(spheres),
        "box_min": min_corner,
        "box_max": max_corner,
        "fraction_inside": frac,
        "estimated_volume": vol
    }

# ------------------------
# Simple asserts for basic tests
# ------------------------
def _self_test():
    # simple box and sphere
    minc, maxc = define_box(0, 10, 0, 10, 0, 10)
    c, r = np.array([5.0, 5.0, 5.0]), 2.0
    res = single_sphere_test(c, r, minc, maxc, n_points=50000, rng_seed=1)
    assert res["rel_err"] < 0.05, "Monte Carlo single sphere test worse than 5% (increase n_points)"
    # pi estimate sanity
    pi_est = estimate_pi_via_montecarlo_in_square(n_points=50000, rng_seed=1)
    assert abs(pi_est - math.pi) < 0.05, "Pi estimate too far (increase n_points)"
    print("Self-tests passed")

if __name__ == "__main__":
    _self_test()

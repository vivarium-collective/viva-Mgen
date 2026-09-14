"""Shared per-site chromosome state (bin-indexed) — the viva-native stand-in for
Karr 2012's CircularSparseMat. Phase 1 implements the LESION layer: a lesion map
{bin(str) -> count(float)} that DNADamage adds to and DNARepair clears from,
shared as one structure. Later phases add bound-protein footprints, per-region
linking number, and polymerized regions on the same bin index (see
docs/superpowers/specs/2026-09-14-unified-chromosome-design.md).
"""
from __future__ import annotations


def empty_lesion_map(n_bins: int) -> dict:
    """Pre-seeded all-zero lesion map (so additive map[float] deltas land)."""
    return {str(b): 0.0 for b in range(int(n_bins))}


def add_lesions(rng, n_bins: int, count) -> dict:  # rng: numpy.random.Generator
    """Additive delta placing ``count`` independent lesions at random bins."""
    count = int(count)
    if count <= 0:
        return {}
    delta: dict = {}
    for b in rng.integers(0, int(n_bins), size=count):
        delta[str(int(b))] = delta.get(str(int(b)), 0.0) + 1.0
    return delta


def n_lesions(lesion_map) -> float:
    return float(sum(float(v) for v in (lesion_map or {}).values()))


def repair_sites(lesion_map, capacity, rng):  # rng: numpy.random.Generator
    """Repair up to ``capacity`` lesion instances from bins with count>0, weighted
    by count. Returns (negative additive delta, number repaired). Never repairs
    more than present, never drives a bin below zero."""
    present = {b: float(v) for b, v in (lesion_map or {}).items() if float(v) > 0.0}
    total = sum(present.values())
    cap = int(min(float(capacity), total))
    if cap <= 0:
        return {}, 0.0
    # expand to per-instance bin list, sample `cap` without replacement
    bins = []
    for b, v in present.items():
        bins.extend([b] * int(v))
    chosen = rng.choice(len(bins), size=cap, replace=False)
    delta: dict = {}
    for i in chosen:
        b = bins[int(i)]
        delta[b] = delta.get(b, 0.0) - 1.0
    return delta, float(cap)


def lesion_positions(lesion_map, genome_length_bp: float, n_bins: int) -> list[float]:
    """bp coordinates (bin start coordinates) of currently-damaged bins."""
    bp_per_bin = float(genome_length_bp) / int(n_bins)
    return [int(b) * bp_per_bin for b, v in (lesion_map or {}).items() if float(v) > 0.0]

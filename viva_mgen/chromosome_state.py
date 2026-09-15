"""Shared per-site chromosome state (bin-indexed) — the viva-native stand-in for
Karr 2012's CircularSparseMat. Phase 1 implements the LESION layer: a lesion map
{bin(str) -> count(float)} that DNADamage adds to and DNARepair clears from,
shared as one structure. Later phases add bound-protein footprints, per-region
linking number, and polymerized regions on the same bin index (see
docs/superpowers/specs/2026-09-14-unified-chromosome-design.md).
"""
from __future__ import annotations

# Default bin resolution of the shared chromosome structure, matching
# ChromosomeDynamics' n_bins.
N_CHROMOSOME_BINS = 580


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


# ---------------------------------------------------------------------------
# LINKING-NUMBER layer (phase 2): per-region superhelical density on the same
# chromosome, the viva-native stand-in for CircularSparseMat's per-region
# linking numbers. Karr partitions the chromosome into topological regions
# (bounded by the replication forks and bound proteins) each with its own linking
# number that gyrase/topoisomerases relax independently; the genome-wide σ the
# DNASupercoiling observable reports is the mean over regions. Until replication/
# transcription wire per-region perturbation in (staged later phases), the regions
# are homogeneous, so mean == the former single-pool σ (no behavior change).
# ---------------------------------------------------------------------------
N_SUPERCOIL_REGIONS = 20  # topological regions the chromosome is partitioned into


# ---------------------------------------------------------------------------
# POLYMERIZED-REGIONS layer (phase 4): which bins the replication fork has copied,
# the CircularSparseMat "polymerized" mask. Two forks diverge bidirectionally from
# oriC to terC; a bin is polymerized once a fork has passed it. Σ(mask)/n_bins is
# the replicated_fraction the Replication submodel already reports, so the mask is
# a per-site VIEW of the same scalar (no behaviour change to replicated_fraction).
# ---------------------------------------------------------------------------
def empty_polymerized_map(n_bins: int) -> dict:
    return {str(b): 0.0 for b in range(int(n_bins))}


def fork_polymerized(fraction: float, n_bins: int, oriC_bin: int = 0) -> dict:
    """Per-bin polymerized mask {bin -> 1.0} for a bidirectional oriC→terC fork
    that has copied ``fraction`` of the chromosome: the ``round(fraction·n_bins)``
    bins closest to oriC (by circular distance) are marked replicated. Σ of the
    mask is ``round(fraction·n_bins)`` so mask-fraction ≈ ``fraction``."""
    nb = int(n_bins)
    n_done = max(0, min(nb, int(round(float(fraction) * nb))))
    if n_done <= 0:
        return {}
    order = sorted(range(nb), key=lambda b: min((b - oriC_bin) % nb, (oriC_bin - b) % nb))
    return {str(b): 1.0 for b in order[:n_done]}


def empty_linking_map(n_regions: int, sigma0: float = 0.0) -> dict:
    """Per-region superhelical density {region_index(str) -> σ(float)}, all
    initialized to ``sigma0`` (pre-seeded so additive/overwrite deltas land)."""
    return {str(r): float(sigma0) for r in range(int(n_regions))}


def mean_sigma(linking_map) -> float:
    """Genome-wide σ = mean of the per-region linking numbers (the scalar the
    ``superhelical_density`` observable reports)."""
    vals = [float(v) for v in (linking_map or {}).values()]
    return sum(vals) / len(vals) if vals else 0.0


def relax_regions(linking_map, setpoint, total_acts, turns_per_region, rng=None):
    """Distribute ``total_acts`` gyrase supercoiling acts across the regions and
    return the per-region σ delta (overwrite semantics: new σ per region). Each
    act moves one region's σ by 1/turns_per_region toward ``setpoint``; acts are
    spread over regions that are still short of the setpoint, so no region is
    driven past it. Returns ``{region: new_sigma}`` for regions that changed."""
    regions = {r: float(v) for r, v in (linking_map or {}).items()}
    if not regions or total_acts <= 0:
        return {}
    short = [r for r, s in regions.items() if abs(setpoint - s) > 1e-12]
    if not short:
        return {}
    per_turn = 1.0 / float(turns_per_region) if turns_per_region else 0.0
    acts_left = float(total_acts)
    order = list(short)
    if rng is not None:
        rng.shuffle(order)
    out: dict = {}
    # round-robin one act at a time keeps regions balanced; cap per region so σ
    # never overshoots the setpoint.
    while acts_left > 0 and order:
        for r in list(order):
            if acts_left <= 0:
                break
            s = regions[r]
            gap = setpoint - s
            if abs(gap) <= 1e-12:
                order.remove(r)
                continue
            step = min(per_turn, abs(gap)) * (1.0 if gap > 0 else -1.0)
            regions[r] = s + step
            out[r] = regions[r]
            acts_left -= 1.0
            if abs(setpoint - regions[r]) <= 1e-12:
                order.remove(r)
    return out

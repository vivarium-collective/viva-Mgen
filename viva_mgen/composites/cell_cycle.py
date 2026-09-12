"""Composite generator for the cell-cycle-regulation study (Fig 4)."""

from __future__ import annotations

from process_bigraph.composite_generator import composite_generator


@composite_generator(
    name="fig4_cell_cycle",
    description="Fig 4 — emergent cell-cycle regulation: dNTP-buffered three-phase replication (initiation → replication → done) where initiation-phase dNTP surplus controls replication speed.",
    parameters={
        "initial_dnaA": {"type": "float", "default": 0.0,
                         "description": "DnaA monomers present at cell birth (shortens initiation)"},
        "initial_dntp": {"type": "float", "default": 0.0,
                         "description": "dNTP molecules present at cell birth"},
        "seed": {"type": "integer", "default": 0, "description": "Random seed"},
    },
    emitters=[{"address": "local:ParquetEmitter",
               "paths": ["stores/replicated_fraction", "stores/dntp_pool",
                         "stores/phase_code", "stores/dnaA_complex"]}],
)
def fig4_cell_cycle(core=None, *, initial_dnaA=0.0, initial_dntp=0.0, seed=0):
    stores = {
        "dntp_synthesis_scale": 1.0,
        "replicated_fraction": 0.0,
        "dntp_pool": 0.0,
        "dnaA_complex": 0.0,
        "phase_code": 0.0,
        "chromosome_copy": 1.0,
        "initiation_duration": 0.0,
        "replication_duration": 0.0,
        "dntp_at_replication_start": 0.0,
    }
    return {
        "stores": stores,
        "replication": {
            "_type": "process",
            "address": "local:ReplicationReproductionProcess",
            "config": {"initial_dnaA": initial_dnaA, "initial_dntp": initial_dntp, "seed": seed},
            "interval": 1.0,
            "inputs": {"dntp_synthesis_scale": ["stores", "dntp_synthesis_scale"]},
            "outputs": {
                "replicated_fraction": ["stores", "replicated_fraction"],
                "dntp_pool": ["stores", "dntp_pool"],
                "dnaA_complex": ["stores", "dnaA_complex"],
                "phase_code": ["stores", "phase_code"],
                "chromosome_copy": ["stores", "chromosome_copy"],
                "initiation_duration": ["stores", "initiation_duration"],
                "replication_duration": ["stores", "replication_duration"],
                "dntp_at_replication_start": ["stores", "dntp_at_replication_start"],
            },
        },
        "emitter": {
            "_type": "step",
            "address": "local:RAMEmitter",
            "config": {"emit": {
                "replicated_fraction": "float", "dntp_pool": "float",
                "dnaA_complex": "float", "phase_code": "float",
                "chromosome_copy": "float", "initiation_duration": "float",
                "replication_duration": "float", "dntp_at_replication_start": "float",
            }},
            "inputs": {k: ["stores", k] for k in (
                "replicated_fraction", "dntp_pool", "dnaA_complex", "phase_code",
                "chromosome_copy", "initiation_duration", "replication_duration",
                "dntp_at_replication_start")},
        },
    }

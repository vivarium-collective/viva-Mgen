"""Parameter calculator (ParCa) as a process — reproduction of Karr 2012
``FitConstants``, exposed so it can be a study's baseline composite that PRECEDES
the figure studies (which consume its fitted parameters).

Unlike the other reproduction processes (which simulate a submodel over time), this
one runs the parameter *fit* and reports its summary + validation: it computes the
per-gene expression/decay/synthesis panel from the real observed data
(:mod:`viva_mgen.parca`), the metabolic demand the expression implies, and the
expression↔metabolism feasibility closure. Its outputs are the report-card scalars
the ``parca-parameter-fitting`` study evaluates.
"""
from __future__ import annotations

import statistics

from process_bigraph import Process

from .. import parca
from ..kb import load_gene_expression, load_genes


class ParameterCalculatorReproductionProcess(Process):
    """Parameter calculator — reproduction of Karr 2012 FitConstants: runs the fit and reports its summary + metabolic-feasibility closure.

    Contract — in: none (reads the knowledge-base-derived datasets). out (snapshots):
    n_genes, median_mrna_synthesis_rate, mrna_halflife_min_mean, rrna_halflife_min,
    nmp_au_fraction, closed_loop_feasible, nmp_supply_flux, aa_supply_flux,
    growth_baseline.
    Fidelity: reproduces the FitConstants initialization + analytic QP (RNA-mass,
    DnaA/FtsZ held, net-supercoiling) + the metabolic-demand feasibility closure,
    from natively-decoded real knowledge-base data.
    """

    description = (
        "Parameter calculator (ParCa) — reproduction of Karr 2012 FitConstants.\n"
        "Computes the per-gene expression/decay/synthesis panel from the real observed "
        "knowledge-base data, fits it under the FitConstants linear constraints "
        "(RNA-mass balance, DnaA/FtsZ expression held, net-supercoiling zero), and "
        "closes the expression↔metabolism loop against the iPS189 FBA network.\n"
        "Contract — in: none (reads the decoded knowledge-base datasets). out (snapshots): "
        "n_genes, median_mrna_synthesis_rate, mrna_halflife_min_mean, rrna_halflife_min, "
        "nmp_au_fraction, closed_loop_feasible, nmp_supply_flux, aa_supply_flux, growth_baseline.\n"
        "Fidelity: genuine inputs (natively-decoded real KB data) + the FitConstants "
        "initialization and analytic QP + the metabolic-feasibility closure; the nonlinear "
        "metabolism resource coupling and heuristic refinement are the reduced part."
    )

    config_schema = {"seed": {"_type": "integer", "_default": 0}}

    def inputs(self):
        return {}

    def outputs(self):
        return {
            "n_genes": "overwrite[float]",
            "median_mrna_synthesis_rate": "overwrite[float]",
            "mrna_halflife_min_mean": "overwrite[float]",
            "rrna_halflife_min": "overwrite[float]",
            "nmp_au_fraction": "overwrite[float]",
            "closed_loop_feasible": "overwrite[float]",
            "nmp_supply_flux": "overwrite[float]",
            "aa_supply_flux": "overwrite[float]",
            "growth_baseline": "overwrite[float]",
        }

    def initial_state(self):
        return {}

    def update(self, state, interval):
        panel = parca.calculate_parameters()
        expr = load_gene_expression()
        types = {}
        for g in load_genes():
            key = (g.get("symbol") or "").strip() or g["gene_id"]
            types[key] = (expr.get(g["gene_id"], {}).get("rna_type") or "").strip()
        synth = sorted(v[0] for v in panel.values())
        mrna_hl = [v[1] / 60.0 for k, v in panel.items() if types.get(k) == "mRNA"]
        rrna_hl = [v[1] / 60.0 for k, v in panel.items() if types.get(k) == "rRNA"]
        demand = parca.metabolic_demand()
        au = (demand.get("nmp", {}).get("A", 0.0) + demand.get("nmp", {}).get("U", 0.0)) if demand else 0.0
        loop = parca.close_metabolic_loop()
        return {
            "n_genes": float(len(panel)),
            "median_mrna_synthesis_rate": float(synth[len(synth) // 2]) if synth else 0.0,
            "mrna_halflife_min_mean": float(statistics.mean(mrna_hl)) if mrna_hl else 0.0,
            "rrna_halflife_min": float(min(rrna_hl)) if rrna_hl else 0.0,
            "nmp_au_fraction": float(au),
            "closed_loop_feasible": 1.0 if loop.get("feasible") else 0.0,
            "nmp_supply_flux": float(loop.get("nmp_supply_flux", 0.0)),
            "aa_supply_flux": float(loop.get("aa_supply_flux", 0.0)),
            "growth_baseline": float(loop.get("growth_baseline", 0.0)),
        }

"""Metabolism submodel — clean-room reproduction.

Reproduces the Karr et al. 2012 whole-cell model's ``Metabolism`` submodel:
a flux-balance-analysis (FBA) process that, each 1 s timestep, maximizes the
biomass objective over the metabolic network subject to enzyme/directionality
bounds, sets the cell's instantaneous growth rate, and reports the reaction
fluxes and energy (ATP/GTP) production used by the rest of the cell.

This is **not** the original MATLAB code. The metabolic network is the
published Suthers et al. 2009 *M. genitalium* reconstruction **iPS189**
(BioModels MODEL1507180052) — the same reconstruction the original metabolism
submodel was built on — solved with :mod:`cobra`.

Gene disruptions (Fig 6) are modelled by knocking out the disrupted genes'
reactions; kinetic-parameter changes (Fig 7) by scaling a reaction's flux
bound. Both use cobra's reverting ``with model:`` context so the cached model
is never mutated.
"""

from __future__ import annotations

import functools

from process_bigraph import Process

from ..kb import load_metabolic_model, normalize_gene_id


@functools.lru_cache(maxsize=8)
def _wt_growth(sbml_path=None) -> float:
    model = load_metabolic_model(sbml_path)
    sol = model.optimize()
    return float(sol.objective_value or 0.0)


class MetabolismFbaReproductionProcess(Process):
    """FBA metabolism over iPS189 (reproduction of Karr 2012 Metabolism).

    Inputs
    ------
    nutrient_scale : float
        Multiplier (default 1.0) on the model's carbon-exchange uptake upper
        bounds — lets a sibling process throttle nutrient availability.

    Outputs
    -------
    growth_rate : overwrite[float]
        Current biomass-objective flux (raw iPS189 units). A current-value
        sensor reading, hence ``overwrite`` (only this process writes it).
    growth_fraction : overwrite[float]
        growth_rate / wild-type growth_rate — the calibrated, unit-free growth
        the mass submodel integrates (1.0 == unperturbed).
    atp_production : overwrite[float]
        ATP synthase (ATPS4r) flux — ATP made per unit time.
    gtp_production : overwrite[float]
        GTP-producing flux (nucleoside-diphosphate-kinase family), current.
    feasible : overwrite[float]
        1.0 if the LP was optimal with positive growth, else 0.0.
    """

    description = (
        "Flux-balance-analysis metabolism — reproduction of the Karr 2012 Metabolism "
        "submodel.\n"
        "Each 1 s step it maximizes the iPS189 biomass objective\n"
        "    max cᵀv   s.t.   S·v = 0,   lb ≤ v ≤ ub\n"
        "over the published Suthers 2009 M. genitalium reconstruction (351 reactions, "
        "346 metabolites, 126 genes; solved with COBRA). Gene disruptions (Fig 6) knock "
        "out the disrupted genes' reactions; a kinetic-parameter change (Fig 7) scales a "
        "reaction's flux bound as a kcat/Vmax proxy. Both use a reverting model context so "
        "the cached network is never mutated.\n"
        "Contract — in: nutrient_scale (carbon-uptake multiplier a sibling can throttle). "
        "out: growth_rate (biomass flux, raw iPS189 units), growth_fraction "
        "(growth_rate ÷ wild-type, the calibrated unit-free growth the mass submodel "
        "integrates), atp_production (ATP-synthase ATPS4r flux), gtp_production "
        "(summed GTP-linked kinase/transport flux: NDPK1/NDPK2 nucleoside-diphosphate "
        "kinases, GK1 guanylate kinase, GTPtp transport), feasible ∈ {0,1}.\n"
        "Fidelity: FULL — the genuine published reconstruction the WCM's metabolism "
        "submodel was built on."
    )

    config_schema = {
        "sbml_path": {"_type": "string", "_default": ""},
        # gene ids to disrupt (Fig 6); MG_001 / MG001 both accepted
        "disrupted_genes": {"_type": "list[string]", "_default": []},
        # reaction_id -> multiplicative scale on its flux upper bound (Fig 7)
        "reaction_bound_scale": {"_type": "map[float]", "_default": {}},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._model = load_metabolic_model(self.config["sbml_path"] or None)
        self._wt = _wt_growth(self.config["sbml_path"] or None) or 1.0
        # normalized-id -> model gene id
        self._gene_index = {normalize_gene_id(g.id): g.id for g in self._model.genes}

    def inputs(self):
        return {"nutrient_scale": "float"}

    def outputs(self):
        return {
            "growth_rate": "overwrite[float]",
            "growth_fraction": "overwrite[float]",
            "atp_production": "overwrite[float]",
            "gtp_production": "overwrite[float]",
            "feasible": "overwrite[float]",
        }

    def initial_state(self):
        return {"nutrient_scale": 1.0}

    def _flux(self, sol, rxn_id, default=0.0):
        try:
            return float(sol.fluxes.get(rxn_id, default))
        except Exception:
            return default

    def update(self, state, interval):
        model = self._model
        scale = float(state.get("nutrient_scale", 1.0) or 1.0)
        with model:  # reverting context — cached model stays pristine
            # gene disruptions
            for gid in self.config.get("disrupted_genes", []) or []:
                mid = self._gene_index.get(normalize_gene_id(gid))
                if mid is not None:
                    model.genes.get_by_id(mid).knock_out()
            # kinetic-parameter (flux-bound) scaling
            for rid, s in (self.config.get("reaction_bound_scale", {}) or {}).items():
                if rid in model.reactions:
                    r = model.reactions.get_by_id(rid)
                    if r.upper_bound > 0:
                        r.upper_bound = r.upper_bound * float(s)
                    if r.lower_bound < 0:
                        r.lower_bound = r.lower_bound * float(s)
            # nutrient throttle on carbon uptake (negative lower bounds)
            if scale != 1.0:
                for r in model.reactions:
                    if r.boundary and r.lower_bound < 0:
                        r.lower_bound = r.lower_bound * scale
            sol = model.optimize()

        growth = float(sol.objective_value or 0.0) if sol.status == "optimal" else 0.0
        atp = abs(self._flux(sol, "ATPS4r")) if sol.status == "optimal" else 0.0
        # GTP production: nucleoside diphosphate kinase (GTP) family
        gtp = 0.0
        if sol.status == "optimal":
            for rid in ("NDPK1", "NDPK2", "GK1", "GTPtp"):
                gtp += abs(self._flux(sol, rid))
        return {
            "growth_rate": growth,
            "growth_fraction": growth / self._wt if self._wt else 0.0,
            "atp_production": atp,
            "gtp_production": gtp,
            "feasible": 1.0 if growth > 1e-6 else 0.0,
        }

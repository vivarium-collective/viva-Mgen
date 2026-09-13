"""Cytokinesis & host-interaction submodels — clean-room reproductions (reduced).

Reproduces three Karr 2012 *M. genitalium* whole-cell model submodels that
together carry the cell through the final phase of the cell cycle and report on
its interaction with the host:

* ``FtsZPolymerization`` — in the cytosol FtsZ cycles through states (inactive
  monomer, GDP-bound / deactivated, GTP-bound / activated, nucleated dimer, and
  elongated filament of ≥3 activated subunits). GTP-driven activation, nucleation
  and elongation assemble FtsZ-GTP into the contractile Z-ring at the division
  site (Surovtsev 2008 kinetic scheme).
* ``Cytokinesis`` — the Z-ring pinches the membrane by a cycle of filament
  binding, GTP-hydrolysis-driven bending, and dissociation (Li 2007). Each cycle
  shrinks the pinched (septum) diameter; the cell divides when the diameter
  reaches zero. In the original whole-cell model **reaching septum diameter 0 is
  the event that ends the simulation** — it is the division trigger of the whole
  cell.
* ``HostInteraction`` — qualitative, reporting-only account of *M. genitalium*
  adhesion to the host urogenital epithelium, mediated by the terminal organelle
  and its adhesin lipoproteins (MG_191/MgPa, MG_192, MG_217, MG_318, …).

These are reduced but mechanistically genuine reproductions — the FtsZ-GTP →
filament → ring-contraction → division chain is real — not the full
ODE/discretization + polygon-pinching + Boolean-host machinery of the originals.
"""

from __future__ import annotations

import numpy as np
from process_bigraph import Process, Step

from .. import constants as C


class FtsZPolymerizationReproductionProcess(Process):
    """GTP-driven assembly of FtsZ monomers into the contractile Z-ring.

    Inputs (current-value sensors)
    ------------------------------
    ftsz_monomer : float
        Free (inactive) FtsZ monomer copy number available to be activated.
    gtp : float
        Free GTP molecule count — limits activation and elongation.

    Outputs
    -------
    ftsz_ring_filaments : overwrite[float]
        Current number of FtsZ subunits assembled into Z-ring filaments
        (the Z-ring pool that Cytokinesis contracts), snapshot.
    gtp : float
        Negative delta — GTP consumed binding newly activated FtsZ this step.
    """

    description = (
        "FtsZ polymerization — reduced reproduction of Karr 2012 FtsZPolymerization "
        "(Surovtsev 2008 kinetic scheme).\n"
        "FtsZ cycles through states: inactive monomer → GTP-activated monomer (binds a GTP) → "
        "nucleated/elongated filament. Reduced mechanism: activation converts free FtsZ to "
        "FtsZ-GTP at a rate limited by both FtsZ monomer and free GTP (1 GTP bound per "
        "activation); FtsZ-GTP then elongates into the Z-ring filament pool; a GDP↔GTP "
        "exchange / dissociation flux returns a fraction of the ring to the monomer cycle.\n"
        "Contract — in (sensors): ftsz_monomer, gtp. out: ftsz_ring_filaments (Z-ring subunit "
        "count snapshot), gtp (negative delta, consumed).\n"
        "Fidelity: REDUCED mechanism (genuine GTP-driven state-cycle → ring assembly; not the "
        "full multi-mer ODE + stochastic discretization of the original)."
    )

    config_schema = {
        # k_on for the conserved free⇌ring relaxation = real FtsZ-GTP forward
        # activation rate (Karr 2012 parameters.json FtsZPolymerization.activationFwd = 1.1 /s).
        "activation_rate": {"_type": "float", "_default": 1.1},
        # fraction of the FtsZ-GTP pool elongating into ring filaments per second.
        # (The real mass-action nucleation/elongation/exchange fwd+rev constants —
        # nucleationFwd 4.2e6 / nucleationRev 40, elongationFwd 5.1e6 / elongationRev 2.9,
        # exchangeFwd 1e4 / exchangeRev 5e3 — are in different (bimolecular) units with no
        # slot in this reduced conserved-pool relaxation; they are available via
        # viva_mgen.kb.karr_process_params("FtsZPolymerization") but not forced into the math.)
        "elongation_rate": {"_type": "float", "_default": 0.5},
        # k_off for the free⇌ring relaxation = real FtsZ-GTP reverse activation rate
        # (Karr 2012 parameters.json FtsZPolymerization.activationRev = 0.01 /s).
        "dissociation_rate": {"_type": "float", "_default": 0.01},
        "initial_ftsz_gtp": {"_type": "float", "_default": 0.0},
        "initial_ring": {"_type": "float", "_default": 0.0},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._ring = float(self.config["initial_ring"])          # subunits in ring filaments
        self._free = None                                        # free FtsZ (set from supply)

    def inputs(self):
        return {"ftsz_monomer": "float", "gtp": "float"}

    def outputs(self):
        return {
            "ftsz_ring_filaments": "overwrite[float]",
            "gtp": "float",
        }

    def initial_state(self):
        return {"ftsz_monomer": 0.0, "gtp": 0.0}

    def update(self, state, interval):
        supply = max(float(state.get("ftsz_monomer", 0.0)), 0.0)
        gtp = max(float(state.get("gtp", 0.0)), 0.0)

        # FtsZ is a CONSERVED pool partitioned between free subunits and the ring;
        # ``ftsz_monomer`` is the total FtsZ available. Modelling the exchange on a
        # conserved pool (rather than an unbounded first-order source) makes the
        # steady-state ring size ≈ supply·k_on/(k_on+k_off) — independent of the
        # timestep — instead of blowing up or collapsing with Δt.
        if self._free is None:
            self._free = supply
        total = self._free + self._ring
        if supply > total:                 # newly synthesised FtsZ enters the free pool
            self._free += supply - total
            total = supply

        # Exact two-state relaxation of the conserved pool (free ⇌ ring) over the
        # interval: ring → ring_eq = total·k_on/(k_on+k_off) with time constant
        # 1/(k_on+k_off). Solving the ODE closed-form makes the result identical at
        # any Δt (no operator-splitting error), so the ring reaches its functional
        # size (≈ full pool, k_on ≫ k_off) and drives the septum to close.
        k_on = max(float(self.config["activation_rate"]), 0.0)
        k_off = max(float(self.config["dissociation_rate"]), 0.0)
        ksum = k_on + k_off
        ring_eq = total * (k_on / ksum) if ksum > 0 else self._ring
        ring_new = ring_eq + (self._ring - ring_eq) * np.exp(-ksum * interval)

        # GTP cost: one GTP per net new ring subunit; cap growth if GTP-limited.
        added = max(0.0, ring_new - self._ring)
        gtp_consumed = min(added, gtp)
        if added > gtp:
            ring_new = self._ring + gtp
        self._ring = ring_new
        self._free = total - self._ring

        return {
            "ftsz_ring_filaments": max(self._ring, 0.0),
            "gtp": -gtp_consumed,  # negative delta: GTP consumed this step
        }


class CytokinesisReproductionProcess(Process):
    """Z-ring contraction pinches the septum to zero → cell divides.

    Inputs (current-value sensors)
    ------------------------------
    ftsz_ring_filaments : float
        FtsZ subunits assembled into the Z-ring (from FtsZPolymerization); the
        contractile machinery. No ring ⇒ no contraction.
    replicated_fraction : float
        Chromosome-replication progress in [0, 1]; contraction only proceeds once
        replication (and hence segregation) is complete, matching the original's
        ``chromosome.segregated`` guard.

    Outputs
    -------
    septum_diameter : overwrite[float]
        Current pinched (septum) diameter in nm, shrinking from the cell width
        toward 0, snapshot.
    contraction_progress : overwrite[float]
        Fraction of the initial width already pinched, in [0, 1], snapshot.
    divided : overwrite[float]
        1.0 once the septum diameter reaches 0 (division), else 0.0.
    """

    description = (
        "Cytokinesis — reduced reproduction of Karr 2012 Cytokinesis (Li 2007 bind/bend/"
        "dissociate model).\n"
        "The FtsZ-GTP Z-ring contracts by a cycle of filament binding, GTP-hydrolysis-driven "
        "bending, and dissociation; each cycle shrinks the pinched (septum) diameter. Reduced "
        "mechanism: once replication/segregation is complete (replicated_fraction ≥ threshold) "
        "and Z-ring filaments are present, the septum diameter shrinks from the M. genitalium "
        "cell width (~200 nm) toward 0 at a contraction rate scaled by ring availability; the "
        "cell divides when the diameter reaches 0.\n"
        "The ORIGINAL whole-cell model terminates the simulation when the septum diameter "
        "reaches 0 — this is THE division trigger of the whole-cell model.\n"
        "Contract — in (sensors): ftsz_ring_filaments, replicated_fraction. out: septum_diameter "
        "(nm snapshot), contraction_progress (0..1 snapshot), divided ∈ {0,1}.\n"
        "Fidelity: REDUCED mechanism (genuine ring-driven pinch-to-division; not the polygon-"
        "edge binding/bending/dissociation cycle of the original)."
    )

    config_schema = {
        # M. genitalium cell width — the initial septum diameter to pinch shut (nm)
        "initial_cell_width_nm": {"_type": "float", "_default": 200.0},
        # contraction rate (nm/s); default pinches ~200 nm over the fitted
        # cytokinesis duration (~3869 s) at full ring availability.
        # (The original's per-filament bind/bend/dissociate rate constants —
        # Karr 2012 parameters.json Cytokinesis.rateFilamentBindingMembrane 0.7,
        # rateFilamentDissociation 0.7, rateFtsZGtpHydrolysis 0.15 (all /s) — describe
        # the polygon-edge cycle the reduced pinch-rate formulation lumps into this
        # single nm/s speed; they have no direct slot here and are available via
        # viva_mgen.kb.karr_process_params("Cytokinesis").)
        "contraction_rate_nm_per_s": {"_type": "float", "_default": 200.0 / C.CYTOKINESIS_DURATION_S},
        # replication must be essentially complete before the ring can pinch
        "replication_complete_threshold": {"_type": "float", "_default": 0.999},
        # ring subunits needed for a fully effective (rate = 1x) contraction
        "ring_subunits_for_full_rate": {"_type": "float", "_default": 320.0},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._diameter = float(self.config["initial_cell_width_nm"])
        self._divided = 0.0

    def inputs(self):
        return {"ftsz_ring_filaments": "float", "replicated_fraction": "float"}

    def outputs(self):
        return {
            "septum_diameter": "overwrite[float]",
            "contraction_progress": "overwrite[float]",
            "divided": "overwrite[float]",
        }

    def initial_state(self):
        return {"ftsz_ring_filaments": 0.0, "replicated_fraction": 0.0}

    def update(self, state, interval):
        ring = max(float(state.get("ftsz_ring_filaments", 0.0)), 0.0)
        repl = float(state.get("replicated_fraction", 0.0))
        width = float(self.config["initial_cell_width_nm"])

        # Guard: contraction only proceeds when replication/segregation is complete
        # (original: `if ~this.chromosome.segregated; return; end`) and a ring exists.
        if repl >= self.config["replication_complete_threshold"] and ring > 0.0 and self._diameter > 0.0:
            # ring availability scales the pinch rate (bind/bend cycle throughput);
            # full rate once enough subunits are present to line the polygon edges.
            avail = min(1.0, ring / self.config["ring_subunits_for_full_rate"])
            self._diameter -= self.config["contraction_rate_nm_per_s"] * avail * interval

        # Division fires when the septum diameter reaches 0 (original: pinched
        # diameter < one filament length -> set to 0 -> division ends the sim).
        if self._diameter <= 0.0:
            self._diameter = 0.0
            self._divided = 1.0

        progress = 0.0 if width <= 0.0 else min(1.0, max(0.0, (width - self._diameter) / width))
        return {
            "septum_diameter": self._diameter,
            "contraction_progress": progress,
            "divided": self._divided,
        }


class HostInteractionReproductionProcess(Process):
    """Terminal-organelle-mediated adhesion of M. genitalium to host epithelium.

    Design choice — Process vs Step: the original ``HostInteraction`` is purely
    reporting (Boolean host-state rules, no substrate/mass flux), which would
    argue for a ``Step``. It is kept a ``Process`` here because it has a genuine
    time-varying input — ``terminal_organelle_fraction`` rises as the terminal
    organelle assembles over the cell cycle — so adhesion evolves in time rather
    than being a one-shot derived quantity. (``Step`` is imported and available
    should a caller prefer the report-only framing.)

    Inputs (current-value sensor)
    -----------------------------
    terminal_organelle_fraction : float
        Assembly / availability of the terminal-organelle core + adhesin
        lipoproteins, in [0, 1]; 1.0 == fully assembled with all adhesins present.

    Outputs
    -------
    adhesion_strength : overwrite[float]
        Qualitative adhesion strength to the host urogenital epithelium, in
        [0, 1], snapshot.
    """

    description = (
        "Host interaction — reduced reproduction of Karr 2012 HostInteraction (qualitative, "
        "reporting).\n"
        "The original decides adherence by a Boolean rule: adherent iff the terminal organelle "
        "is fully assembled AND all adhesin lipoproteins (MG_191/MgPa, MG_192, MG_217, MG_318) "
        "are present. Reduced mechanism: that all-or-nothing AND is softened to a cooperative "
        "(Hill-like) saturating function of terminal_organelle_fraction, so adhesion rises "
        "steeply toward 1 only as assembly approaches completion.\n"
        "Kept a Process (not a Step) because terminal_organelle_fraction is a genuine "
        "time-varying input, so adhesion_strength evolves over the cell cycle.\n"
        "Contract — in (sensor): terminal_organelle_fraction (0..1). out: adhesion_strength "
        "(0..1 snapshot).\n"
        "Fidelity: REDUCED / qualitative (as the original is qualitative; continuous softening "
        "of the Boolean adherence rule)."
    )

    config_schema = {
        # Hill coefficient — cooperativity of the softened AND rule (higher = more switch-like)
        "cooperativity": {"_type": "float", "_default": 4.0},
        # fraction at which adhesion is half-maximal
        "half_max_fraction": {"_type": "float", "_default": 0.7},
    }

    def inputs(self):
        return {"terminal_organelle_fraction": "float"}

    def outputs(self):
        return {"adhesion_strength": "overwrite[float]"}

    def initial_state(self):
        return {"terminal_organelle_fraction": 0.0}

    def update(self, state, interval):
        frac = min(max(float(state.get("terminal_organelle_fraction", 0.0)), 0.0), 1.0)
        # Hill function softens the Boolean "all adhesins present" AND rule: adhesion
        # approaches 1 only as the terminal organelle + adhesins near full assembly.
        n = self.config["cooperativity"]
        k = self.config["half_max_fraction"]
        strength = (frac ** n) / (frac ** n + k ** n) if frac > 0.0 else 0.0
        return {"adhesion_strength": min(max(strength, 0.0), 1.0)}

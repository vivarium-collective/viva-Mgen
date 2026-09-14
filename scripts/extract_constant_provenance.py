#!/usr/bin/env python
"""Write datasets/constant_provenance.json from viva_mgen.provenance.audit_constants()."""
from __future__ import annotations
import json
from collections import Counter
from pathlib import Path
from viva_mgen.provenance import audit_constants, SOURCE_TIERS

OUT = Path(__file__).resolve().parents[1] / "datasets" / "constant_provenance.json"

def main():
    a = audit_constants()
    tier = Counter()
    n = 0
    for consts in a.values():
        for rec in consts.values():
            tier[rec["source_tier"]] += 1
            n += 1
    payload = {"source": "viva_mgen.provenance.audit_constants() — kinetic-constant provenance vs Karr 2012 KB + supplement",
               "n_constants": n, "by_tier": {t: tier.get(t, 0) for t in SOURCE_TIERS},
               "constants": a}
    OUT.write_text(json.dumps(payload, indent=2))
    print(f"wrote {OUT}: {n} constants, by_tier={dict(tier)}")

if __name__ == "__main__":
    main()

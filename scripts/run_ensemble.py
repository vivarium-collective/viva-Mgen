#!/usr/bin/env python
"""Run an ensemble of the whole-cell composite and write an aggregated observable
(mean +/- spread + end-of-cycle distribution) — a demonstration of viva_mgen's
single-cell variation infrastructure.

  python scripts/run_ensemble.py --n 8 --hours 1.0 --observable mass \
      --out reports/ensemble_mass.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from viva_mgen.ensemble import run_ensemble, aggregate, final_values


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--hours", type=float, default=1.0)
    ap.add_argument("--observable", default="mass")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    ens = run_ensemble(args.n, args.hours * 3600.0)
    agg = aggregate(ens["cells"], args.observable)
    fv = [v for v in final_values(ens["cells"], args.observable) if v is not None]
    payload = {**agg, "seeds": ens["seeds"], "final_values": fv,
               "hours": args.hours}
    out = Path(args.out) if args.out else Path(f"reports/ensemble_{args.observable}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2))
    if fv:
        print(f"{args.observable}: {len(fv)} cells, final mean={np.mean(fv):.4g} "
              f"std={np.std(fv):.4g} (CV={np.std(fv)/np.mean(fv):.2%})" if np.mean(fv) else
              f"{args.observable}: {len(fv)} cells, final mean=0")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

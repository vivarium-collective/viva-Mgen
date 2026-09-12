# Contributing to viva-Mgen

## Development setup

uv is required. Install with `brew install uv` or `pip install uv`.

    uv venv .venv
    source .venv/bin/activate
    uv pip install -e ".[dev]"
    pytest

## What this repo is

viva-Mgen is a **clean-room reproduction** (`--reproduce`) of the Karr et al.
2012 *Mycoplasma genitalium* whole-cell model, re-expressed as native
[process-bigraph](https://github.com/vivarium-collective/process-bigraph)
Processes. It is **not** the original MATLAB code and does not bridge to it.
Every process class is named `*ReproductionProcess` and its module docstring
cites the original submodel it reproduces. Known divergences from the original
are listed in the README.

## Releasing to PyPI

Tag a commit with `git tag v<VERSION>` and push the tag. The
`.github/workflows/release.yml` workflow publishes to PyPI via trusted
publishing (no tokens after initial setup). See
https://docs.pypi.org/trusted-publishers/.

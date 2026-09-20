"""Assemble a self-contained pbg_parsimony 3D viewer app for a packed cell.

The "see the results" half of the structural-model investigation: given a
``build_mgen_pack`` result (a ``pack.json`` + sidecar + a ``meshes/`` dir), copy
the bundled ``pbg_parsimony`` viewer next to the pack, rewrite the ingredient
mesh URLs to be app-relative, and drop the pack/sidecar into ``data/`` so the
viewer (WebGL + WebXR/VR) renders the cell. Mirrors ``pbg_parsimony.demo``'s
``build_demo`` assembly, but for an arbitrary pack instead of the demo.

CLI:
    python -m viva_mgen.structural.viewer \
        --pack workspace/studies/s01-maritan-baseline/pack/mgen_maritan_baseline.pack.json \
        --app  workspace/studies/s01-maritan-baseline/viz/viewer
Then: ``python -m http.server -d <app>`` and open the printed URL.
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pbg_parsimony

VIEWER_SRC = Path(pbg_parsimony.__file__).resolve().parent / "viewer"
# All static viewer assets worth shipping (VR helpers included). Tests excluded.
_VIEWER_FILES = ("index.html", "viewer.js", "obj-worker.js", "vr.js", "vr-helpers.js")


def relativize_pack_urls(pack_path: str | Path) -> Path:
    """Rewrite a pack's LOD mesh URLs from absolute paths to ``meshes/<file>``
    (in place), so the workbench's built-in Parsimony Viewer — which serves the
    pack statically from ``<study>/viz/3d/`` and resolves meshes relative to it
    — can fetch each mesh from the sibling ``meshes/`` dir. Idempotent."""
    pack_path = Path(pack_path)
    pack = json.loads(pack_path.read_text())
    changed = False
    for ing in pack.get("ingredients", []):
        for lod in ing.get("shape", {}).get("lods", []):
            rel = "meshes/" + Path(lod["url"]).name
            if lod["url"] != rel:
                lod["url"] = rel
                changed = True
    if changed:
        pack_path.write_text(json.dumps(pack))
    return pack_path


def assemble_viewer_app(pack_path: str | Path, app_dir: str | Path, *,
                        meshes_dir: str | Path | None = None,
                        sidecar_path: str | Path | None = None,
                        link_meshes: bool = True) -> Path:
    """Assemble a self-contained viewer app at ``app_dir`` for ``pack_path``.

    ``meshes_dir`` defaults to ``<pack dir>/meshes``; ``sidecar_path`` defaults
    to the pack's ``*.meta.json`` sibling. The pack is written to the viewer's
    zero-config default path (``data/demo.pack.json`` + ``data/demo.meta.json``),
    so the stock pbg_parsimony viewer loads it with no query string.

    ``link_meshes`` (default) symlinks ``app/meshes`` → the source meshes dir
    instead of copying it (the mesh set is hundreds of MB); pass ``False`` to
    copy each mesh for a fully portable/publishable app. Returns the app dir.
    """
    pack_path = Path(pack_path)
    pack_dir = pack_path.parent
    meshes_dir = Path(meshes_dir) if meshes_dir else pack_dir / "meshes"
    if sidecar_path is None:
        cand = pack_path.with_name(pack_path.name.replace(".pack.json", ".meta.json"))
        sidecar_path = cand if cand.exists() else None

    app = Path(app_dir)
    (app / "data").mkdir(parents=True, exist_ok=True)

    for fn in _VIEWER_FILES:
        src = VIEWER_SRC / fn
        if src.exists():
            shutil.copy(src, app / fn)

    pack = json.loads(pack_path.read_text())
    for ing in pack.get("ingredients", []):
        for lod in ing.get("shape", {}).get("lods", []):
            lod["url"] = "meshes/" + Path(lod["url"]).name

    # Meshes: symlink the whole dir (default) or copy each file (portable).
    app_meshes = app / "meshes"
    if link_meshes:
        if app_meshes.is_symlink() or app_meshes.exists():
            if app_meshes.is_symlink():
                app_meshes.unlink()
        if meshes_dir.exists() and not app_meshes.exists():
            app_meshes.symlink_to(meshes_dir.resolve(), target_is_directory=True)
    else:
        app_meshes.mkdir(exist_ok=True)
        for ing in pack.get("ingredients", []):
            for lod in ing.get("shape", {}).get("lods", []):
                fn = Path(lod["url"]).name
                src = meshes_dir / fn
                if src.exists():
                    shutil.copy(src, app_meshes / fn)

    # Zero-config default the stock viewer loads (data/demo.pack.json + sidecar).
    (app / "data" / "demo.pack.json").write_text(json.dumps(pack))
    if sidecar_path and Path(sidecar_path).exists():
        shutil.copy(sidecar_path, app / "data" / "demo.meta.json")
    return app


def main() -> int:
    ap = argparse.ArgumentParser(description="Assemble a pbg_parsimony viewer app for a pack.")
    ap.add_argument("--pack", required=True, help="path to the *.pack.json")
    ap.add_argument("--app", required=True, help="output app dir")
    a = ap.parse_args()
    app = assemble_viewer_app(a.pack, a.app)
    print(f"viewer app assembled at {app}")
    print(f"serve with:  python -m http.server -d {app}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

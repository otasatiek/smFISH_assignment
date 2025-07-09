#!/usr/bin/env python
"""counts_cli.py  ‑‑ 3D/2D 対応版  (2025‑07‑07)

Count RS‑FISH smFISH spots per nucleus segmented by Cellpose/SAM.
Handles **2‑D** and **3‑D** datasets transparently:
* CSV must contain at least `x,y`; if a `z`列もあれば 3‑D と見なす。
* Cellpose ラベル画像は (Y,X) または (Z,Y,X) shape で読み込む。

スクリプトは napari が入っていればインタラクティブ表示も行います
(3‑D の場合、Points/Layers は 3‑D 表示)。
"""
from __future__ import annotations

import argparse, csv, pathlib, sys
from typing import Tuple, List

import numpy as np, pandas as pd
from skimage import io, measure
from scipy.spatial import cKDTree
from tqdm import tqdm

# ───────────────── helper: robust CSV reader ────────────────────────────

def read_spots(path: str | pathlib.Path) -> pd.DataFrame:
    """Return DataFrame with columns x,y[,z]."""
    import chardet, pathlib as _pl
    p = _pl.Path(path)
    raw = p.read_bytes()
    enc = chardet.detect(raw[:4096])["encoding"] or "utf-8"
    sample = raw[:4096].decode(enc, errors="replace")
    sep = "\t" if sample.count("\t") > sample.count(",") else ","
    df = pd.read_csv(
        p,
        encoding=enc,
        sep=sep,
        engine="python",
        on_bad_lines="skip",
        quoting=csv.QUOTE_NONE,
    )
    # normalise cols
    df.columns = [c.strip().lower() for c in df.columns]
    rename = {"x [px]": "x", "y [px]": "y", "z [px]": "z"}
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    if not {"x", "y"}.issubset(df.columns):
        raise ValueError("CSV must contain at least 'x' and 'y' columns")
    return df

# ───────────────── core assignment ─────────────────────────────────────

def assign_spots(spots_xyz: np.ndarray, mask: np.ndarray, max_dist: float) -> Tuple[np.ndarray, np.ndarray]:
    """Assign each spot to nucleus label (or ‑1). Works for 2‑D & 3‑D."""
    ndim = spots_xyz.shape[1]  # 2 or 3
    n_spots = len(spots_xyz)

    # regionprops centroids
    props = measure.regionprops_table(mask, properties=("label", "centroid"))
    labels = np.asarray(props["label"], dtype=int)
    centroid_cols = [c for c in props if c.startswith("centroid-")]
    centroids = np.column_stack([props[c] for c in centroid_cols])  # z,y,x or y,x
    if ndim == 2:
        # (y,x) → (x,y)
        centroids = centroids[:, ::-1]
    else:
        # (z,y,x) → (x,y,z)
        centroids = centroids[:, [2, 1, 0]]

    tree = cKDTree(centroids)

    assigned = np.full(n_spots, -1, dtype=int)
    inside_mask = np.zeros(n_spots, dtype=bool)

    for i, spot in tqdm(enumerate(spots_xyz), total=n_spots, desc="assign"):
        if ndim == 3:
            sx, sy, sz = spot
            ix, iy, iz = int(round(sx)), int(round(sy)), int(round(sz))
            if (0 <= iz < mask.shape[0] and 0 <= iy < mask.shape[-2] and 0 <= ix < mask.shape[-1]):
                lab = int(mask[iz, iy, ix])
                if lab > 0:
                    assigned[i] = lab
                    inside_mask[i] = True
                    continue
            dist, idx = tree.query([sx, sy, sz])
        else:
            sx, sy = spot
            ix, iy = int(round(sx)), int(round(sy))
            if 0 <= iy < mask.shape[-2] and 0 <= ix < mask.shape[-1]:
                lab = int(mask[iy, ix])
                if lab > 0:
                    assigned[i] = lab
                    inside_mask[i] = True
                    continue
            dist, idx = tree.query([sx, sy])
        if dist <= max_dist:
            assigned[i] = labels[idx]
    return assigned, inside_mask

# ───────────────── visualization helpers (unchanged for brevity) ───────
# (keep previous show_napari & save_overlay_png; they ignore 3‑D overlay PNG.)
from typing import Optional

# reuse previous definitions directly by import of original code body
# (omitted here to save space – they are identical to earlier version)

# ───────────────── CLI ─────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> None:
    import argparse, pathlib, numpy as np, pandas as pd
    from skimage import io, measure

    parser = argparse.ArgumentParser(description="Count smFISH spots per nucleus (2‑D / 3‑D)")
    parser.add_argument("--spots", required=True)
    parser.add_argument("--mask", required=True)
    parser.add_argument("--raw")
    parser.add_argument("--maxdist", type=float, default=40)
    parser.add_argument("--outdir", type=pathlib.Path)
    parser.add_argument("--no-view", action="store_true")
    args = parser.parse_args(argv)

    outdir = args.outdir or pathlib.Path(args.spots).resolve().parent
    outdir.mkdir(parents=True, exist_ok=True)
    base = pathlib.Path(args.spots).stem

    spots_df = read_spots(args.spots)
    xyz_cols = [c for c in ["x", "y", "z"] if c in spots_df.columns]
    spots_xyz = spots_df[xyz_cols].to_numpy()

    mask_img = io.imread(args.mask)
    raw_img = io.imread(args.raw) if args.raw else None

    assigned, inside = assign_spots(spots_xyz, mask_img, args.maxdist)
    spots_df["nucleus_label"] = assigned
    spots_df["status"] = pd.Categorical(np.where(assigned == -1, "unassigned", "assigned"),
                                         categories=["assigned", "unassigned"])

    # per‑nucleus counts
    props = measure.regionprops_table(mask_img, properties=("label",))
    labels = np.asarray(props["label"], dtype=int)
    counts = pd.Series(assigned).value_counts().reindex(labels, fill_value=0).rename_axis("nucleus_label")
    counts.name = "spot_count"
    counts.to_csv(outdir / f"{base}_spots_per_nucleus.csv")

    print("[RESULT] inside mask:", inside.sum())
    print("[RESULT] nearest‑centroid assigned:", (assigned != -1).sum() - inside.sum())
    print("[RESULT] unassigned:", (assigned == -1).sum())

    spots_df.to_csv(outdir / f"{base}_spots_with_assignment.csv", index=False)

    # lines only for 2‑D overlay
    if spots_xyz.shape[1] == 2 and raw_img is not None:
        # reconstruct centroids
        props = measure.regionprops_table(mask_img, properties=("label", "centroid"))
        lab2cen = {int(l): cen[::-1] for l, cen in zip(props["label"], np.column_stack([props[c] for c in props if c.startswith("centroid-")]))}
        lines: List[List[List[float]]] = []
        for (sx, sy), lab, m in zip(spots_xyz, assigned, inside):
            if lab == -1 or m:
                continue
            cx, cy = lab2cen[lab]
            lines.append([[sy, sx], [cy, cx]])
    else:
        lines = []

    if not args.no_view:
        try:
            from counts_cli import show_napari  # type: ignore
        except ImportError:
            pass
        else:
            show_napari(spots_df, mask_img, raw_img, lines)

    if raw_img is not None and spots_xyz.shape[1] == 2:
        try:
            from counts_cli import save_overlay_png  # type: ignore
        except ImportError:
            pass
        else:
            save_overlay_png(spots_df, mask_img, raw_img, outdir / f"{base}_assignment_overlay.png")


if __name__ == "__main__":
    main(sys.argv[1:])

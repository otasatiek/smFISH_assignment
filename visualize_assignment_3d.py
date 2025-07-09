#!/usr/bin/env python
"""
visualize_assignment_3d.py
--------------------------
Show 3-D nuclei mask, smFISH spots, and (optionally) assignment vectors.

Inputs
------
--mask      Cellpose 3-D ラベル画像  (Z,Y,X)  e.g. nuclei_labels.tif
--spotscsv  counts_cli.py が出力した *_spots_with_assignment.csv
--raw       元の smFISH 3-D stack (任意, 背景に表示)
--vectors   True なら核外スポット → 重心への矢印を表示
"""
import argparse, numpy as np, pandas as pd
from skimage import io
import napari

p = argparse.ArgumentParser()
p.add_argument("--mask",      required=True)
p.add_argument("--spotscsv",  required=True)
p.add_argument("--raw")
p.add_argument("--vectors",   action="store_true")
args = p.parse_args()

# ---------- load data ----------
mask = io.imread(args.mask)                     # (Z,Y,X)
spots = pd.read_csv(args.spotscsv)              # x,y,z,nucleus_label,status など

# napari は (zyx) 座標で扱うので列順を揃える
coords = spots[["z", "y", "x"]].to_numpy(float)

# ---------- build viewer ----------
viewer = napari.Viewer(ndisplay=3)

if args.raw:
    raw = io.imread(args.raw)
    viewer.add_image(raw, name="raw", colormap="magenta",
                     blending="additive", opacity=0.6)

viewer.add_labels(mask, name="nuclei", opacity=0.3, rendering="translucent")

# --- Points layer with assignment info ---
viewer.add_points(
    coords,
    name="spots",
    features={
        "status":   spots["status"].to_numpy(),        # assigned / unassigned
        "nucleus":  spots["nucleus_label"].to_numpy()
    },
    face_color="status",
    face_color_cycle={"assigned":"dodgerblue",
                      "unassigned":"gray"},
    border_color        = "status",
    border_color_cycle  = {"assigned":"dodgerblue","unassigned":"gray"},
    size=4
)

# --- optional arrows: spot → centroid ---
if args.vectors:
    # 重心辞書を作る
    from skimage.measure import regionprops_table
    props = regionprops_table(mask, properties=("label", "centroid"))
    label2cen = {lab: np.array([cz, cy, cx])
                 for lab, cz, cy, cx in zip(
                     props["label"],
                     props["centroid-0"], props["centroid-1"], props["centroid-2"])
                 }

    vec_list = []
    for (z, y, x), lab, stat in zip(coords,
                                    spots["nucleus_label"],
                                    spots["status"]):
        if stat == "assigned" and lab in label2cen:
            start = np.array([z, y, x], dtype=float)
            end   = label2cen[lab]
            direction = end - start
            vec_list.append([start, direction])
    vec_arr = np.asarray(vec_list, dtype=float)    # shape (N, 2, 3)

    if vec_list:
        viewer.add_vectors(
        vec_arr,
        edge_color="yellow",
        edge_width=0.4,
        vector_style='arrow',
        name="assign_vectors"        
    )

napari.run()

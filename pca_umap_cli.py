#!/usr/bin/env python
"""pca_umap_cli.py

Run PCA and UMAP on per‑nucleus spot‑count table and plot results.
Colouring priority:
1. If --colour COL is given and COL exists in CSV → use it.
2. Else if a 'pattern' column exists → use it.
3. Else if clustering (KMeans) requested via --k K → colour by cluster ID.
4. Otherwise plot all points in one colour.

Example
-------
python pca_umap_cli.py --csv summarized1.csv --threshold 1 --k 5
"""
import argparse
import subprocess, sys
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns

# ---------- Parse CLI ----------
parser = argparse.ArgumentParser(description="PCA/UMAP visualisation of spot counts")
parser.add_argument("--csv", required=True, help="Input CSV file")
parser.add_argument("--genes", nargs="*", default=None,
                    help="Gene column names (default: auto-detect numeric columns)")
parser.add_argument("--threshold", "-t", type=float, default=None,
                    help="Apply log1p then >=THRESHOLD binarisation to generate a 'pattern' column")
parser.add_argument("--colour", "-c", default=None,
                    help="Existing column name to colour by (overrides others)")
parser.add_argument("--k", type=int, default=None,
                    help="If set, perform KMeans with K clusters and colour by cluster")
parser.add_argument("--out", default=None, help="Output PNG path (default: same dir, pca_umap.png)")
args = parser.parse_args()

# ---------- Load data ----------
path = Path(args.csv)
df = pd.read_csv(path)

# detect gene columns
if args.genes:
    gene_cols = args.genes
else:
    numeric = df.select_dtypes("number").columns
    blacklist = {"nucleus_label", "x", "y", "z"}
    gene_cols = [c for c in numeric if c.lower() not in blacklist]
    if len(gene_cols) < 2:
        parser.error("Cannot find ≥2 gene columns automatically; use --genes to specify.")

# ---------- Preprocess ----------
X = np.log1p(df[gene_cols].to_numpy())
X = StandardScaler().fit_transform(X)

# ---------- Optionally compute 'pattern' (binary combination) ----------
if args.threshold is not None:
    presence = (df[gene_cols] >= args.threshold).astype(int)
    pattern = (presence * (1 << np.arange(len(gene_cols)))).sum(axis=1)
    df["pattern"] = pattern

# ---------- Determine colouring column ----------
colour_col = None
legend_title = None
if args.colour and args.colour in df.columns:
    colour_col = args.colour
    legend_title = colour_col
elif "pattern" in df.columns:
    colour_col = "pattern"
    legend_title = "pattern"
elif args.k:
    labels = KMeans(n_clusters=args.k, random_state=0).fit_predict(X)
    df["cluster"] = labels
    colour_col = "cluster"
    legend_title = f"kmeans (k={args.k})"
else:
    colour_col = None

# ---------- Dimensionality reduction ----------
try:
    import umap
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=0)
    emb = reducer.fit_transform(X)
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "umap-learn"])
    import umap
    emb = umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=0).fit_transform(X)

pca = PCA(n_components=2, random_state=0).fit_transform(X)

# ---------- Plot ----------
fig, axs = plt.subplots(1, 2, figsize=(10, 4))

def scatter(ax, coords, title):
    if colour_col:
        unique = sorted(df[colour_col].unique())
        palette = sns.color_palette("tab10", len(unique))
        lut = dict(zip(unique, palette))
        col = df[colour_col].map(lut)
        ax.scatter(coords[:,0], coords[:,1], c=col, s=8, alpha=0.8)
        # legend
        handles = [plt.Line2D([0],[0],marker='o',color='w',label=str(u),
                              markerfacecolor=lut[u],markersize=6) for u in unique]
        ax.legend(handles=handles,title=legend_title,fontsize=8,frameon=False)
    else:
        ax.scatter(coords[:,0], coords[:,1], c="gray", s=8, alpha=0.7)
    ax.set_title(title)
    ax.set_xlabel("Dim1"); ax.set_ylabel("Dim2")

scatter(axs[0], pca, "PCA")
scatter(axs[1], emb, "UMAP")
plt.tight_layout()

if args.out is None:
    args.out = path.with_suffix("").parent / "pca_umap.png"
plt.savefig(args.out, dpi=300)
print(f"Saved figure to {args.out}")

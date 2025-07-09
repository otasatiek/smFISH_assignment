#!/usr/bin/env python
"""plot_upset.py
Generate an UpSet plot (PNG) **and** a table of co‑expression pattern counts (CSV)
from a per‑nucleus spot‑count table **excluding cells negative for all genes**.

Usage
-----
python plot_upset.py \
    --csv summarized1.csv            # 必須: 入力カウント表
    --threshold 2                    # 発現判定カットオフ (既定 1)
    --genes GeneA GeneB GeneC        # 手動で遺伝子列を指定 (省略可)
    --png   result_upset.png         # 図の保存先 (省略可)
    --counts result_counts.csv       # パターン集計 CSV (省略可)
    --no-show                        # 図をポップアップしない

変更点
-------
* 発現指標行列 `presence` で **全列 False の行を除外**。
  これにより「すべてネガティブ (pattern=0...0)」の細胞は
  UpSet プロットと集計から完全に除外される。
"""
import argparse, sys, subprocess, pathlib
import pandas as pd, matplotlib.pyplot as plt

# ---------- 0. CLI ----------
parser = argparse.ArgumentParser(description="UpSet plot + pattern count exporter (negatives excluded)")
parser.add_argument("--csv", "-i", required=True, help="Input CSV file")
parser.add_argument("--threshold", "-t", type=int, default=1,
                    help="≥ THRESHOLD spots ⇒ expressed (default: 1)")
parser.add_argument("--genes", "-g", nargs="*", default=None,
                    help="Gene column names (default: auto-detect numeric columns)")
parser.add_argument("--png", default=None, help="Output PNG path (default: <csv>_upset.png)")
parser.add_argument("--counts", default=None,
                    help="Output pattern-count CSV path (default: <csv>_pattern_counts.csv)")
parser.add_argument("--no-show", action="store_true", help="Do not open figure window")
args = parser.parse_args()

csv_path  = pathlib.Path(args.csv)
threshold = args.threshold
png_path  = pathlib.Path(args.png) if args.png else csv_path.with_suffix("").with_name(csv_path.stem + "_upset.png")
counts_out = pathlib.Path(args.counts) if args.counts else csv_path.with_suffix("").with_name(csv_path.stem + "_pattern_counts.csv")
gene_cols = args.genes  # None → auto detect later

# ---------- 1. upsetplot import ----------
try:
    import upsetplot
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "upsetplot>=0.8", "-q"])
    import upsetplot

from_func = upsetplot.from_indicators  # ≥0.8 API
UpSet     = upsetplot.UpSet

# ---------- 2. Load & choose gene columns ----------
df = pd.read_csv(csv_path)
if gene_cols is None:
    numeric = df.select_dtypes("number").columns
    blacklist = {"nucleus_label", "x", "y", "z"}
    gene_cols = [c for c in numeric if c.lower() not in blacklist]
    if len(gene_cols) < 2:
        parser.error("Gene columns could not be auto‑detected; use --genes.")

# ---------- 3. Build presence matrix (True/False) ----------
presence = df[gene_cols] >= threshold

# *** NEW ***: drop cells negative for all genes
presence = presence[presence.any(axis=1)]
if presence.empty:
    sys.exit("Error: no cells express any gene under the given threshold.")

# ---------- 4. UpSet data & plot ----------
data = from_func(presence)
up = UpSet(data, subset_size="count", show_percentages=True, orientation="horizontal")
up.plot()
plt.suptitle(f"UpSet (≥{threshold}, negatives excluded)", fontsize=14)
plt.tight_layout()
plt.savefig(png_path, dpi=300)
print(f"[saved] {png_path}")
if not args.no_show:
    plt.show()
else:
    plt.close()

# ---------- 5. Pattern count CSV ----------
pattern_series = presence.apply(lambda r: ''.join('1' if v else '0' for v in r), axis=1)
counts = pattern_series.value_counts().sort_index()
counts.name = "n_cells"
counts.index.name = "pattern (bit order: " + ','.join(gene_cols) + ")"
counts.to_csv(counts_out)
print(f"[saved] {counts_out}")

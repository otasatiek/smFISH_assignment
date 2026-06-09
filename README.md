# smFISH_assignment

`smFISH_assignment` is a small command-line toolkit for assigning RS-FISH-detected smFISH spots to Cellpose/SAM-segmented nuclei, summarizing spot counts per nucleus, and visualizing co-expression patterns in 2D or 3D datasets.

The repository is intended for workflows such as:

1. Detect smFISH spots with RS-FISH and export spot coordinates as CSV.
2. Segment nuclei with Cellpose/SAM and export a label image as TIFF.
3. Assign each spot to the nucleus containing it, or to the nearest nuclear centroid within a distance cutoff.
4. Summarize spot counts per nucleus.
5. Visualize gene-set intersections with Venn diagrams or UpSet plots.
6. Explore per-nucleus expression patterns with PCA/UMAP or napari-based 3D visualization.

## Repository contents

| File | Purpose |
| --- | --- |
| `counts_cli.py` | Assigns smFISH spot coordinates to nuclei in 2D or 3D and exports per-nucleus spot counts plus per-spot assignment results. |
| `plot_upset_nonega.py` | Generates an UpSet plot and a co-expression pattern count table after excluding cells negative for all selected genes. |
| `venn_upset.py` | Generates a Venn diagram for 2 or 3 sets, or attempts an UpSet-style plot when the number of columns exceeds the Venn threshold. |
| `pca_umap_cli.py` | Runs PCA and UMAP on a per-nucleus spot-count table and saves a combined PCA/UMAP figure. |
| `visualize_assignment_3d.py` | Opens a napari 3D viewer showing the nuclei mask, assigned spots, and optionally spot-to-centroid assignment vectors. |

## Installation

Use a Python environment with scientific image-analysis packages installed. Python 3.10 or later is recommended.

```bash
conda create -n smfish-assignment python=3.10 -y
conda activate smfish-assignment
pip install pandas numpy scipy scikit-image tqdm matplotlib matplotlib-venn upsetplot scikit-learn seaborn umap-learn chardet
```

For interactive 3D visualization, also install napari:

```bash
pip install "napari[all]"
```

## Input data

### Spot CSV

`counts_cli.py` expects a CSV or TSV file exported from RS-FISH or a similar spot detector.

Required coordinate columns:

- 2D data: `x`, `y`
- 3D data: `x`, `y`, `z`

The script also accepts RS-FISH-style column names and normalizes them automatically:

- `x [px]` → `x`
- `y [px]` → `y`
- `z [px]` → `z`

### Nucleus mask

The nucleus mask should be a labeled TIFF image, typically exported from Cellpose/SAM or equivalent segmentation software.

Expected dimensions:

- 2D: `(Y, X)`
- 3D: `(Z, Y, X)`

Background should be label `0`. Each nucleus should have a positive integer label.

### Raw image

A raw image stack can optionally be supplied for visualization. It is not required for the assignment itself.

## Basic workflow

### 1. Assign smFISH spots to nuclei

```bash
python counts_cli.py \
  --spots path/to/rsfish_spots.csv \
  --mask path/to/nuclei_labels.tif \
  --maxdist 40 \
  --outdir results \
  --no-view
```

Optional raw image:

```bash
python counts_cli.py \
  --spots path/to/rsfish_spots.csv \
  --mask path/to/nuclei_labels.tif \
  --raw path/to/raw_image.tif \
  --maxdist 40 \
  --outdir results \
  --no-view
```

Assignment rule:

1. If a spot coordinate falls inside a labeled nucleus, that nuclear label is assigned directly.
2. If the spot is outside all nuclei, the nearest nuclear centroid is queried.
3. If the nearest centroid is within `--maxdist` pixels, that nucleus is assigned.
4. Otherwise, the spot is marked as unassigned with `nucleus_label = -1`.

Outputs:

| Output | Description |
| --- | --- |
| `<spot_csv_stem>_spots_per_nucleus.csv` | Per-nucleus spot count table. |
| `<spot_csv_stem>_spots_with_assignment.csv` | Original spot table plus `nucleus_label` and `status` columns. |

Example output columns in `*_spots_with_assignment.csv`:

```text
x,y,z,nucleus_label,status
```

`status` is either `assigned` or `unassigned`.

## Combining multiple gene channels

For multi-gene analysis, run `counts_cli.py` separately for each smFISH channel, then merge the per-nucleus count tables by `nucleus_label`.

A merged count table should look like this:

```text
nucleus_label,GeneA_spot_count,GeneB_spot_count,GeneC_spot_count
1,3,0,5
2,0,2,1
3,4,4,0
```

The downstream scripts generally treat numeric columns other than `nucleus_label`, `x`, `y`, and `z` as gene/count columns.

## Co-expression visualization

### UpSet plot excluding all-negative cells

Use `plot_upset_nonega.py` when you want to remove nuclei that are negative for all selected genes before plotting co-expression patterns.

```bash
python plot_upset_nonega.py \
  --csv results/merged_counts.csv \
  --threshold 2 \
  --genes GeneA_spot_count GeneB_spot_count GeneC_spot_count \
  --png results/upset.png \
  --counts results/pattern_counts.csv \
  --no-show
```

If `--genes` is omitted, numeric columns are auto-detected after excluding common coordinate/label columns.

Outputs:

| Output | Description |
| --- | --- |
| `*_upset.png` | UpSet plot. |
| `*_pattern_counts.csv` | Binary co-expression pattern counts. |

The pattern bit order is written into the index name of the output CSV.

### Venn diagram or UpSet-style plot

Use `venn_upset.py` for quick Venn diagrams when analyzing 2 or 3 expression sets.

```bash
python venn_upset.py \
  --csv results/merged_counts.csv \
  --output results/venn_or_upset.png \
  --out_csv results/coexpression_counts.csv \
  --bin_threshold 2 \
  --venn_threshold 4
```

Interpretation:

- Columns other than `nucleus_label` are binarized.
- Values `>= --bin_threshold` are treated as positive.
- If the number of target columns is 2 or 3 and below `--venn_threshold`, a Venn diagram is generated.
- For larger numbers of columns, use `plot_upset_nonega.py` as the more reliable UpSet workflow.

## PCA and UMAP

Use `pca_umap_cli.py` to visualize per-nucleus expression profiles in reduced dimensions.

```bash
python pca_umap_cli.py \
  --csv results/merged_counts.csv \
  --genes GeneA_spot_count GeneB_spot_count GeneC_spot_count \
  --threshold 1 \
  --k 5 \
  --out results/pca_umap.png
```

Processing steps:

1. Select gene/count columns.
2. Apply `log1p` transformation.
3. Standardize features with `StandardScaler`.
4. Compute PCA and UMAP coordinates.
5. Color points by one of the following, in priority order:
   - a user-specified `--colour` column,
   - a generated `pattern` column if `--threshold` is given,
   - KMeans cluster ID if `--k` is given,
   - a single default color if none of the above applies.

## 3D visualization in napari

After running `counts_cli.py` on 3D data, inspect the assignment interactively:

```bash
python visualize_assignment_3d.py \
  --mask path/to/nuclei_labels_3d.tif \
  --spotscsv results/rsfish_spots_with_assignment.csv \
  --raw path/to/raw_stack.tif \
  --vectors
```

The viewer displays:

- raw smFISH image stack, if supplied,
- labeled nuclei,
- assigned and unassigned spots,
- optional assignment vectors from each assigned spot to the centroid of the assigned nucleus.

## Suggested end-to-end example

```bash
# 1. Assign spots for each gene/channel
python counts_cli.py --spots GeneA_spots.csv --mask nuclei_labels.tif --outdir results --maxdist 40 --no-view
python counts_cli.py --spots GeneB_spots.csv --mask nuclei_labels.tif --outdir results --maxdist 40 --no-view
python counts_cli.py --spots GeneC_spots.csv --mask nuclei_labels.tif --outdir results --maxdist 40 --no-view

# 2. Merge *_spots_per_nucleus.csv files by nucleus_label outside this toolkit
#    Example merged file: results/merged_counts.csv

# 3. Plot co-expression patterns
python plot_upset_nonega.py --csv results/merged_counts.csv --threshold 2 --png results/upset.png --counts results/pattern_counts.csv --no-show

# 4. Explore expression profiles
python pca_umap_cli.py --csv results/merged_counts.csv --threshold 1 --k 5 --out results/pca_umap.png

# 5. Inspect 3D assignment if using 3D images
python visualize_assignment_3d.py --mask nuclei_labels_3d.tif --spotscsv results/GeneA_spots_with_assignment.csv --raw GeneA_raw_stack.tif --vectors
```

## Notes and current limitations

- `counts_cli.py` currently contains placeholders for some 2D napari/overlay helper functions. The command-line assignment and CSV export are the primary supported functions.
- `visualize_assignment_3d.py` assumes the spot assignment table contains `z`, `y`, `x`, `nucleus_label`, and `status` columns.
- `plot_upset_nonega.py` excludes nuclei that are negative for all selected genes. This is useful when focusing on positive co-expression patterns but should be considered when calculating absolute frequencies across all segmented nuclei.
- `pca_umap_cli.py` may try to install `umap-learn` automatically if it is missing. For reproducible environments, install all dependencies explicitly before running.
- Coordinate convention differs between tools: spot tables use `x,y,z`, while napari/image stacks use `z,y,x`. The scripts convert these internally where needed.


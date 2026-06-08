if "snakemake" in locals():
    coembedding_path = snakemake.input["csv"]
    rna_annotations_path = snakemake.input["rna_annotations"]
    atac_h5ad_path = snakemake.input["atac_h5ad"]
    plot_path = snakemake.output["plot"]
    coords_path = snakemake.output["coords"]
    annotated_plot_path = snakemake.output["annotated_plot"]
else:
    coembedding_path = "results/glue/coembedding.csv"
    rna_annotations_path = "config/biopsydata_withkpmp_final_annotations(in).csv"
    atac_h5ad_path = "data/ATAC/peak_mat2.h5ad"
    plot_path = "plots/glue/coembedding_umap.pdf"
    coords_path = "results/glue/coembedding_umap.csv"
    annotated_plot_path = "plots/glue/coembedding_umap_annotations.pdf"

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


UNKNOWN_LABEL = "Unannotated"
PALETTE_21 = [
    "#e6194B", "#4363d8", "#ffe119", "#3cb44b", "#f58231", "#911eb4",
    "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990", "#dcbeff",
    "#9A6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",
    "#000075", "#a9a9a9", "#000000",
]
DISTINCT_COLORS = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#999999", "#A65628", "#FF69B4",
    "#87CEEB", "#228B22", "#FFD700", "#6A5ACD", "#DC143C",
    "#00CED1", "#8B0000", "#B8860B", "#9400D3", "#40E0D0",
]
UNANNOTATED_COLOR = "#D3D3D3"
RNA_ANNOTATION_POINT_SIZE = 0.275
ATAC_ANNOTATION_POINT_SIZE = 0.75


def parse_modality(index):
    modalities = []
    cell_ids = []
    for value in index.astype(str):
        if value.startswith("rna_"):
            modalities.append("RNA")
            cell_ids.append(value.removeprefix("rna_"))
        elif value.startswith("atac_"):
            modalities.append("ATAC")
            cell_ids.append(value.removeprefix("atac_"))
        else:
            raise ValueError(
                "Coembedding row names must start with 'rna_' or 'atac_', "
                f"but found: {value}"
            )
    return modalities, cell_ids


def add_rna_annotations(coords):
    annotations = pd.read_csv(rna_annotations_path)
    required_columns = {"barcode_combined_id", "main_class"}
    missing = required_columns.difference(annotations.columns)
    if missing:
        raise ValueError(
            "RNA annotation table is missing required columns: "
            f"{', '.join(sorted(missing))}."
        )
    if annotations["barcode_combined_id"].duplicated().any():
        raise ValueError("RNA annotation table contains duplicated cell IDs.")

    annotation_map = annotations.set_index("barcode_combined_id")["main_class"]
    rna_mask = coords["modality"] == "RNA"
    coords.loc[rna_mask, "rna_annotation"] = (
        coords.loc[rna_mask, "cellid"].map(annotation_map).fillna(UNKNOWN_LABEL)
    )
    missing_count = int((coords.loc[rna_mask, "rna_annotation"] == UNKNOWN_LABEL).sum())
    if missing_count:
        print(f"RNA annotations missing for {missing_count} coembedded RNA cells.")
    return coords


def add_atac_annotations(coords):
    atac = sc.read_h5ad(atac_h5ad_path, backed="r")
    if "cluster_label" in atac.obs.columns:
        label_column = "cluster_label"
    elif "cluster" in atac.obs.columns:
        label_column = "cluster"
        print(
            "ATAC obs does not contain 'cluster_label'; using existing "
            "'cluster' column for the ATAC annotation panel."
        )
    else:
        raise ValueError(
            "ATAC AnnData must contain 'cluster_label' or fallback column 'cluster'."
        )

    annotation_map = atac.obs[label_column].astype(str)
    if annotation_map.index.has_duplicates:
        raise ValueError("ATAC annotation table contains duplicated cell IDs.")

    atac_mask = coords["modality"] == "ATAC"
    coords.loc[atac_mask, "atac_annotation"] = (
        coords.loc[atac_mask, "cellid"].map(annotation_map).fillna(UNKNOWN_LABEL)
    )
    missing_count = int(
        (coords.loc[atac_mask, "atac_annotation"] == UNKNOWN_LABEL).sum()
    )
    if missing_count:
        print(f"ATAC annotations missing for {missing_count} coembedded ATAC cells.")
    return coords, label_column


def main_class_sort_key(label):
    try:
        return (0, int(str(label).split("_", 1)[0]), str(label))
    except ValueError:
        return (1, str(label))


def color_map(labels, palette, include_unannotated=False, numeric_prefix_order=False):
    levels = sorted(pd.Index(labels.dropna().astype(str)).unique())
    if numeric_prefix_order:
        levels = sorted(levels, key=main_class_sort_key)
    if include_unannotated and UNKNOWN_LABEL in levels:
        levels = [level for level in levels if level != UNKNOWN_LABEL] + [UNKNOWN_LABEL]
    if len(levels) > len(palette) + int(include_unannotated):
        raise ValueError(
            f"Need colors for {len(levels)} labels, but only {len(palette)} "
            "palette colors are available."
        )

    colors = {}
    palette_index = 0
    for level in levels:
        if include_unannotated and level == UNKNOWN_LABEL:
            colors[level] = UNANNOTATED_COLOR
        else:
            colors[level] = palette[palette_index]
            palette_index += 1
    return colors


def scatter_by_label(
    ax,
    data,
    label_column,
    title,
    palette,
    point_size,
    include_unannotated=False,
    numeric_prefix_order=False,
):
    colors = color_map(
        data[label_column],
        palette=palette,
        include_unannotated=include_unannotated,
        numeric_prefix_order=numeric_prefix_order,
    )
    for label in colors:
        group = data.loc[data[label_column].astype(str) == label]
        if group.empty:
            continue
        ax.scatter(
            group["UMAP1"],
            group["UMAP2"],
            s=point_size,
            alpha=0.65,
            c=[colors[str(label)]],
            label=f"{label} (n={len(group):,})",
            linewidths=0,
            rasterized=True,
        )
    ax.set_title(title)
    ax.set_xlabel("GLUE UMAP 1")
    ax.set_ylabel("GLUE UMAP 2")
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0,
        frameon=False,
        markerscale=8,
        fontsize=6,
    )


def main():
    coembedding = pd.read_csv(coembedding_path, index_col=0)
    if coembedding.empty:
        raise ValueError("Coembedding table is empty.")
    if coembedding.index.has_duplicates:
        raise ValueError("Coembedding table contains duplicated row names.")

    modalities, cell_ids = parse_modality(coembedding.index)
    adata = ad.AnnData(
        X=coembedding.to_numpy(dtype=np.float32),
        obs=pd.DataFrame(
            {"modality": modalities, "cellid": cell_ids},
            index=coembedding.index,
        ),
    )

    sc.pp.neighbors(adata, n_neighbors=15, use_rep="X", random_state=0)
    sc.tl.umap(adata, min_dist=0.5, random_state=0)

    coords = pd.DataFrame(
        {
            "cellid": adata.obs["cellid"].to_numpy(),
            "modality": adata.obs["modality"].to_numpy(),
            "UMAP1": adata.obsm["X_umap"][:, 0],
            "UMAP2": adata.obsm["X_umap"][:, 1],
        },
        index=adata.obs_names,
    )
    coords = add_rna_annotations(coords)
    coords, atac_label_column = add_atac_annotations(coords)
    coords.to_csv(coords_path, index_label="coembedding_id")

    colors = {"RNA": "#2F6B4F", "ATAC": "#D8843A"}
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    for modality in ["RNA", "ATAC"]:
        mask = coords["modality"] == modality
        ax.scatter(
            coords.loc[mask, "UMAP1"],
            coords.loc[mask, "UMAP2"],
            s=0.2,
            alpha=0.45,
            c=colors[modality],
            label=f"{modality} (n={int(mask.sum()):,})",
            linewidths=0,
            rasterized=True,
        )

    ax.set_xlabel("GLUE UMAP 1")
    ax.set_ylabel("GLUE UMAP 2")
    ax.set_title("RNA and ATAC GLUE coembedding")
    ax.legend(markerscale=10, frameon=False)
    ax.set_aspect("equal", adjustable="datalim")
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)

    annotated = coords.loc[
        coords["rna_annotation"].notna() | coords["atac_annotation"].notna()
    ].copy()
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), constrained_layout=True)
    scatter_by_label(
        axes[0],
        annotated.loc[annotated["modality"] == "RNA"],
        "rna_annotation",
        "RNA coembedding by main_class",
        palette=PALETTE_21,
        point_size=RNA_ANNOTATION_POINT_SIZE,
        include_unannotated=True,
        numeric_prefix_order=True,
    )
    scatter_by_label(
        axes[1],
        annotated.loc[annotated["modality"] == "ATAC"],
        "atac_annotation",
        f"ATAC coembedding by {atac_label_column}",
        palette=DISTINCT_COLORS,
        point_size=ATAC_ANNOTATION_POINT_SIZE,
    )
    fig.savefig(annotated_plot_path, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()

if "snakemake" in locals():
    coembedding_path = snakemake.input["csv"]
    rna_annotations_path = snakemake.input["rna_annotations"]
    atac_annotations_path = snakemake.input["atac_annotations"]
    clusters_path = snakemake.input["clusters"]
    glue_annotations_path = snakemake.input["glue_annotations"]
    plot_path = snakemake.output["plot"]
    coords_path = snakemake.output["coords"]
    annotated_plot_path = snakemake.output["annotated_plot"]
else:
    coembedding_path = "results/glue/coembedding.csv"
    rna_annotations_path = "results/seurat/subset_high_quality/metadata_annotations.csv"
    atac_annotations_path = "data/ATAC/label_transfer_RNA_ATAC_V3_clean.csv"
    clusters_path = "results/glue/coembedding_clusters.csv"
    glue_annotations_path = "config/GLUE_annotations.csv"
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
SHARED_ANNOTATION_POINT_SIZE = 0.2
GLUE_CLUSTER_COLUMN = "leiden_res_0_7"
GLUE_ANNOTATION_COLORS = {
    "PT": "#0072B2",
    "aPT": "#D55E00",
    "dPT": "#009E73",
    "TAL": "#E69F00",
    "DCT": "#CC79A7",
    "CNT_PC": "#56B4E9",
    "IC": "#F0E442",
    "EC": "#7F3C8D",
    "Immune": "#11A579",
    "Stromal": "#3969AC",
    "POD": "#F2B701",
    "PEC": "#E73F74",
}
ATAC_CLUSTER_ORDER = [
    "atlas_C1",
    "atlas_C12",
    "atlas_C14",
    "atlas_C3",
    "atlas_C4",
    "pc_cnt_dct_C1",
    "pc_cnt_dct_C2",
    "pt_C1",
    "pt_C2",
    "pt_C3",
    "pt_C4",
    "pt_C5",
    "pt_C6",
    "pt_C8",
    "vsm_p_fib_C1",
    "vsm_p_fib_C2",
    "vsm_p_fib_C4",
    "vsm_p_fib_C5",
]
ATAC_CLUSTER_COLORS = dict(zip(ATAC_CLUSTER_ORDER, DISTINCT_COLORS))


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


def require_columns(data, columns, label):
    missing = set(columns).difference(data.columns)
    if missing:
        raise ValueError(
            f"{label} is missing required columns: "
            f"{', '.join(sorted(missing))}."
        )


def add_rna_annotations(coords):
    annotations = pd.read_csv(rna_annotations_path)
    require_columns(annotations, {"cellid", "level2"}, "RNA annotation table")
    if annotations["cellid"].duplicated().any():
        raise ValueError("RNA annotation table contains duplicated cell IDs.")

    annotation_map = annotations.set_index("cellid")["level2"]
    rna_mask = coords["modality"] == "RNA"
    coords.loc[rna_mask, "rna_annotation"] = (
        coords.loc[rna_mask, "cellid"].map(annotation_map).fillna(UNKNOWN_LABEL)
    )
    missing_count = int((coords.loc[rna_mask, "rna_annotation"] == UNKNOWN_LABEL).sum())
    if missing_count:
        print(f"RNA annotations missing for {missing_count} coembedded RNA cells.")
    return coords


def add_atac_annotations(coords):
    annotations = pd.read_csv(atac_annotations_path, index_col=0)
    require_columns(
        annotations,
        {"cluster_label", "cluster"},
        "ATAC annotation table",
    )
    if annotations.index.has_duplicates:
        raise ValueError("ATAC annotation table contains duplicated cell IDs.")
    annotations = annotations.loc[:, ["cluster_label", "cluster"]].astype(str)

    atac_mask = coords["modality"] == "ATAC"
    coords.loc[atac_mask, "atac_cluster_label"] = (
        coords.loc[atac_mask, "cellid"]
        .map(annotations["cluster_label"])
        .fillna(UNKNOWN_LABEL)
    )
    coords.loc[atac_mask, "atac_cluster"] = (
        coords.loc[atac_mask, "cellid"]
        .map(annotations["cluster"])
        .fillna(UNKNOWN_LABEL)
    )
    missing_count = int(
        (coords.loc[atac_mask, "atac_cluster_label"] == UNKNOWN_LABEL).sum()
    )
    if missing_count:
        print(
            f"ATAC cluster_label annotations missing for "
            f"{missing_count} coembedded ATAC cells."
        )
    missing_count = int((coords.loc[atac_mask, "atac_cluster"] == UNKNOWN_LABEL).sum())
    if missing_count:
        print(
            f"ATAC cluster annotations missing for "
            f"{missing_count} coembedded ATAC cells."
        )
    return coords


def add_glue_annotations(coords):
    clusters = pd.read_csv(clusters_path)
    required_columns = {"coembedding_id", GLUE_CLUSTER_COLUMN}
    missing = required_columns.difference(clusters.columns)
    if missing:
        raise ValueError(
            "Coembedding clusters are missing required columns: "
            f"{', '.join(sorted(missing))}."
        )
    if clusters["coembedding_id"].duplicated().any():
        raise ValueError("Coembedding clusters contain duplicated coembedding IDs.")

    annotations = pd.read_csv(glue_annotations_path)
    required_columns = {GLUE_CLUSTER_COLUMN, "annotation"}
    missing = required_columns.difference(annotations.columns)
    if missing:
        raise ValueError(
            "GLUE annotations are missing required columns: "
            f"{', '.join(sorted(missing))}."
        )
    annotations[GLUE_CLUSTER_COLUMN] = annotations[GLUE_CLUSTER_COLUMN].astype(str)
    annotations["annotation"] = annotations["annotation"].astype(str)
    if annotations[GLUE_CLUSTER_COLUMN].duplicated().any():
        raise ValueError("GLUE annotations contain duplicated Leiden cluster IDs.")

    cluster_map = clusters.set_index("coembedding_id")[GLUE_CLUSTER_COLUMN].astype(str)
    annotation_map = annotations.set_index(GLUE_CLUSTER_COLUMN)["annotation"]
    observed_clusters = set(cluster_map.dropna())
    missing_annotations = sorted(observed_clusters.difference(annotation_map.index))
    if missing_annotations:
        raise ValueError(
            "GLUE annotations are missing Leiden cluster IDs: "
            f"{', '.join(missing_annotations)}."
        )

    coords["glue_cluster"] = coords.index.astype(str).map(cluster_map)
    coords["glue_annotation"] = (
        coords["glue_cluster"].map(annotation_map).fillna(UNKNOWN_LABEL)
    )
    missing_count = int((coords["glue_annotation"] == UNKNOWN_LABEL).sum())
    if missing_count:
        print(f"GLUE annotations missing for {missing_count} coembedded cells.")
    return coords


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


def configured_color_map(labels, colors, include_unannotated=False):
    levels = sorted(pd.Index(labels.dropna().astype(str)).unique())
    if include_unannotated and UNKNOWN_LABEL in levels:
        levels = [level for level in levels if level != UNKNOWN_LABEL] + [UNKNOWN_LABEL]
    missing = [
        level for level in levels
        if not (include_unannotated and level == UNKNOWN_LABEL) and level not in colors
    ]
    if missing:
        raise ValueError(
            f"Missing configured colors for labels: {', '.join(sorted(missing))}."
        )

    mapped = {}
    for level in levels:
        if include_unannotated and level == UNKNOWN_LABEL:
            mapped[level] = UNANNOTATED_COLOR
        else:
            mapped[level] = colors[level]
    return mapped


def scatter_by_label(
    ax,
    data,
    label_column,
    title,
    palette,
    point_size,
    include_unannotated=False,
    numeric_prefix_order=False,
    colors=None,
):
    if colors is None:
        colors = color_map(
            data[label_column],
            palette=palette,
            include_unannotated=include_unannotated,
            numeric_prefix_order=numeric_prefix_order,
        )
    else:
        colors = configured_color_map(
            data[label_column],
            colors=colors,
            include_unannotated=include_unannotated,
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
    ax.set_xticks([])
    ax.set_yticks([])
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
    coords = add_atac_annotations(coords)
    coords = add_glue_annotations(coords)
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
    ax.set_xticks([])
    ax.set_yticks([])
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)

    annotated = coords.loc[
        coords["rna_annotation"].notna()
        | coords["atac_cluster_label"].notna()
        | coords["glue_annotation"].notna()
        | coords["atac_cluster"].notna()
    ].copy()
    fig, axes = plt.subplots(2, 2, figsize=(16, 12), constrained_layout=True)
    scatter_by_label(
        axes[0, 0],
        annotated.loc[annotated["modality"] == "RNA"],
        "rna_annotation",
        "RNA coembedding by level2",
        palette=PALETTE_21,
        point_size=RNA_ANNOTATION_POINT_SIZE,
        include_unannotated=True,
    )
    scatter_by_label(
        axes[0, 1],
        annotated.loc[annotated["modality"] == "ATAC"],
        "atac_cluster_label",
        "ATAC coembedding by cluster_label",
        palette=DISTINCT_COLORS,
        point_size=ATAC_ANNOTATION_POINT_SIZE,
        include_unannotated=True,
    )
    scatter_by_label(
        axes[1, 0],
        annotated,
        "glue_annotation",
        "Shared coembedding by GLUE annotation",
        palette=None,
        point_size=SHARED_ANNOTATION_POINT_SIZE,
        include_unannotated=True,
        colors=GLUE_ANNOTATION_COLORS,
    )
    scatter_by_label(
        axes[1, 1],
        annotated.loc[annotated["modality"] == "ATAC"],
        "atac_cluster",
        "ATAC coembedding by cluster",
        palette=None,
        point_size=ATAC_ANNOTATION_POINT_SIZE,
        include_unannotated=True,
        colors=ATAC_CLUSTER_COLORS,
    )
    fig.savefig(annotated_plot_path, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()

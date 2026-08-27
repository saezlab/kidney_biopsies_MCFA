MARKER_ULM_INPUTS = {
    "level1": "results/seurat/subset_high_quality/markers/findallmarkers_level1.csv",
    "level2": "results/seurat/subset_high_quality/markers/findallmarkers_level2.csv",
}

MARKER_ULM_LABELS = {
    "level1": "Level1 marker signature",
    "level2": "Level2 marker signature",
}

MARKER_ULM_LEVELS = "|".join(MARKER_ULM_INPUTS.keys())


rule atac_marker_ulm:
    """Score ATAC gene scores against RNA marker signatures with ULM."""
    input:
        h5ad="data/ATAC/gene_score_matrix.h5ad",
        markers_csv=lambda wildcards: MARKER_ULM_INPUTS[wildcards.level]
    output:
        scores="results/ATAC/{level}_marker_ulm_scores.csv",
        pvals="results/ATAC/{level}_marker_ulm_pvals.csv"
    wildcard_constraints:
        level=MARKER_ULM_LEVELS
    params:
        signature_label=lambda wildcards: wildcards.level
    container:
        config["containers"]["scanpy"]
    script:
        "../scripts/ATAC/marker_ulm.py"


rule plot_atac_marker_ulm_umap:
    """Plot ATAC UMAPs colored by marker ULM enrichment scores."""
    input:
        scores="results/ATAC/{level}_marker_ulm_scores.csv",
        pvals="results/ATAC/{level}_marker_ulm_pvals.csv",
        coords="data/ATAC/coordinates.csv"
    output:
        plot="plots/ATAC/{level}_marker_ulm_umap.pdf"
    wildcard_constraints:
        level=MARKER_ULM_LEVELS
    params:
        signature_label=lambda wildcards: MARKER_ULM_LABELS[wildcards.level]
    container:
        config["containers"]["scanpy"]
    script:
        "../scripts/ATAC/plot_marker_ulm_umap.py"


rule summarize_atac_marker_ulm_by_cluster:
    """Summarize ATAC marker ULM scores by ATAC cluster."""
    input:
        scores="results/ATAC/{level}_marker_ulm_scores.csv",
        annotations="data/ATAC/label_transfer_RNA_ATAC_V3_clean.csv"
    output:
        cell_zscores="results/ATAC/{level}_marker_ulm_cell_zscores.csv",
        cluster_mean_scores="results/ATAC/{level}_marker_ulm_cluster_mean_scores.csv"
    wildcard_constraints:
        level=MARKER_ULM_LEVELS
    params:
        signature_label=lambda wildcards: wildcards.level
    container:
        config["containers"]["scanpy"]
    script:
        "../scripts/ATAC/summarize_marker_ulm_by_cluster.py"


rule plot_atac_marker_ulm_cluster_heatmap:
    """Plot cluster-level heatmaps for ATAC marker ULM mean scores."""
    input:
        cluster_mean_scores="results/ATAC/{level}_marker_ulm_cluster_mean_scores.csv"
    output:
        heatmap="plots/ATAC/{level}_marker_ulm_cluster_mean_scores_heatmap.pdf"
    wildcard_constraints:
        level=MARKER_ULM_LEVELS
    params:
        signature_label=lambda wildcards: MARKER_ULM_LABELS[wildcards.level]
    container:
        config["containers"]["scanpy"]
    script:
        "../scripts/ATAC/plot_marker_ulm_cluster_heatmap.py"


rule PT_aPT_differential_accessibility:
    """Compare pt_C6 (aPT) against pt_C3 (PT) differential accessibility."""
    input:
        project_path="data/ATAC/proj_peaks2"
    output:
        marker_peaks="results/ATAC/PT_aPT_differential_accessibility_marker_peaks.rds",
        apt_open_enrich_motifs="results/ATAC/aPT_open_motif_enrichment.rds",
        pt_open_enrich_motifs="results/ATAC/PT_open_motif_enrichment.rds"
    params:
        group_by="cluster",
        use_group="pt_C6",
        background_group="pt_C3"
    threads: 24
    container:
        config["containers"]["archr"]
    script:
        "../scripts/ATAC/PT_aPT_differential_accessibility.R"


rule motif_peak_matches:
    """Map enriched motifs to direction-matched PT/aPT differentially accessible peaks."""
    input:
        project_path="data/ATAC/proj_peaks2",
        marker_peaks="results/ATAC/PT_aPT_differential_accessibility_marker_peaks.rds",
        apt_open_enrich_motifs="results/ATAC/aPT_open_motif_enrichment.rds",
        pt_open_enrich_motifs="results/ATAC/PT_open_motif_enrichment.rds"
    output:
        rds="results/ATAC/motif_peak_matches.rds",
        csv="results/ATAC/motif_peak_matches.csv"
    threads: 1
    container:
        config["containers"]["archr"]
    script:
        "../scripts/ATAC/motif_peak_matches.R"


rule export_PT_aPT_differential_accessibility_marker_peaks:
    """Export PT/aPT differential accessibility statistics with peak coordinates."""
    input:
        project_path="data/ATAC/proj_peaks2",
        marker_peaks="results/ATAC/PT_aPT_differential_accessibility_marker_peaks.rds"
    output:
        csv="results/ATAC/PT_aPT_differential_accessibility_marker_peaks.csv"
    threads: 1
    container:
        config["containers"]["archr"]
    script:
        "../scripts/ATAC/export_pt_apt_differential_accessibility_marker_peaks.R"


rule export_rela_motif_peak_matches:
    """Export ArchR peaks with RELA motif matches."""
    input:
        project_path="data/ATAC/proj_peaks2"
    output:
        csv="results/ATAC/rela_motif_peak_matches.csv"
    params:
        motif_prefix="RELA"
    threads: 1
    container:
        config["containers"]["archr"]
    script:
        "../scripts/ATAC/export_rela_motif_peak_matches.R"


rule export_ap1_motif_peak_matches:
    """Export ArchR peaks with AP-1 family motif matches."""
    input:
        project_path="data/ATAC/proj_peaks2"
    output:
        csv="results/ATAC/ap1_motif_peak_matches.csv"
    params:
        motif_prefixes=["JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "BATF"]
    threads: 1
    container:
        config["containers"]["archr"]
    script:
        "../scripts/ATAC/export_ap1_motif_peak_matches.R"

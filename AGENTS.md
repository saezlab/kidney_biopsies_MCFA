## Snakemake

- Any Snakemake behaviour must be compatible with Snakemake 9, and updated if necessary.
- Use the standard script interface patterns:
  - R scripts should use `if (exists("snakemake")) { ... } else { ... }`.
  - Python scripts should use `if "snakemake" in locals(): ... else: ...`.
  - In the `if` branch, read all `snakemake.input`, `snakemake.output`, and `snakemake.params`.
  - In the `else` branch, keep local defaults for direct script execution.
  - This workflow assumes Snakemake handles creation of parent directories for declared outputs; do not add manual parent-directory creation in scripts just to satisfy declared outputs.
  - For Seurat/SeuratObject v5 code, do not use `slot=` with `GetAssayData()`; use `layer=` instead (for example `layer = "counts"`).
  - For merged Seurat v5 assays, `layer=` alone may still fail if the assay contains multiple per-sample layers. Inspect with `Layers(...)` and, when needed, join relevant layers with `JoinLayers(...)` before extracting data with `LayerData(...)`.
- For Python scripts run via Snakemake `script:`, do **not** use `from __future__ import ...`.
  - Snakemake prepends runtime code, so future imports may no longer be the first statement and can fail with `SyntaxError`.
- When suggesting workflow execution, always provide exact copy-paste Snakemake commands.
- Prefer explicit targets and include key flags: `--cores` and `--software-deployment-method apptainer`.
- Command format defaults:
  - Dry-run: `snakemake --cores 1 <target> --software-deployment-method apptainer --rerun-triggers mtime -n`
  - Run: `snakemake --cores 1 <target> --software-deployment-method apptainer --rerun-triggers mtime`
- Prefer `--rerun-triggers mtime` for workflow execution in this repository.
- Do not use the cache/home wrapper variant by default in suggested commands.
- Only provide the cache/home wrapper variant if the user explicitly asks for it or if plain `snakemake ...` fails in this environment.
- When touching containerized rules or related config, prefer container references that Snakemake/Apptainer can pull directly from Docker Hub.
  - Prefer `container:` with values sourced from `config["containers"][...]` for new or touched rules.
  - Do not introduce new machine-specific absolute `.sif` paths as the canonical container source.
  - Exception: the rule `run_glue` in `workflow/rules/glue.smk` should prefer `conda:` with `workflow/envs/glue.yaml` instead of containerized execution, because GLUE runs are expected to use GPU and should avoid being coupled to `snakemake` + `singularity/apptainer`.
- Keep Snakemake naming simple when reuse is likely:
  - Prefer generic I/O key names like `rds`, `csv`, `plot`, `degs_plot`, `top_markers_plot`, `markers_csv`.
  - Avoid branch-specific key names (for example `subset_rds`, `annotation_rds`, `mesimm_degs_plot`) unless the rule/script truly needs multiple variants of the same type.
  - Reusable scripts should accept simple keys first (for example `rds`) and only keep legacy aliases temporarily for backward compatibility.
- Keep Snakemake rule files easy to read:
  - Prefer direct literal paths in `input:` and `output:` when a path is owned by a single rule.
  - Do not create top-of-file path constants for one-off outputs just to avoid repeating a literal once.
  - Centralize only config-backed values, shared wildcard expansions, or aggregates reused across multiple rules or across a rule file and `Snakefile`.
  - Avoid: `FOO_PLOT = "plots/geomx/foo.pdf"` when only one rule uses it.
  - Prefer: `output: plot="plots/geomx/foo.pdf"` directly in the owning rule.
- For reusable rules/scripts, keep interfaces minimal and explicit:
  - Inputs: only the files required by that script.
  - Params: only behavioral controls (for example `cluster_column_prefix`, `resolution`, `markers_group`).
  - Outputs: stable generic names that map cleanly to file paths.
- Keep `config.yaml` focused on genuinely variable settings:
  - Prefer hard-coded values in rules/scripts for stable analysis constants that are not intended to vary between runs.
  - Only add config entries for user-facing settings, values reused across workflows, or settings that are genuinely expected to vary across runs/projects.
- For reproducible Snakemake scheduling in this clean tree, prefer the repo-local workflow profile defined by `profiles/config.yaml`.

Python Snakemake skeleton:
```python
if "snakemake" in locals():
    in_path = snakemake.input["in_path"]
    out_path = snakemake.output["out_path"]
else:
    in_path = "path/to/local/input"
    out_path = "path/to/local/output"
```

R Snakemake skeleton:
```r
if (exists("snakemake")) {
  in_path <- snakemake@input[["in_path"]]
  out_path <- snakemake@output[["out_path"]]
} else {
  in_path <- "path/to/local/input"
  out_path <- "path/to/local/output"
}
```

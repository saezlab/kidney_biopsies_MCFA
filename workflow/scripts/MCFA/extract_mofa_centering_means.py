# 07.07.2026 Charlotte Boys

# %% Import libraries
import h5py
import pandas as pd


def _decode_array(values):
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


# %% Parse snakemake object
if "snakemake" in locals():
    model = snakemake.input["model"]
    csv = snakemake.output["csv"]
else:
    model = "results/MCFA/model~hvg~5.h5ad"
    csv = "results/MCFA/centering_means~hvg~5.csv"


# %% Extract MOFA centering means
rows = []

with h5py.File(model, "r") as handle:
    views = _decode_array(handle["views/views"][:])

    for view in views:
        features = _decode_array(handle[f"features/{view}"][:])
        intercept_groups = list(handle[f"intercepts/{view}"].keys())

        for group in intercept_groups:
            centering_means = handle[f"intercepts/{view}/{group}"][:]

            if len(features) != len(centering_means):
                raise ValueError(
                    f"Feature and centering mean lengths differ for view {view}, "
                    f"group {group}: {len(features)} != {len(centering_means)}"
                )

            rows.append(
                pd.DataFrame(
                    {
                        "view": view,
                        "variable": features,
                        "group": group,
                        "centering_mean": centering_means,
                    }
                )
            )

centering_means = pd.concat(rows, ignore_index=True)


# %% Save centering means
centering_means.to_csv(csv, index=False)

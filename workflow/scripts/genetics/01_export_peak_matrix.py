# (C) Sergio Alejandro Gómez
# Aim: Sum the ArchR PeakMatrix (peak_mat2.h5ad, nuclei x peaks, Tn5 insertion counts capped at 4) into a cluster x library
#   pseudobulk with integer arithmetic, and write it in a form base R reads without an HDF5 reader. All nuclei are kept;
#   the cohort and the pt_C5 filters are applied later in R.
# Inputs: the h5ad file (obs columns sample, diseaseNames, cluster; X = CSR sparse matrix).
# Outputs (in --out_dir): pseudobulk_counts_int32_le.bin (int32 little-endian, column-major peaks x groups),
#   pseudobulk_peaks.tsv, pseudobulk_groups.tsv (group = cluster|library), cells.tsv.gz, export_summary.json.
# Usage: python workflow/scripts/genetics/01_export_peak_matrix.py --h5ad=data/atac/peak_mat2.h5ad --out_dir=results/genetics/full_cohort/peak_matrix_export
import argparse, gzip, io, json, os, sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
import h5py
import numpy as np


def decode(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a], dtype=object)


def obs_column(f, col):
    o = f["obs"][col]
    if isinstance(o, h5py.Group):  # anndata categorical
        codes = o["codes"][:]
        if np.any(codes < 0):
            raise SystemExit(f"obs.{col} has missing values")
        return decode(o["categories"][:])[codes]
    return decode(o[:])


def attr_str(obj, key, default=None):
    v = obj.attrs.get(key, default)
    return v.decode() if isinstance(v, bytes) else v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--sample_col", default="sample")
    ap.add_argument("--disease_col", default="diseaseNames")
    ap.add_argument("--cluster_col", default="cluster")
    ap.add_argument("--rows_per_chunk", type=int, default=4000)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    with h5py.File(a.h5ad, "r") as f:
        peaks = decode(f["var"][attr_str(f["var"], "_index", "_index")][:])
        barcodes = decode(f["obs"][attr_str(f["obs"], "_index", "_index")][:])
        sample = obs_column(f, a.sample_col)
        disease = obs_column(f, a.disease_col)
        cluster = obs_column(f, a.cluster_col)
        X = f["X"]
        if not isinstance(X, h5py.Group) or not (attr_str(X, "encoding-type") == "csr_matrix" or attr_str(X, "h5sparse_format") == "csr"):
            raise SystemExit("X is not a CSR sparse matrix (nuclei x peaks)")
        indptr = X["indptr"][:].astype(np.int64)
        n_cells, n_peaks = len(indptr) - 1, len(peaks)
        if n_cells != len(barcodes):
            raise SystemExit("obs and X disagree on the number of nuclei")

        # one pseudobulk group per cluster|library, sorted
        key = np.array([f"{c}|{s}" for c, s in zip(cluster, sample)], dtype=object)
        groups, gidx = np.unique(key, return_inverse=True)
        acc = np.zeros((len(groups), n_peaks), dtype=np.int64)  # row g = column g of the peaks x groups matrix
        spectrum, total, row_sum = {}, 0, np.zeros(n_cells, dtype=np.int64)
        for r0 in range(0, n_cells, a.rows_per_chunk):
            r1 = min(r0 + a.rows_per_chunk, n_cells)
            lo, hi = indptr[r0], indptr[r1]
            d = X["data"][lo:hi]
            j = X["indices"][lo:hi].astype(np.int64)
            if d.size:
                if np.any(d != np.round(d)) or d.min() < 1 or d.max() > 4:
                    raise SystemExit("matrix values are not integers in 1-4")
                d = d.astype(np.int64)
                rows = np.repeat(np.arange(r0, r1), np.diff(indptr[r0:r1 + 1]))
                np.add.at(row_sum, rows, d)
                np.add.at(acc, (gidx[rows], j), d)
                v, c = np.unique(d, return_counts=True)
                for vi, ci in zip(v, c):
                    spectrum[int(vi)] = spectrum.get(int(vi), 0) + int(ci)
                total += int(d.sum())
            if (r1 // a.rows_per_chunk) % 4 == 0 or r1 == n_cells:
                print(f"rows {r1}/{n_cells}", flush=True)
        nnz = int(indptr[-1])

    if int(acc.sum()) != total or acc.max() >= 2 ** 31:
        raise SystemExit("pseudobulk total differs from the matrix total, or a count exceeds int32")
    acc.astype("<i4").tofile(os.path.join(a.out_dir, "pseudobulk_counts_int32_le.bin"))
    with open(os.path.join(a.out_dir, "pseudobulk_peaks.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("peak\n")
        fh.writelines(f"{p}\n" for p in peaks)
    disease_of = {}
    for s, dz in zip(sample, disease):
        if disease_of.setdefault(s, dz) != dz:
            raise SystemExit(f"library {s} has more than one disease group")
    n_nuclei = np.bincount(gidx, minlength=len(groups))
    with open(os.path.join(a.out_dir, "pseudobulk_groups.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("group\tcluster\tsample\tdiseaseNames\tn_nuclei\n")
        for g, n in zip(groups, n_nuclei):
            c, s = g.split("|", 1)
            fh.write(f"{g}\t{c}\t{s}\t{disease_of[s]}\t{int(n)}\n")
    with gzip.open(os.path.join(a.out_dir, "cells.tsv.gz"), "wt", encoding="utf-8", newline="\n") as fh:
        fh.write("barcode\tsample\tdiseaseNames\tcluster\tinsertion_row_sum\n")
        for b, s, dz, c, r in zip(barcodes, sample, disease, cluster, row_sum):
            fh.write(f"{b}\t{s}\t{dz}\t{c}\t{int(r)}\n")
    summary = {"input": os.path.abspath(a.h5ad), "n_nuclei": int(n_cells), "n_peaks": int(n_peaks), "n_groups": int(len(groups)),
               "nnz": nnz, "total_count": int(total), "value_spectrum": {str(k): v for k, v in sorted(spectrum.items())},
               "python": sys.version.split()[0], "h5py": h5py.__version__, "numpy": np.__version__}
    with open(os.path.join(a.out_dir, "export_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k != "value_spectrum"}, indent=1))


if __name__ == "__main__":
    main()

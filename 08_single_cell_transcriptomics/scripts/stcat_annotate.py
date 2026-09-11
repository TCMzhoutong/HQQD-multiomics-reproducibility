import argparse
import os
import warnings

import pandas as pd
import scipy.io
import scanpy as sc
import STCAT


warnings.filterwarnings("ignore")


def load_from_export_dir(export_dir: str) -> sc.AnnData:
    """Assemble AnnData from R sparse matrix export directory."""
    # Stored as genes x cells in mtx; transpose to cells x genes
    counts = scipy.io.mmread(os.path.join(export_dir, "counts.mtx")).T.tocsr()
    data   = scipy.io.mmread(os.path.join(export_dir, "data.mtx")).T.tocsr()

    barcodes = pd.read_csv(
        os.path.join(export_dir, "barcodes.tsv"), header=None
    )[0].tolist()
    features = pd.read_csv(
        os.path.join(export_dir, "features.tsv"), header=None
    )[0].tolist()

    obs = pd.read_csv(os.path.join(export_dir, "obs_metadata.csv"), index_col=0)

    # STCAT requires adata.X to be raw counts (integer matrix)
    adata = sc.AnnData(X=counts)
    adata.obs_names = barcodes
    adata.var_names = features
    adata.layers["counts"] = counts
    adata.layers["data"] = data
    adata.obs = obs.reindex(barcodes)

    umap_path = os.path.join(export_dir, "umap.csv")
    if os.path.exists(umap_path):
        umap_df = pd.read_csv(umap_path, index_col=0)
        umap_cols = [c for c in umap_df.columns if c.lower() != "barcode"]
        adata.obsm["X_umap"] = umap_df.reindex(barcodes)[umap_cols].values

    return adata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run STCAT on T-cell data.")
    parser.add_argument(
        "--input", required=True,
        help="R export directory (counts.mtx / data.mtx / barcodes.tsv ...) or h5ad file"
    )
    parser.add_argument("--output", required=True, help="Output CSV path for predictions")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if os.path.isdir(args.input):
        adata = load_from_export_dir(args.input)
    else:
        adata = sc.read_h5ad(args.input)

    adata = STCAT.STCAT(adata)

    cols = ["Prediction", "Uncertainty score", "Cluster"]
    keep_cols = [c for c in cols if c in adata.obs.columns]
    if not keep_cols:
        raise RuntimeError("STCAT output is missing expected columns in adata.obs")

    out_df = adata.obs[keep_cols].copy()
    out_df.index.name = "Barcode"
    out_df.to_csv(args.output)
    print(f"STCAT finished. Saved predictions to {args.output}")


if __name__ == "__main__":
    import traceback
    try:
        main()
    except Exception:
        traceback.print_exc()
        raise SystemExit(1)

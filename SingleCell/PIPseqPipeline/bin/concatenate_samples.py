#!/usr/bin/env python3
"""Concatenate per-subsample DRAGEN scRNA matrices into supersample-level AnnData files.

Memory strategy
---------------
Each subsample is read into memory once, immediately written to an uncompressed temporary
.h5ad (plus a CRISPR-features-only temporary .h5ad) and then freed. The temporary files are
then concatenated with ``anndata.experimental.concat_on_disk``, which streams row blocks from
the inputs into the output instead of holding every subsample (and the concatenated copy)
in memory at once. Peak memory therefore scales with the largest single subsample rather
than with the sum of all subsamples, which is what an in-memory ``anndata.concat`` needs.

The resulting files are identical (obs, var, X) to what ``anndata.concat(..., join='outer',
merge='same', label='subsample_id', keys=..., index_unique='_')`` followed by ``write_h5ad``
would produce; this equivalence was verified on the PotC_Pilot_2 CRISPR sublibraries.
"""

import argparse
import shutil
import sys
import tempfile
from pathlib import Path

import anndata as ad
import h5py
import scanpy as sc
import scipy.sparse as sp
from anndata.experimental import concat_on_disk
from anndata.io import read_elem, write_elem

sc.settings.n_jobs = -1

CRISPR_FEATURE_TYPE = 'CRISPR Direct Capture'
SUBSAMPLE_LABEL = 'subsample_id'
INDEX_UNIQUE = '_'
# Number of array elements copied per slice when recompressing X (~80-160 MB per slice).
RECOMPRESS_CHUNK_ELEMS = 20_000_000


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Concatenate multiple subsamples using AnnData"
    )
    parser.add_argument(
        "--matrices",
        type=str,
        nargs='+',
        required=True,
        help="Paths to filtered matrix files"
    )
    parser.add_argument(
        "--barcodes",
        type=str,
        nargs='+',
        required=True,
        help="Paths to filtered barcodes files"
    )
    parser.add_argument(
        "--features",
        type=str,
        nargs='+',
        required=True,
        help="Paths to filtered features files"
    )
    parser.add_argument(
        "--subsample-ids",
        type=str,
        required=True,
        help="Comma-separated list of subsample IDs"
    )
    parser.add_argument(
        "--supersample-basename",
        type=str,
        required=True,
        help="Basename for output files"
    )
    return parser.parse_args()


def read_10x(matrix, barcodes, features):
    """Read one subsample's 10x-format triplet (all feature types) as a CSR float32 AnnData."""
    # scanpy expects the three files in one directory with these exact names
    with tempfile.TemporaryDirectory() as tmpdir:
        shutil.copy(matrix, f"{tmpdir}/matrix.mtx.gz")
        shutil.copy(barcodes, f"{tmpdir}/barcodes.tsv.gz")
        shutil.copy(features, f"{tmpdir}/features.tsv.gz")
        adata = sc.read_10x_mtx(tmpdir, var_names='gene_symbols', gex_only=False, cache=False)
    # read_10x_mtx returns CSC (it transposes the genes-x-cells .mtx); concat_on_disk can only
    # concatenate CSR along the cell axis, and CSR is the natural layout for cells-x-features
    # anyway, so the published .h5ad files are now CSR (they used to be CSC).
    adata.X = sp.csr_matrix(adata.X)
    return adata


def concat_h5ad_files_on_disk(in_files, keys, out_file):
    """On-disk equivalent of
    ad.concat(adatas, join='outer', merge='same', label=SUBSAMPLE_LABEL, keys=keys, index_unique=INDEX_UNIQUE).

    A single input is copied through unchanged (no label column, no index suffix), matching
    what the in-memory code path used to do when it skipped concatenation for one subsample.

    anndata 0.11.x has a bug in concat_on_disk(merge='same'): it re-reads `var` as a mapping and
    overwrites `/var` in the output with a broken object. So we concatenate with merge=None (which
    writes a `var` with the right index but no columns) and then compute the merged `var` columns
    ourselves, with the exact same semantics, by running the in-memory ad.concat on var-only
    AnnData objects (tiny: no X) and writing the result into the output file.
    """
    if len(in_files) == 1:
        shutil.copy(in_files[0], out_file)
        return

    concat_on_disk(
        in_files, out_file,
        join='outer', merge=None, label=SUBSAMPLE_LABEL, keys=keys, index_unique=INDEX_UNIQUE,
    )

    var_only = []
    for f in in_files:
        with h5py.File(f, 'r') as h:
            var_only.append(ad.AnnData(var=read_elem(h['var'])))
    merged_var = ad.concat(var_only, join='outer', merge='same').var

    with h5py.File(out_file, 'a') as h:
        out_var_names = read_elem(h['var']).index
        del h['var']
        write_elem(h, 'var', merged_var.loc[out_var_names])


def write_compressed_copy(src, dst, chunk=RECOMPRESS_CHUNK_ELEMS):
    """Copy h5ad `src` to `dst`, gzip-compressing the sparse X arrays slice by slice.

    concat_on_disk has no compression option, so its output is written uncompressed; this pass
    produces the same gzip-compressed file layout that write_h5ad(compression='gzip') used to,
    without ever loading X into memory (the CSR data/indices/indptr arrays are copied in bounded
    slices). Everything other than X (obs, var, ...) is small, so it is simply read and re-written
    through anndata's own element writer with compression.
    """
    dataset_kwargs = dict(compression='gzip')
    with h5py.File(src, 'r') as fi, h5py.File(dst, 'w') as fo:
        fo.attrs.update(fi.attrs)
        for key in fi:
            if key != 'X':
                write_elem(fo, key, read_elem(fi[key]), dataset_kwargs=dataset_kwargs)
        x_in = fi['X']
        x_out = fo.create_group('X')
        x_out.attrs.update(x_in.attrs)
        for name in x_in:
            d_in = x_in[name]
            d_out = x_out.create_dataset(
                name, shape=d_in.shape, dtype=d_in.dtype, maxshape=(None,), **dataset_kwargs
            )
            d_out.attrs.update(d_in.attrs)
            for start in range(0, d_in.shape[0], chunk):
                d_out[start:start + chunk] = d_in[start:start + chunk]


def main():
    """Main execution function."""
    args = parse_args()
    try:
        print(f"Concatenating {len(args.matrices)} subsamples...")

        subsample_ids = args.subsample_ids.split(',')

        if len(args.matrices) != len(args.barcodes) or len(args.matrices) != len(args.features):
            raise ValueError("Number of matrix, barcode, and feature files must match")

        if len(args.matrices) != len(subsample_ids):
            raise ValueError("Number of files must match number of subsample IDs")

        if len(subsample_ids) == 0:
            raise ValueError("No subsamples found")

        output_path = f"{args.supersample_basename}.h5ad"
        crispr_output_path = f"{args.supersample_basename}.crispr.h5ad"

        # Temporary files go into the task's working directory (not /tmp), which is where the
        # scratch space allocated to this process lives.
        with tempfile.TemporaryDirectory(dir='.', prefix='concat_tmp_') as tmpdir:
            tmpdir = Path(tmpdir)
            full_files, crispr_files = [], []

            # Read each subsample once, spill it to disk, free it
            for subsample_id, matrix, barcode, feature in zip(
                subsample_ids, args.matrices, args.barcodes, args.features
            ):
                print(f"Reading subsample {subsample_id} from {matrix}...")
                adata = read_10x(matrix, barcode, feature)
                print(f"subsample {subsample_id}: {adata.n_obs} cells, {adata.n_vars} features")

                full_file = tmpdir / f"{subsample_id}.h5ad"
                crispr_file = tmpdir / f"{subsample_id}.crispr.h5ad"
                adata.write_h5ad(full_file)
                crispr_mask = (adata.var['feature_types'] == CRISPR_FEATURE_TYPE).values
                adata[:, crispr_mask].copy().write_h5ad(crispr_file)
                full_files.append(full_file)
                crispr_files.append(crispr_file)
                del adata

            print(f"\nConcatenating {len(full_files)} subsamples on disk...")
            concat_file = tmpdir / "concatenated.h5ad"
            concat_crispr_file = tmpdir / "concatenated.crispr.h5ad"
            concat_h5ad_files_on_disk(full_files, subsample_ids, concat_file)
            concat_h5ad_files_on_disk(crispr_files, subsample_ids, concat_crispr_file)

            with h5py.File(concat_file, 'r') as h:
                n_obs, n_vars = h['X'].attrs['shape']
            print(f"Concatenated dataset: {n_obs} cells, {n_vars} features")

            print(f"Writing data to {output_path}...")
            write_compressed_copy(concat_file, output_path)
            print(f"Writing CRISPR data to {crispr_output_path}...")
            write_compressed_copy(concat_crispr_file, crispr_output_path)

        print("\n✓ Concatenation completed successfully")
        return 0

    except Exception as e:
        print(f"\n✗ Error during concatenation: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())

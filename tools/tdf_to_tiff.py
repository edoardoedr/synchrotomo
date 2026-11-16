"""Estrae le immagini di un file TDF (dataset `exchange/data`) in TIFF singoli.

Ogni immagine e' data[:, i, :] (indice sulla seconda dimensione), salvata come
<prefix>_0000.tif, <prefix>_0001.tif, ...

Uso:
    python tools/tdf_to_tiff.py /path/to/file.tdf /path/to/output_dir [--prefix img]
"""
from __future__ import print_function

import argparse
import os

import h5py
import tifffile

DATASET_PATH = "exchange/data"


def extract_stack_to_tif(tdf_path, output_dir, prefix="tomo"):
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    with h5py.File(tdf_path, "r") as f:
        if DATASET_PATH not in f:
            f.visit(print)
            raise KeyError("Dataset '{}' non trovato in {}".format(DATASET_PATH, tdf_path))

        dset = f[DATASET_PATH]
        if dset.ndim != 3:
            raise ValueError("Atteso stack 3D, trovato shape {}".format(dset.shape))

        n_images = dset.shape[1]
        print("Trovate {} immagini".format(n_images))
        for i in range(n_images):
            out = os.path.join(output_dir, "{}_{:04d}.tif".format(prefix, i))
            tifffile.imwrite(out, dset[:, i, :]) if hasattr(tifffile, "imwrite") \
                else tifffile.imsave(out, dset[:, i, :])

    print("Conversione completata.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("tdf_path")
    parser.add_argument("output_dir")
    parser.add_argument("--prefix", default="tomo")
    args = parser.parse_args()

    if not os.path.exists(args.tdf_path):
        raise IOError("File non trovato: {}".format(args.tdf_path))
    extract_stack_to_tif(args.tdf_path, args.output_dir, args.prefix)

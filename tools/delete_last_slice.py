"""Rimuove l'ultima slice lungo l'asse y (seconda dimensione) di un file TDF, in place.

Utile per allineare le dimensioni di volumi acquisiti a dosi diverse.

Uso:
    python tools/delete_last_slice.py /path/to/file.tdf
"""
import argparse

import h5py

DATASET_PATH = "exchange/data"


def delete_last_slice(file_path):
    with h5py.File(file_path, "r+") as f:
        data = f[DATASET_PATH][:]
        del f[DATASET_PATH]
        f[DATASET_PATH] = data[:, :-1, :]
    print("Shape {} -> {}".format(data.shape, (data.shape[0], data.shape[1] - 1, data.shape[2])))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("tdf_path")
    delete_last_slice(parser.parse_args().tdf_path)

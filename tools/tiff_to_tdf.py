# -*- coding: utf-8 -*-
"""Esegue solo la conversione TIFF -> TDF (sezione `tiff2tdf` del file YAML).

Uso (dalla radice della repo):
    python tools/tiff_to_tdf.py config_reconstruction.yaml
"""
from __future__ import print_function

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main_reconstruction import load_config, run_tiff2tdf_from_config  # noqa: E402

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python tools/tiff_to_tdf.py <config_reconstruction.yaml>")
        sys.exit(1)
    config = load_config(sys.argv[1])
    config.setdefault('tiff2tdf', {})['enabled'] = True
    print(run_tiff2tdf_from_config(config))

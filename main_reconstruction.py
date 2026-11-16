# -*- coding: utf-8 -*-
"""
Pipeline di ricostruzione tomografica guidata da un file YAML.

Step (ognuno abilitabile/disabilitabile dal file di configurazione):
    1. tiff2tdf          - conversione di proiezioni/sinogrammi TIFF in un file TDF (HDF5)
    2. phase_retrieval   - phase retrieval (TIE-Hom / CTF) sul file TDF
    3. reconstruct       - ricostruzione delle slice (FBP, iterativi, GRIDREC, MR-FBP, FISTA-TV)

Uso:
    python main_reconstruction.py config_reconstruction.yaml
"""
from __future__ import print_function

import os
import sys

import yaml

from src_reconstruction import exec_tiff2tdf
from src_reconstruction import exec_phaseretrieval
from src_reconstruction import exec_reconstruct


def load_config(config_file):
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)


def run_tiff2tdf_from_config(config):
    """Step 1: TIFF -> TDF. Restituisce il percorso del file TDF da usare negli step successivi.

    Se `tiff2tdf.enabled` e' False la conversione viene saltata e `tiff2tdf.output_file`
    viene usato come file TDF di input gia' esistente.
    """
    cfg = config.get('tiff2tdf', {})
    output_file = cfg.get('output_file', '')

    if not cfg.get('enabled', True):
        print("[1/3] TIFF -> TDF disabilitato. Uso il TDF esistente: {}".format(output_file))
        return output_file

    argv = [
        str(cfg.get('from', 0)),
        str(cfg.get('to', -1)),
        cfg.get('input_path', ''),
        output_file,
        str(cfg.get('crop_top', 0)),
        str(cfg.get('crop_bottom', 0)),
        str(cfg.get('crop_left', 0)),
        str(cfg.get('crop_right', 0)),
        cfg.get('tomo_prefix', 'tomo'),
        cfg.get('flat_prefix', 'flat'),
        cfg.get('dark_prefix', 'dark'),
        str(cfg.get('projection_order', True)),
        str(cfg.get('privilege_sino', True)),
        str(cfg.get('compression', 5)),
        str(cfg.get('nr_threads', 4)),
        cfg.get('log_file', 'conversion.log'),
    ]

    print("[1/3] Conversione TIFF -> TDF: {} -> {}".format(argv[2], output_file))
    exec_tiff2tdf.main(argv)

    if not os.path.exists(output_file):
        raise IOError("Conversione TIFF -> TDF fallita: {} non trovato".format(output_file))
    print("[1/3] Conversione TIFF -> TDF completata.")
    return output_file


def run_phaseretrieval_from_config(config, input_tdf_file):
    """Step 2: phase retrieval. Restituisce il TDF da passare alla ricostruzione."""
    cfg = config.get('phase_retrieval', {})

    if not cfg.get('enabled', False):
        print("[2/3] Phase retrieval disabilitato. Skip.")
        return input_tdf_file

    output_file = cfg.get('output_file', input_tdf_file.replace('.tdf', '_phrt.tdf'))

    argv = [
        str(cfg.get('from', config.get('from', 0))),
        str(cfg.get('to', config.get('to', -1))),
        input_tdf_file,
        output_file,
        str(cfg.get('method', 0)),
        str(cfg.get('param1', 1.0)),
        str(cfg.get('param2', 2000.0)),
        str(cfg.get('energy', 22.0)),
        str(cfg.get('distance', 300.0)),
        str(cfg.get('pixsize', 2.2)),
        str(cfg.get('pad', True)),
        str(cfg.get('nr_threads', config.get('nr_threads', 4))),
        cfg.get('log_file', 'phaseretrieval.log'),
    ]

    print("[2/3] Phase retrieval: {} -> {}".format(input_tdf_file, output_file))
    exec_phaseretrieval.main(argv)
    print("[2/3] Phase retrieval completato.")
    return output_file


def run_reconstruct_from_config(config, input_tdf_file):
    """Step 3: ricostruzione tomografica delle slice dal file TDF."""
    cfg = config.get('reconstruct', {})

    if not cfg.get('enabled', False):
        print("[3/3] Ricostruzione disabilitata. Skip.")
        return

    nr_threads = cfg.get('nr_threads', config.get('tiff2tdf', {}).get('nr_threads', 4))

    # L'ordine degli argomenti e' quello atteso da exec_reconstruct.main (42 parametri).
    argv = [
        str(cfg.get('from', 0)),
        str(cfg.get('to', -1)),
        input_tdf_file,
        cfg.get('output_path', './reconstructed/'),
        str(cfg.get('angles', 3.1416)),
        str(cfg.get('offset', 0.0)),
        str(cfg.get('param1', 'shepp-logan')),
        str(cfg.get('scale', 1.0)),
        str(cfg.get('overpad', True)),
        str(cfg.get('logtransform', True)),
        str(cfg.get('circle', True)),
        cfg.get('output_prefix', 'recon'),
        str(cfg.get('preprocessing_required', False)),
        str(cfg.get('flat_end', False)),
        str(cfg.get('half_half', False)),
        str(cfg.get('half_half_line', 0)),
        str(cfg.get('ext_fov', False)),
        str(cfg.get('ext_fov_rot_right', False)),
        str(cfg.get('ext_fov_overlap', 0)),
        str(cfg.get('ext_fov_normalize', False)),
        str(cfg.get('ext_fov_average', False)),
        str(cfg.get('ext_fov_proj', 0)),
        str(cfg.get('norm_sx', 0)),
        str(cfg.get('norm_dx', 0)),
        str(cfg.get('skip_ringrem', True)),
        cfg.get('ringrem', 'none:'),
        str(cfg.get('zerone_mode', False)),
        str(cfg.get('corr_offset', 0.0)),
        cfg.get('method', 'FBP_CUDA'),
        str(cfg.get('decim_factor', 1)),
        str(cfg.get('downsc_factor', 1)),
        str(cfg.get('postprocess_required', False)),
        cfg.get('polarfilt_opt', 'none:'),
        cfg.get('convert_opt', 'none:'),
        cfg.get('crop_opt', '0:0:0:0'),
        str(cfg.get('angles_projfrom', 0)),
        str(cfg.get('angles_projto', 1799)),
        str(cfg.get('rolling', False)),
        str(cfg.get('roll_shift', 0)),
        str(cfg.get('dynamic_ff', False)),
        str(nr_threads),
        cfg.get('log_file', 'reconstruction.log'),
    ]

    print("[3/3] Ricostruzione ({}): {} -> {}".format(argv[28], input_tdf_file, argv[3]))
    exec_reconstruct.main(argv)
    print("[3/3] Ricostruzione completata.")


def main(config_file):
    config = load_config(config_file)
    tdf_file = run_tiff2tdf_from_config(config)
    phrt_file = run_phaseretrieval_from_config(config, tdf_file)
    run_reconstruct_from_config(config, phrt_file or tdf_file)
    print("\n=== Pipeline completata ===")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python main_reconstruction.py <config_reconstruction.yaml>")
        sys.exit(1)
    main(sys.argv[1])

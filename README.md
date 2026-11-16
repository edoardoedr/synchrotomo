# synchrotomo

A Python pipeline for the **tomographic reconstruction of synchrotron CT data** (e.g. phase-contrast breast CT).

This project is a **Python adaptation of the original STP library** (STP-Core, part of the SYRMEP Tomo Project by Elettra – Sincrotrone Trieste). The original code has been extracted from the STP toolbox, adapted to synchrotron acquisitions, and wrapped in a single pipeline driven by a YAML configuration file.

Starting from projections in TIFF format, the pipeline produces the reconstructed slices, going through the TDF (HDF5) format and an optional phase-retrieval step.

```
 TIFF (projections)  ──►  [1] tiff2tdf  ──►  .tdf  ──►  [2] phase retrieval  ──►  _phrt.tdf  ──►  [3] reconstruct  ──►  TIFF slices
```

Each step can be enabled or disabled independently (`enabled: True/False`).

> ⚠️ **Work in progress.** This project is still under active development. The current code is a working adaptation of STP-Core, which still relies on Python 2.7 and on legacy dependencies. The long-term goal is to **fully re-implement the pipeline from scratch**, possibly in **Rust**, to get a modern, faster and easier-to-maintain codebase.

---

## Repository layout

```
.
├── main_reconstruction.py       # entry point: runs the pipeline from the YAML file
├── config_reconstruction.yaml   # example configuration, with comments
├── src_reconstruction/          # reconstruction core (Python adaptation of STP-Core)
│   ├── exec_tiff2tdf.py         #   step 1: TIFF -> TDF
│   ├── exec_phaseretrieval.py   #   step 2: phase retrieval
│   ├── exec_reconstruct.py      #   step 3: reconstruction
│   ├── reconstruct/             #   algorithms (ASTRA, GRIDREC, MR-FBP, FISTA-TV, SIRT-FBP)
│   ├── phaseretrieval/          #   TIE-Hom 2002/2020, CTF, BAC
│   ├── preprocess/              #   flat fielding (including dynamic), ring removal, extended FOV
│   ├── postprocess/             #   polar filters, conversion/crop
│   ├── stpio/                   #   TDF and EDF I/O
│   └── utils/                   #   padding, caching, center finding
├── packages_reconstruction/     # compiled extensions: tomopy (gridrec), pymrfbp, pytvtomo
├── tools/                       # standalone utilities for TDF files
├── environment.yml              # conda environment (Python 2.7)
├── install_dependencies.sh      # full installation (conda + ASTRA + local packages)
└── vendor/STP_1.6.3/            # original STP toolbox, kept for reference only (not used by the code)
```

---

## Installation

The code runs on **Python 2.7** and requires an **NVIDIA GPU with CUDA** for the `*_CUDA` methods (ASTRA Toolbox). Tested on Linux.

```bash
./install_dependencies.sh
conda activate synchrotomo
```

The script:
1. creates the `synchrotomo` conda environment (Python 2.7.15);
2. installs ASTRA Toolbox 2.1.2 from the `astra-toolbox` channel;
3. installs `numpy`, `scipy 0.12`, `h5py`, `pyfftw`, `opencv-python`, `pywavelets`, `tifffile`, `pyyaml`, ... via pip;
4. builds and installs the local packages in `packages_reconstruction/` (`pymrfbp`, `pytvtomo`, `gridrec`).

Alternatively, `environment.yml` pins the versions of the environment in use: `conda env create -f environment.yml` (ASTRA and the local packages must then be installed separately, as in the script).

Check the installation:

```bash
python -c "import astra, mrfbp, tvtomo; from _gridrec import paralrecon; print('OK')"
```

---

## Usage

```bash
python main_reconstruction.py config_reconstruction.yaml
```

Copy `config_reconstruction.yaml`, replace the placeholder paths (`/path/to/...`) and adapt the parameters to your acquisition. Each step writes its own log file (`log_file`).

Typical scenarios:
- **Reconstruction from TIFF, with phase retrieval**: all three steps `enabled: True`.
- **Reconstruction from an existing TDF**: `tiff2tdf.enabled: False` and `tiff2tdf.output_file` = path of the input TDF.
- **Data with phase retrieval already applied**: `phase_retrieval.enabled: False`.

---

## Parameters

### 1. `tiff2tdf` — TIFF → TDF conversion

| Parameter | Description |
|---|---|
| `from`, `to` | first/last image index (`-1` = all) |
| `input_path` | folder containing the TIFF files |
| `output_file` | output TDF file |
| `crop_top/bottom/left/right` | pixels to remove at the borders |
| `tomo_prefix`, `flat_prefix`, `dark_prefix` | file prefixes; `"-"` = no flat/dark |
| `projection_order` | `True` = images are projections, `False` = sinograms |
| `privilege_sino` | HDF5 layout optimized for sinogram access |
| `compression` | GZIP level 0–9 |
| `nr_threads` | number of parallel processes |

### 2. `phase_retrieval` — phase retrieval (optional)

| Parameter | Description |
|---|---|
| `method` | `0` TIE-Hom (Paganin 2002), `1` generalized TIE-Hom (Paganin 2020), `≥2` Projected CTF (Moosmann 2011) |
| `param1`, `param2` | TIE-Hom: `beta`, `delta` (the δ/β ratio is what matters). CTF: threshold, regularization |
| `energy` | beam energy [keV] |
| `distance` | sample–detector distance [mm] |
| `pixsize` | pixel size [µm] |
| `pad` | padding before the FFT (recommended) |
| `output_file` | TDF after phase retrieval (input of the reconstruction) |

### 3. `reconstruct` — reconstruction

**Geometry and method**

| Parameter | Description |
|---|---|
| `from`, `to` | range of slices (sinograms) to reconstruct |
| `angles` | total rotation angle [rad] (π = 180°) |
| `offset` | center-of-rotation offset [pixels] |
| `angles_projfrom`, `angles_projto` | projections to use |
| `method` | see the methods table |
| `param1` | FBP: filter · iterative: number of iterations · FISTA-TV: `"lambda:fgpiter:tviter"` |
| `scale`, `overpad`, `circle` | scale factor, extra padding, circular mask |
| `logtransform` | apply `-log`; `False` if data are already in absorption (e.g. after phase retrieval) |
| `output_path`, `output_prefix` | where/how slices are saved (`recon_0000.tif`, ...) |

**Available methods**

| `method` | Algorithm | Backend |
|---|---|---|
| `FBP_CUDA` | Filtered Back-Projection (filters `ram-lak`, `shepp-logan`, `cosine`, `hamming`, `hann`, ...) | ASTRA |
| `SIRT_CUDA`, `SART_CUDA`, `CGLS_CUDA` | algebraic iterative methods | ASTRA |
| `MR-FBP_CUDA` | Minimum-Residual FBP | `pymrfbp` + ASTRA |
| `FISTA-TV_CUDA` | FISTA with Total Variation regularization | `pytvtomo` + ASTRA |
| `GRIDREC` | gridding (CPU) | `tomopy` gridrec |

**On-the-fly pre-processing** (`preprocessing_required: True`)

- flat/dark fielding, including flats at the end of the scan (`flat_end`) and half-half acquisitions (`half_half`, `half_half_line`);
- dynamic flat fielding (`dynamic_ff`);
- extended field of view (`ext_fov*`), border normalization (`norm_sx`, `norm_dx`);
- **ring removal** (`skip_ringrem: False`, `ringrem: "method:params"`): `rivers`, `raven`, `oimoen`, `muenchetal`, `boinhaibel`, `sijberspostnov` — e.g. `"rivers:11;0"`; `"none:"` = disabled.

**Advanced**: `zerone_mode`, `corr_offset`, `decim_factor` (projection decimation), `downsc_factor` (detector downscaling), `rolling`/`roll_shift` (360° acquisitions).

**Post-processing** (`postprocess_required: True`)

| Parameter | Format |
|---|---|
| `polarfilt_opt` | filter in polar coordinates, e.g. `"homomorphic:p1;p2"`; `"none:"` = disabled |
| `convert_opt` | `"linear8:min;max"` → uint8, `"linear:min;max"` → uint16; `"none:"` = float32 |
| `crop_opt` | `"top:bottom:left:right"` |

---

## Tools

Standalone utilities for TDF files (HDF5 dataset `exchange/data`, shape `(rows, n_images, columns)`):

```bash
python tools/tiff_to_tdf.py config_reconstruction.yaml            # runs only step 1 of the pipeline
python tools/tdf_to_tiff.py file.tdf out_dir/ --prefix img          # extracts data[:, i, :] into single TIFFs
python tools/delete_last_slice.py file.tdf                          # removes the last image (axis 1), in place
```

---

## TDF format

TDF is an HDF5 file with data stored in `exchange/data` (plus optional `exchange/data_white` / `exchange/data_dark` for flats and darks) and `min`/`max` attributes. It can be optimized for projection or sinogram access (`privilege_sino`). Read/write functions are in [`src_reconstruction/stpio/tdf.py`](src_reconstruction/stpio/tdf.py).

---

## Roadmap

- [x] Extract and adapt the STP-Core reconstruction code into a standalone Python pipeline
- [x] Single YAML-driven entry point (TIFF → TDF → phase retrieval → reconstruction)
- [ ] Full re-implementation of the pipeline, possibly in Rust
- [ ] Remove the dependency on Python 2.7 and on legacy packages
- [ ] Tests and example datasets

---

## Credits and license

The code in `src_reconstruction/` is adapted from **STP-Core** (SYRMEP Tomo Project), © Elettra – Sincrotrone Trieste S.C.p.A., author Francesco Brun, released under the **GNU GPL v3**. The packages in `packages_reconstruction/` (`pymrfbp`, `pytvtomo`, `gridrec` from TomoPy) keep their respective licenses. The original STP 1.6.3 toolbox is kept in `vendor/` for reference only.

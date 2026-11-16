###########################################################################
# (C) 2016 Elettra - Sincrotrone Trieste S.C.p.A.. All rights reserved.   #
#                                                                         #
#                                                                         #
# This file is part of STP-Core, the Python core of SYRMEP Tomo Project,  #
# a software tool for the reconstruction of experimental CT datasets.     #
#                                                                         #
# STP-Core is free software: you can redistribute it and/or modify it     #
# under the terms of the GNU General Public License as published by the   #
# Free Software Foundation, either version 3 of the License, or (at your  #
# option) any later version.                                              #
#                                                                         #
# STP-Core is distributed in the hope that it will be useful, but WITHOUT #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or   #
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License    #
# for more details.                                                       #
#                                                                         #
# You should have received a copy of the GNU General Public License       #
# along with STP-Core. If not, see <http://www.gnu.org/licenses/>.        #
#                                                                         #
###########################################################################

#
# Author: Francesco Brun
# Last modified: July, 8th 2016
#

import astra
import numpy as np

# MR-FBP weight function (Pelt, 2013)
def mr_weight_filter(n):
    """Return the Fourier-space MR-FBP weight filter."""
    # High-pass regularization weighting
    freq = np.fft.fftfreq(n).astype(np.float32)
    return np.abs(freq)

def recon_mr_fbp(im, angles, cor_shift):
    """
    Reconstruct a sinogram with a simplified Minimum Residual FBP (MR-FBP)
    implementation using ASTRA Toolbox.
    """

    # ---------------------------
    # 1) ASTRA geometries
    # ---------------------------
    det_count = im.shape[1]
    angle_vec = np.linspace(0, angles, im.shape[0], endpoint=False)

    vol_geom = astra.create_vol_geom(det_count, det_count)
    proj_geom = astra.create_proj_geom(
        'parallel',
        1.0,
        det_count,
        angle_vec
    )

    # ---------------------------
    # 2) MR-FBP filtering
    # ---------------------------
    # Apply MR-FBP frequency-domain filter to sinogram
    filt = mr_weight_filter(det_count)
    im_f = np.fft.fft(im, axis=1)
    im_f *= filt
    im_filtered = np.real(np.fft.ifft(im_f, axis=1))

    # ---------------------------
    # 3) Standard ASTRA FBP reconstruction
    # ---------------------------
    sinogram_id = astra.data2d.create('-sino', proj_geom, im_filtered)
    vol_id = astra.data2d.create('-vol', vol_geom)

    cfg = astra.astra_dict('FBP')
    cfg['ReconstructionDataId'] = vol_id
    cfg['ProjectionDataId'] = sinogram_id

    alg_id = astra.algorithm.create(cfg)
    astra.algorithm.run(alg_id)

    rec = astra.data2d.get(vol_id)

    # Cleanup
    astra.algorithm.delete(alg_id)
    astra.data2d.delete(sinogram_id)
    astra.data2d.delete(vol_id)

    return rec.astype(np.float32)
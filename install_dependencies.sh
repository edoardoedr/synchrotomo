#!/bin/bash

# Installa le dipendenze della pipeline di ricostruzione CT
# Crea un nuovo ambiente conda chiamato synchrotomo con Python 2.7.15

set -e  # Esci se qualche comando fallisce

echo "=== Installazione dipendenze synchrotomo ==="
echo "Data: $(date)"
echo

# Aggiorna conda
echo "Aggiornamento conda..."
conda update -n base -c defaults conda -y

# Nome dell'ambiente
ENV_NAME="synchrotomo"

# Controlla se l'ambiente esiste già
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "L'ambiente '${ENV_NAME}' esiste già."
    read -p "Vuoi rimuoverlo e ricrearlo? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Rimozione ambiente esistente..."
        conda env remove -n $ENV_NAME -y
    else
        echo "Operazione annullata."
        exit 1
    fi
fi

# Crea il nuovo ambiente con Python 2.7.15
echo "Creazione ambiente conda '$ENV_NAME' con Python 2.7.15..."
conda create -n $ENV_NAME python=2.7.15 -y

# Attiva l'ambiente
echo "Attivazione ambiente '$ENV_NAME'..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

echo "Ambiente conda attivo: $CONDA_DEFAULT_ENV"

# Aggiungi i canali necessari nell'ordine corretto di priorità
echo "Configurazione canali conda..."
conda config --add channels conda-forge
conda config --add channels astra-toolbox
conda config --prepend channels conda-forge

# Installa le dipendenze conda dai canali ufficiali
echo "Installazione dipendenze conda..."

# Pacchetti dal canale principale (defaults)
conda install -c defaults -y \
    _libgcc_mutex=0.1 \
    _openmp_mutex=5.1 \
    blas=1.0 \
    ca-certificates=2025.11.4 \
    certifi=2020.6.20 \
    intel-openmp=2025.0.0 \
    libedit=3.1.20230828 \
    libffi=3.2.1 \
    libgcc=15.2.0 \
    libgcc-ng=15.2.0 \
    libgfortran-ng=7.5.0 \
    libgfortran4=7.5.0 \
    libgomp=15.2.0 \
    libstdcxx=15.2.0 \
    libstdcxx-ng=15.2.0 \
    libxcb=1.17.0 \
    libzlib=1.3.1 \
    mkl=2020.2 \
    ncurses=6.5 \
    openssl=1.1.1w \
    pthread-stubs=0.3 \
    python=2.7.15 \
    readline=7.0 \
    setuptools=44.0.0 \
    sqlite=3.33.0 \
    tk=8.6.15 \
    wheel=0.37.1 \
    xorg-libx11=1.8.12 \
    xorg-libxau=1.0.12 \
    xorg-libxdmcp=1.1.5 \
    xorg-xorgproto=2024.1 \
    zlib=1.3.1


# Installa ASTRA Toolbox dal suo canale specifico
echo "Installazione ASTRA Toolbox..."
conda install -c astra-toolbox -y \
    astra-toolbox=2.1.2 \
    libastra=2.1.2

# Installa i pacchetti Python tramite pip
echo "Installazione pacchetti pip..."
pip install --upgrade pip==9.0.3

pip install \
    cython==0.29.36 \
    enum34==1.1.10 \
    futures==3.4.0 \
    h5py==2.10.0 \
    numpy==1.16.0 \
    opencv-python==3.4.2.17 \
    pathlib==1.0.1 \
    pillow==6.2.2 \
    pyfftw==0.12.0 \
    pywavelets==1.0.3 \
    pyyaml==5.4.1 \
    scipy==0.12.0 \
    six==1.17.0 \
    tifffile==2018.10.18

echo "Dipendenze conda e pip installate con successo!"
echo

# Installazione dei pacchetti locali
echo "=== Installazione pacchetti locali ==="

# Controlla se siamo nella directory corretta
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGES_DIR="$SCRIPT_DIR/packages_reconstruction"

if [ ! -d "$PACKAGES_DIR" ]; then
    echo "Errore: Directory packages non trovata in $PACKAGES_DIR"
    exit 1
fi

cd "$PACKAGES_DIR"

# Installa pymrfbp
echo "Installazione pymrfbp..."
if [ -d "pymrfbp" ]; then
    cd pymrfbp
    python setup.py build
    python setup.py install
    cd ..
    echo "pymrfbp installato con successo!"
else
    echo "Avviso: Directory pymrfbp non trovata"
fi

# Installa pytvtomo
echo "Installazione pytvtomo..."
if [ -d "pytvtomo" ]; then
    cd pytvtomo
    python setup.py build_ext --inplace
    python setup.py install
    cd ..
    echo "pytvtomo installato con successo!"
else
    echo "Avviso: Directory pytvtomo non trovata"
fi

# Installa tomopy/gridrec
echo "Installazione gridrec (tomopy)..."
if [ -d "tomopy" ]; then
    cd tomopy
    python setup.py build_ext --inplace
    python setup.py install
    cd ..
    echo "gridrec installato con successo!"
else
    echo "Avviso: Directory tomopy non trovata"
fi

cd "$SCRIPT_DIR"

echo
echo "=== Installazione completata ==="
echo "Tutte le dipendenze e i pacchetti locali sono stati installati."
echo "Ambiente conda '$ENV_NAME' creato e configurato con successo."
echo
echo "Per attivare l'ambiente in futuro, esegui:"
echo "conda activate $ENV_NAME"
echo
echo "Per verificare l'installazione, attiva l'ambiente e esegui:"
echo "python -c \"import astra; import mrfbp; import tvtomo; from _gridrec import paralrecon; print('Tutti i moduli importati con successo!')\""
echo
echo "Versioni installate:"
echo "Python: $(python --version 2>&1)"
echo "Pip: $(pip --version)"
echo
echo "Se riscontri problemi, controlla i log sopra per eventuali errori."
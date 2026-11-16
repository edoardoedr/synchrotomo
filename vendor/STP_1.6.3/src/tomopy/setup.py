# tifffile_setup.py
# Usage: ``python gridrec_setup.py build_ext --inplace``
from distutils.core import setup, Extension
import numpy, os
#os.environ["CC"]="icc -fast"
#os.environ["CCX"]="icc -fast"
os.environ["CC"]="gcc -O3"
os.environ["CCX"]="gcc -O3"
setup(
    name='gridrec',  # Rimosso l'underscore
    version='1.0.0',  # Aggiunta versione esplicita
    ext_modules=[Extension('_gridrec', ['gridrec.c'], include_dirs=[numpy.get_include()])]
)
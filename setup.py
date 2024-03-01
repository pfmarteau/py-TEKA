from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from glob import glob

# Use python setup.py build_ext --inplace
# to compile


setup(
  name = "teka", version='0.1.2', description='wrapper for TEKA',
  ext_modules = cythonize(["teka.pyx"])
)

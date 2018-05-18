'''To compile, run the following in the terminal:
    python setup.py build_ext --inplace
'''

from distutils.core import setup
from Cython.Build import cythonize

setup(
	name = 'CythonIzotov06Method',
	ext_modules = cythonize(["metal_Izotov06_Cython.pyx", "Izotov06_error_Cython.pyx"]),
)
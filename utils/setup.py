# Build with:
# python setup.py build_ext --inplace

from setuptools import setup
from Cython.Build import cythonize

setup(
        name='lowess',
            ext_modules=cythonize("_smoothers_lowess.pyx"),
                zip_safe=False,
    )

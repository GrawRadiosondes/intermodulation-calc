# Metadata lives in pyproject.toml, this file only declares the C extension.
# Build in place:  python setup.py build_ext --inplace   (or install: pip install .)
from setuptools import Extension, setup

setup(ext_modules=[Extension('imcore', ['imcore.c'])])

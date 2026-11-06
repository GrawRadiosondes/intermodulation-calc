# Build the compiled search core in place:  python setup.py build_ext --inplace
from setuptools import Extension, setup

setup(
    name='intermodulation-calc',
    ext_modules=[Extension('imcore', ['imcore.c'], extra_compile_args=[])],
    py_modules=[],
)

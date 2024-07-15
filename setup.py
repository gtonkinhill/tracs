from setuptools import setup, Extension, find_packages
import platform
from glob import glob
from pybind11.setup_helpers import Pybind11Extension, build_ext
from ctypes.util import find_library
from codecs import open
from os import path
import os
import re
import io


# Get version strip
def read(*names, **kwargs):
    with io.open(os.path.join(os.path.dirname(__file__), *names),
                 encoding=kwargs.get("encoding", "utf8")) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

# deal with openmp etc
openmp = os.getenv("ompy_OpenMP")
if openmp is None and platform.system() == 'Darwin':  # Check if macOS
    if find_library("omp") != None:
        openmp = True
        print("libOMP found, building with OpenMP")
    else:
        print("libOMP not found, building without OpenMP")
elif openmp in (None, True, "True", "true"):
    openmp = True
elif openmp in (False, "False", "false"):
    openmp = False
    print("Building without OpenMP")
else:
    raise ValueError("Env var ompy_OpenMP must be either True or False "
                     "(or not set); use eg. 'export ompy_OpenMP=False'"
                     f"Now it is: {openmp}")

extra_compile_args = ["-O3", "-ffast-math", "-march=native"]
extra_link_args = ["-lz"]
if openmp and platform.system() == 'Darwin':
    extra_compile_args.insert(-1, "-Xpreprocessor")
    extra_compile_args.insert(-1, "-fopenmp")
    extra_link_args.insert(-1, "-lomp")
elif openmp:
    extra_compile_args.insert(-1, "-fopenmp")
    extra_link_args.insert(-1, "-fopenmp")


ext_modules = [
    Pybind11Extension("TRACM",
        ["src/python_bindings.cpp"],
        # Example: passing in the version to the compiled code
        define_macros = [('VERSION_INFO', find_version("tracm/__init__.py"))],
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args),
]


setup(
    name="tracm",
    version=find_version("tracm/__init__.py"),
    author="Gerry Tonkin-Hill",
    description=
    "A fast python and and c++ pipeline for identifying transmission clusters from single genome and metagenomic data",
    long_description_content_type="text/markdown",
    url="https://github.com/gtonkinhill/tracm",
    install_requires=[
        'numpy', 'scipy', 'plotly', 'pyfastx', 'datetime', 'tqdm', 'sourmash', 'joblib', 'ncbi_genome_download'
    ],
    python_requires='>=3.9.0',
    packages=['tracm'],
    keywords='transmission clustering metagenomics',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    entry_points={
        'console_scripts':
        ['tracm = tracm.__main__:main',],
    },
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
    # Currently, build_ext only provides an optional "highest supported C++
    # level" feature, but in the future it may provide more features.
    cmdclass={"build_ext": build_ext},
)




from setuptools import setup, find_packages
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


here = path.abspath(path.dirname(__file__))

setup(
    name="mtrans",
    version=find_version("mtrans/__init__.py"),
    author="Gerry Tonkin-Hill",
    description=
    "A fast python and and c++ pipeline for clustering sequences using the 'transcluster' method",
    long_description_content_type="text/markdown",
    url="https://github.com/gtonkinhill/mtrans",
    install_requires=[
        'numpy', 'scipy', 'plotly', 'pyfastx', 'datetime', 'numba',
        'pysam', 'pysamstats', 'tqdm'
    ],
    python_requires='>=3.6.0',
    packages=['mtrans'],
    keywords='transmission clustering pathogen covid sars-cov-2 epidemiology',
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
        ['mtrans = mtrans.__main__:main',
        ['bampileup = mtrans.bampileup:main']],
    },
)

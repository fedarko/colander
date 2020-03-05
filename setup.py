#!/usr/bin/env python
# NOTE: This file is derived from Qeeseburger's setup.py file.

from setuptools import find_packages, setup

classes = """
    Development Status :: 3 - Alpha
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3.6
    Programming Language :: Python :: 3.7
    Programming Language :: Python :: 3 :: Only
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split("\n") if s]

description = "Strain variation estimation from a metagenomic assembly graph"

with open("README.md") as f:
    long_description = f.read()

setup(
    name="colander",
    version="0.0.0",
    license="MIT",
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Marcus Fedarko",
    author_email="mfedarko@ucsd.edu",
    maintainer="Marcus Fedarko",
    maintainer_email="mfedarko@ucsd.edu",
    url="https://github.com/fedarko/colander",
    packages=find_packages(),
    include_package_data=True,
    install_requires=["networkx"],
    # Based on how Altair splits up its requirements:
    # https://github.com/altair-viz/altair/blob/master/setup.py
    extras_require={
        "dev": ["pytest", "pytest-cov", "flake8", "black"]
    },
    classifiers=classifiers,
    zip_safe=False,
    python_requires=">=3.6",
)

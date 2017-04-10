#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name="signalAlign",
      version="0.1.3",
      description="A library for signal-level analysis of ONT data",
      author="Art Rand",
      author_email="arand@soe.ucsc.edu",
      url="https://github.com/ArtRand/signalAlign",
      package_dir={"": "src"},
      packages=find_packages("src"),
      install_requires=["numpy==1.9.2",
                        "marginAlign==1.1.9",
                        "h5py>=2.2.1",
                        "pysam==0.8.2.1",
                        "pandas==0.18.1",
                        "sonLib==1.1.0",
                        "PyYAML==3.12"]
      )

## SignalAlign

#### MinION signal-level alignment and methylation detection using hidden Markov Models with hierarchical Dirichlet process kmer learning.

### Introduction
Nanopore sequencing is based on the principal of isolating a nanopore in a membrane separating buffered salt solutions, then applying a voltage across the membrane and monitoring the ionic current through the nanopore. The Oxford Nanopore Technologies (ONT) MinION sequences DNA by recording the ionic current as DNA strands are enzymatically guided through the nanopore. **SignalAlign** will align the ionic current from the MinION to a reference sequence using a trainable hidden Markov model (HMM). The emissions model for the HMM can either be the table of parametric normal distributions provided by ONT or a hierarchical Dirichlet process (HDP) mixture of normal distributions. The HDP models enable mapping of methylated bases to your reference sequence. Instructions for usage including building/training HDP models can be found in the [manual](https://github.com/ArtRand/signalAlign/blob/master/Manual.md).

### Requirements
* Python 2.7 (most packages are included in [Anaconda](https://www.continuum.io/downloads) distribution of Python 2.7)
    1. H5Py
    2. Numpy
    3. Pandas
* BWA-MEM (Li H. (2013), instructions can be found (https://github.com/lh3/bwa))
    * Needs to be in path
* GCC 4.4.7 or newer (tested on 4.4.7 and 5.0)

### Installation:
1. Recursively clone this repo `git clone --recursive https://github.com/ArtRand/signalAlign.git`
2. cd into the directory and run `make`
3. Test the installation by running `make test`
4. All of the programs can be found in the `/bin/` directory

*Code in this repo is based on cPecan (https://github.com/benedictpaten/cPecan)*

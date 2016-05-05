## SignalAlign Manual

#### Art Rand
Updated 5/2/16

### Introduction
SignalAlign is a hidden Markov model (HMM) software package for aligning ionic current from the Oxford Nanopore Technologies (ONT) MinION to a reference sequence and inferring properties of the reference sequence.

### Installation
1. Recursively clone this repo `git clone --recursive https://github.com/ArtRand/signalAlign.git`
2. cd into the directory and run `make`
3. Test the installation by running `make test`
4. All of the programs can be found in the `/bin/` directory

### Quick Start
1. Install
2. `cd bin`
3. Run:
```bash
./runSignalAlign -d /path/to/minionReads/ -r /path/to/reference.fasta -o /path/to/output/ 2> /dev/null
```
4. If you want to try the test _E. coli_ files run:
```bash
./runSignalAlign -d ../tests/minion_test_reads/ecoli/ -r ../tests/test_sequences/E.coli_K12.fasta -o ../tests/ 2> ../tests/test_run.err
```

### Description of programs
* runSignalAlign
    * Aligns ionic current events from a directory of base-called MinION reads (.fast5) to a reference sequence. With appropriate inputs, characters in the reference sequence can be flagged as _ambiguous_, meaning multiple bases will be probabilistically aligned to a position. Right now, the program only supports aligning multiple cytosine variants (methylation variants) to a position.
* trainModels
    * Trains the transitions and/or emissions of the HMM. Uses a directory of base-called reads and a reference sequence to learn transition parameters to the model. Once enough assignments have been accumulated it will (optionally) rebuild a hierarchical Dirichlet process (HDP) from these assignments.
* hdp_pipeline
    * Recommended method for building HDP models. Uses a directory of _alignments_ to make a new HDP model.

### runSignalAlign
#### Input
*   A directory of MinION reads (*.fast5) that have been basecalled. Right now, Metrichor versions 1.15.0 and 1.19.0.
*   A reference sequence in FASTA format.
_Optional_
*   A file containing trained HMM transitions parameters.
*   A file containing a HDP model or other emissions (normal distributions) parameters.
*   Target regions file. Only reads that map to these regions will follow on to event-alignment
*   Ambiguity positions file. T

#### Options
#### Output

#### Using substitution and target files
#### Calling methylation

### trainModels
#### Input
#### Options
#### Output

### hdp_pipeline
#### Input
#### Options
#### Output

### File formats

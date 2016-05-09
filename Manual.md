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
* A directory of MinION reads (*.fast5) that have been basecalled. Right now, Metrichor versions 1.15.0 and 1.19.0.
* A reference sequence in FASTA format.
* Output location, a path to use as working directory. A new directory will be made here, so the program won't pollute this directory.

_Optional_
* A file containing trained HMM transitions parameters.
* A file containing a HDP model or other emissions (normal distributions) parameters.
* Target regions file. Only reads that map to these regions will follow on to event-alignment
* Ambiguity positions file, specifies which positions in the reference to flag as 'methylation ambiguous' (used for methylation variant calling)

#### Options
* `--stateMachineType, -smt` HMM to use. Options: `threeState` and `threeStateHdp`.
* `--twoWay`, specify for C vs mC classification. Three-way classification (C, mC, hmC) by default.
* `--threshold, -t`. Minimum posterior match probability threshold (matches below this threshold will not be tabulated). Default: 0.01.
* `--diagonalExpansion, -e` Mumber of diagonals to expand around each anchor, Default: 50.
* `--constraintTrim, -m` Amount to remove from an anchor constraint. Default: 14.
* `--un_banded, -ub` Remove all anchors, compute entire DP matrix. Default: False.
* `--jobs, -j` Number of jobs to run concurrently, Default: 4.
* `--nb_files, -n` Maximum number of files to align. Default: 500.
* `--sparse_output` Sparse output option (see output section)

#### Output

There are two output formats. Normal and sparse. Normal output has the following tab-separated-format:

| Contig | Reference Index | Reference 6-mer | Read File | Strand | Event Index | Event Mean | Event Noise | Event Duration | Aligned 6-mer | Scaled Mean Current | Scaled Noise | Posterior Probability | Descaled Event Mean | Model (ONT) Mean | Path 6-mer |
|--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |

And sparse format:

| Contig | Reference Index | Reference 6-mer | Read File | Strand | Event Index | Aligned 6-mer | Posterior Probability | Path 6-mer |
|--- | --- | --- | --- | --- | --- | --- | --- | --- |

Each read that aligns to the target regions (if you specified that option) will come as a separate file.

#### Mapping cytosine methylation by using substitution and target files
**signalAlign** uses two files to map methylation in a reference sequence. First it uses a *target file* that means only reads that map to regions specified in the file will be aligned on the signal-level. The *target file* has format: `start \t end \t sequence`. An example can be seen here:`/signalAlign/tests/test_regions/test_sites_bal_1.tgt`. The second file is a *label file* that tells signalAlign which bases to flag as ambiguous. This file has the format `X \t position \t ... \n` one line for both the template and complement strand. An example is at `/signalAlign/tests/test_regions/test_labels_bal_1.tsv`.

An example command that would produce 2-way methylation probabilities in _E. coli_
(n.b. to run this example, download the HDP models from the [HDP_models](https://github.com/ArtRand/HDP_models) repo)
```bash
./runSignalAlign \
> -d=../tests/minion_test_reads/ecoli/ \
> -r=../tests/test_sequences/E.coli_K12.fasta \
> -T=../../HDP_models/ecoli_models/template_trained.hmm \
> -C=../../HDP_models/ecoli_models/complement_trained.hmm \
> -tH=../../HDP_models/ecoli_models/template.multisetPrior2.nhdp \
> -cH=../../HDP_models/ecoli_models/complement.multisetPrior2.nhdp \
> --twoWay \
> -smt=threeStateHdp \
> -s \
> -o=../../ \
> 2> ../../a.err
```

### trainModels
#### Input

#### Options
#### Output

### hdp_pipeline
#### Input
#### Options
#### Output

### File formats

## SignalAlign Manual

**UPDATED for R9 refactor**
Updated 12/19/16 - Working to integrate new matrichor base caller, and multi-contig reference sequences

### Introduction
SignalAlign is a hidden Markov model (HMM) software package for aligning the ionic current signal from the Oxford Nanopore Technologies (ONT) MinION to a reference sequence and inferring properties of the reference sequence.

### Installation
1. Recursively clone this repo `git clone --recursive https://github.com/ArtRand/signalAlign.git`
2. Make a virtual environment in the directory `virtualenv venv`
3. Activate it `. venv/bin/activate`
4. Pip install the requirements `pip install -r requirements.txt`
5. run `make`
6. Test the installation by running `make test`
7. All of the programs can be found in the `/bin/` directory
8. If python can't find the modules in the repo, add the directory to the PYTHONPATH like so: `export PYTHONPATH="$PYTHONPATH:$(pwd)"`

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
    * Aligns ionic current events from a directory of basecalled MinION reads (.fast5) to a reference sequence. With appropriate inputs, characters in the reference sequence can be flagged as _ambiguous_, meaning multiple bases will be probabilistically aligned to a position. Right now, the program the program supports aligning cytosine variants (5-mC and 5-hmC) and adenine variants (6-mA) to a position.
* trainModels
    * Trains the transitions and/or emissions of the HMM. Uses a directory of basecalled reads and a reference sequence to learn transition parameters to the model. Once enough assignments have been accumulated it will (optionally) rebuild a hierarchical Dirichlet process (HDP) from these assignments.
* hdp_pipeline
    * Recommended method for building HDP models. Uses a directory of _alignments_ to make a new HDP model, also see the [pipeline examples](https://github.com/ArtRand/signalAlign_notebook/blob/master/experiment_scripts/BuildModels.py#L512) for detailed instructions. 

### runSignalAlign
#### Input
_Required_
* A directory of MinION reads, `-d`, (*.fast5) that have been basecalled. Currently, the following versions of Metrichor basecalled files are supported: 1.15.0, 1.19.0, 1.20.0, 1.22.2, 1.22.4. If you have a more recent or unsupported version open an issue and I'll modify the program (or feel free to implement it yourself and issue a pull request!). 
* A reference sequence, `-r`, in FASTA format.
* Output location, `-o`, a path to use as working directory. A new directory will be made here, so the program won't pollute this directory.

_Optional_
* A file containing trained HMM transitions parameters `-T` (template HMM) `-C` (complement HMM).
* A file containing a HDP model or other emissions (normal distributions) parameters `-tH` (template HDP) `-cH` (complement HDP).
* Target regions file. Only reads that map to these regions will follow on to event-alignment `-q`
* Ambiguity positions file, `-p` flags positions to be aligned to multiple bases (variant calling).

#### Options and flags
* `--in_tempalte_hmm`, `-T` template HMM parameters file
* `--in_complement_hmm`, `-C` complement HMM parameters file 
* `--template_hdp`, `-tH` template HDP model file
* `--complement_hdp`, `-cH` complement HDP model file 
* `--degenerate`, `-x` nucleotide options for degenerate or _ambiguous_ positions. `variant` = {A,C,G,T} `cytosine2` = {CE} `cytosine3` = {CEO} `adenosine` = {AI}. **n.b.** E = 5-methylcytosine, O = 5-hydroxymethylcytosine, I = 6-methyladenine
* `--stateMachineType`, `-smt` HMM to use. Options: `threeState` and `threeStateHdp`. Default: `threeState`.
* `--file_of_files`, `-fofn` a file containing the absolute path to files to align with, one file path per line
* `--threshold`, `-t`. Minimum posterior match probability threshold (matches below this threshold will not be tabulated). Default: 0.01.
* `--diagonalExpansion`, `-e` Mumber of diagonals to expand around each anchor, Default: 50.
* `--constraintTrim`, `-m` Amount to remove from an anchor constraint. Default: 14.
* `--target_regions`, `-q` File containing target regions to align to, if the read doesn't get mapped to a region in this file the signal-level alignment will not precede. The format is `start \t end \t kmer` where `start` and `end` are the genomic coordinates and the `kmer` is the sequence at those coordinates. The `kmer` is mostly historical, so it doesn't have to be perfect. See below for details.
* `--ambiguity_positions`, `p` file containing positions to flag as ambiguous see below for details.
* `--jobs, -j` Number of jobs to run concurrently, Default: 4.
* `--nb_files`, `-n` Maximum number of files to align, will randomly select if you point to a directory/fofn with more files than this so if you want to make sure you use all of the files set a large number (or count the lines in your fofn). Default: 500.
* `--ambig_char`, `-X` in development, will be for specifying specific subsets of ambiguous positions **leave as default**
* `--output_format`, `-f` format of output (this is different than the summary that goes to `stdout`). Options are: `full`, `variantValler`, and `assignments`, see below for description.
* `--output_location`, `-o` place to put results, it will make a new working directory in this location
* `--error_correct`, not implemented, yet.

#### Output

There are three output formats. `full`, `variantCaller`, and `assignments`. Each read that aligns to the target regions (if you specified that option) will come as a separate file.
 
`full` has the following tab-separated-format:

| Contig | Reference Index | Reference k-mer | Read File | Strand | Event Index | Event Mean | Event Noise | Event Duration | Aligned k-mer | Scaled Mean Current | Scaled Noise | Posterior Probability | Descaled Event Mean | Model (ONT) Mean | Path k-mer |
|--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |

`variantCaller` has the following tab-separated format:

| Event Index | Reference Position | Base | Posterior Probability | Strand | Forward Mapped | Read File |
|--- | --- | --- | --- | --- | --- | --- |

finally, `assignments` has the following tab-separated format: 

| k-mer | Read File | Descaled Event Mean | Posterior Probability | 
|--- | --- | --- | --- |

### trainModels
#### Input
_Required_
* A directory of MinION reads, `-d`, (*.fast5) that have been basecalled.
* A reference sequence, `-r`, in FASTA format.
* Output location, `-o`, a path to use as working directory. A new directory will be made here, so the program won't pollute this directory.

_Optional_
* A file containing trained HMM transitions parameters `-T` (template HMM) `-C` (complement HMM).
* A file containing a HDP model or other emissions (normal distributions) parameters `-tH` (template HDP) `-cH` (complement HDP).

#### Options
* `--stateMachineType`, `-smt` HMM to use. Options: `threeState` and `threeStateHdp`. Default: `threeState`.
* `--file_of_files`, `-fofn` a file containing the absolute path to files to align with, one file path per line
* `--iterations, -i` number of iterations to perform. Default: 10.
* `--train_amount, -a` batch size _in bases_. Default: 15 kb (15000)
* `--emissions` flag, train emissions. Default: false (see note below about training HDP emissions)
* `--transitions` flag, train transitions
* `--threshold, -t` minimum posterior match probability threshold (matches below this threshold will not be tabulated). Default: 0.01. If you're just training the transitions (_recommended_) this won't matter.
* `--diagonalExpansion, -e` number of diagonals to expand around each anchor, Default: 50.
* `--constraintTrim, -m` amount to remove from an anchor constraint. Default: 14.
* `--in_tempalte_hmm`, `-T` template HMM parameters file
* `--in_complement_hmm`, `-C` complement HMM parameters file 
* `--template_hdp`, `-tH` template HDP model file
* `--complement_hdp`, `-cH` complement HDP model file 
* `--jobs, -j` Number of jobs to run concurrently, Default: 4.
* `--test` used for internal CI test
* `--ambiguity_positions`, `p` file containing positions to flag as ambiguous see below for details.
* `--ambig_char`, `-X` base to substitute in at positions indicated in _positions file_ this is used for supervised learning, we change the base at the positions for signal-level alignment. example `-X C E` would mean substitute the positions in the first batch of files with `C` and substitute the second batch with `E`. This way we can train the model on multiple base modifications at once. 
**HDP training options** It's recommended to use the pipeline described below instead of these options. 
* `--verbose` flag, verbose Gibbs sampling output, default: False
* `--samples, -s` number of samples to collect during Gibbs sampling
* `--thinning, -th` number of samples to thin (pick every `x` samples)
* `--min_assignments` Do not initiate Gibbs sampling unless this many assignments have been accumulated
#### Output
Trained models, within the working `/tempFiles_expectations/` will have:
* `template_trained.hmm` template HMM with trained parameters
* `complement_trained.hmm` complement HMM with trained parameters
If HDPs were used, a copy of the input HDP will also be here, they have `.nhdp` suffix.

### hdp_pipeline
#### Input
* **Recommended** use a build alignment. There are scripts in the /signalAlign/scripts folder and examples in the [HDP_models](https://github.com/ArtRand/HDP_models) repo to do this.
* Alternatively, you can use alignments in `full` format. In the case of using plain alignments, `-C`, `-mC`, and `-hmC` will globally change all of the cytosine bases to methyl or hydroxy-methyl cytosines. You may specify just `-C` if you don't want to label methylated cytosines. This can be used for canonical base training. 

#### Options
* `--number_of_assignments`, `-n` number of k-mer/event mean assignments to collect for each new label (C, mC, hmC).
* `--build_alignment` a concatenated (multi read) alignment in the `full` format to extract assignments from. This is the recommended way to use this program.
* `--threshold, -t` only take assignments (from the alignment) with posterior probability >= to this
* `--hdp_type` HDP type Options: `Prior`, `Fixed`, `twoWay` Default: `Prior`. `twoWay` is a Prior-type model that uses a gamma-distribution over the hyperparameters (concentration parameters) of the HDP.
* `--template_model, -tM` Template lookup table for priors on k-mers
* `--complement_model, -cM` Template lookup table for priors on k-mers
* `-B`, `-M`, `-L` Base, middle, and leaf concentration parameters, respectively. Default: 1
* `-Ba`, `-Bb,`, `-Ma`, `-Mb`, `-La`, `-Lb` Base alpha(`a`)/beta(`b`) parameters for gamma distribution prior over concentration parameters. Default: 1.
* `--samples, -s` number of samples to collect during Gibbs sampling. Default: 10000.
* `--thinning, -th` number of samples to thin (pick every `x` samples). Default: 100.
* `--verbose` flag, verbose Gibbs sampling output, default: False.
* `--grid_start` sampling grid start (in pA). Default: 30 pA.
* `--grid_end` sampling grid end (in pA). Default: 90 pA.
* `--grid_length` number of nodes in the sampling grid. Default 1200.
* `--out, -o` output location

#### Output
Trained HDP models. Script will make all of the models for 3-way classification and the singleLevel and multiset models for 2-way classification.

Example command:
```bash
./hdp_pipeline --build_alignment=../../HDP_models/ecoli_models/buildAlignment_PCR_andLabeled_b1.tsv \
> -tM ../models/testModel_template.model \
> -cM ../models/testModel_complement.model \
> -Ba 1 -Bb 1 -Ma 1 -Mb 1 -La 1 -Lb 1 -s 15000 --verbose --hdp_type twoWay \
> -o ./hdp/
2> ./hdp/a.err
```

#### Mapping methylation by using substitution and target files
**signalAlign** uses two files to map methylation in a reference sequence. First it uses a *target file* that specifies positions in the reference, only reads that map to regions specified in the file will be aligned on the signal-level. The *target file* has format: `start \t end \t sequence`. An example can be seen here:`/signalAlign/tests/test_regions/test_sites_bal_1.tgt`. The second file is a *label file* that tells signalAlign which bases to flag as ambiguous. This file has the format `X \t position \t ... \n` one line for both the template and complement strand. An example is at `/signalAlign/tests/test_regions/test_labels_bal_1.tsv`.

An example command that would produce methylation call probabilities in _E. coli_
(**n.b.** to run this example, download the HDP models from the [HDP_models](https://github.com/ArtRand/HDP_models) repo)
```bash
./runSignalAlign \
> -d ../tests/minion_test_reads/ecoli/ \
> -r ../tests/test_sequences/E.coli_K12.fasta \
> -T ../../HDP_models/ecoli_r7.3_models/template_trained.hmm \
> -C ../../HDP_models/ecoli_r7.3_models/complement_trained.hmm \
> -tH ../../HDP_models/ecoli_r7.3_models/template.multisetPrior2.nhdp \
> -cH ../../HDP_models/ecoli_r7.3_models/complement.multisetPrior2.nhdp \
> -x cytosine2 \
> -smt=threeStateHdp \
> -q ../tests/test_regions/test_sites_bal_1.tgt \
> -p ../tests/test_regions/test_labels_bal_1.tsv \
> -f variantCaller \
> -o ../../ \
> 2> ../../a.err
```

### Example Pipelines
Pipeline scripts that correctly run all of these programs can be found at the [signalAlign_notebook](https://github.com/ArtRand/signalAlign_notebook/tree/master/experiment_scripts) repo. Specifically, `BuildModels.py` takes as input the pcr-amplified reads, genomic reads, a reference, and will run the full EM pipeline on all `GATC` or `CCWGG` positions in your reference. There are, of course, other options too.

##signalAlign
####Hidden Markov Models with Hierarchal Dirichlet Process kmer learning
Code in this repo is based on cPecan (https://github.com/benedictpaten/cPecan)

### Getting Around:
Most of the scripts and programs you want are in ../sonLib/bin make sure you can find:
>1. runSignalAlign 
>2. hdp_pipeline
>3. compareDistributions

##### runSignalAlign
This is a python script that will produce event-to-nucleotide alignments. Most of the flags are described if you do
`./runSignalAlign -h`
A basic command that would align MinION reads to a reference would be:
```
$ ./runSignalAlign \
> -d=../../signalAlign/tests/minion_test_reads/C/ \
> -r=../../signalAlign/tests/test_sequences/zymo_sequence.fasta \
> -smt=threeState -o=../../temp/ \
> 2> ../../temp/a.err
```
This will align the `.fast5` files in `../../signalAlign/tests/minion_test_reads/C` to the zymo reference sequence
To use an HDP model, the only thing you need change is the -smt (stateMachineType) and add the paths to the HDP files:
```
$ ./runSignalAlign \
> -d=../../signalAlign/tests/minion_test_reads/C/ \
> -r=../../signalAlign/tests/test_sequences/zymo_sequence.fasta \
> -smt=threeStateHdp \
> -tH=/projects/nanopore-working/zymo/signalAlign_results/03_03_hdpPipeline/multisetPrior_11_11_11_10K/template.multisetPrior.nhdp \
> -cH=/projects/nanopore-working/zymo/signalAlign_results/03_03_hdpPipeline/multisetPrior_11_11_11_10K/complement.multisetPrior.nhdp \
> -o=../../temp/ \
> 2> ../../temp/a.err
```
##### hdp_pipeline
This script will build a hdp from alignments.  Again check the flags with`./hdp_pipeline -h`.  A typical use would be:
```
$ ./hdp_pipeline \
> -C=/projects/nanopore-working/zymo/signalAlign_results/02_26_sm3_ecoliTrained_defBand/C/tempFiles_alignment/*.tsv 
> -mC=/projects/nanopore-working/zymo/signalAlign_results/02_26_sm3_ecoliTrained_defBand/mC/tempFiles_alignment/*.tsv 
> -hmC=/projects/nanopore-working/zymo/signalAlign_results/02_26_sm3_ecoliTrained_defBand/hmC/tempFiles_alignment/*.tsv 
> -n=10000 --hdp_type=multisetPrior -Ba=1 -Bb=1 -Ma=1 -Mb=1 -La=1 -Lb=1 -s=10000 -th=200 --verbose --no_train 
> -o=/projects/nanopore-working/zymo/signalAlign_results/03_03_hdpPipeline/multisetPrior_11_11_11_10K/
```
##### compare_distributions
This C program will evaluate the nanopore HDP over a reasonable range and output the kmer densities in a directory. You **should** make a 
new directory for the files (there will be 46,656 of them). Usage:
```
$ ./compareDistributions \
> /projects/nanopore-working/zymo/signalAlign_results/03_03_hdpPipeline/multisetPrior_11_11_11_10K/template.multisetPrior.nhdp \
> /projects/nanopore-working/zymo/signalAlign_results/03_03_hdpPipeline/multisetPrior_11_11_11_10K/distributions_template 
```


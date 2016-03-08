##signalAlign
####Hidden Markov Models with Hierarchal Dirichlet Process kmer learning
Code in this repo is based on cPecan (https://github.com/benedictpaten/cPecan)

### Getting Around:
Most of the scripts and programs you want are in ../sonLib/bin make sure you can find:
    1. runSignalAlign 
    2. hdp_pipeline
    3. compareDistributions

##### runSignalAlign
This is a python script that will produce event-to-nucleotide alignments. Most of the flags are described if you do
`./runSignalAlign -h`
A basic command that would align MinION reads to a reference would be:
```
./runSignalAlign \
> -d=../../signalAlign/tests/minion_test_reads/C/ 
> -r=../../signalAlign/tests/test_sequences/zymo_sequence.fasta \
> -smt=threeState -o=../../temp/ \
> 2> ../../temp/a.err

```


### Quick tutorial 
#### Generating a new HDP from some alignments (HDP or otherwise)

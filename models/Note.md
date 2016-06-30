### Note on files in this directory
There are two kinds of files in this directory, models `*.model`, and HMMs `*.hmm`. All of the files in this directory are special formats meant to work with signalAlign. The format for `model` files is:
```
line 0: [correlation coefficient] [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ] (.../kmer) \n
line 1: [correlation coefficient] [level_mean] [level_sd, scaled] [noise_mean] [noise_sd] [noise_lambda ] (.../kmer) \n
```

The files `template_median68pA.model`, `complement_median68pA.model`, and `complement_median68pA_pop2.model` were extracted from R7.3 MinION reads and have entries for 4096 canonical 6-mers. The `testModel_*` files have entries for 46,656 6-mers (accounting for a 6-letter alphabet that includes methylation). The 6-mer distributions are not informed, however, meaning that the distribution for `AGGCTA` and `AGGETA` where `E` is 5-methylcytosine will be the same in these models. See the manual for more information on learning new 6-mer distributions.

The `hmm` files have the format:
```
line 0: type \t stateNumber \t symbolSetSize \t hasModel \n
line 1: [transitions... \t] likelihood \n
line 2: [Event Model] \n
line 3: [Event Expectations] \n
line 4: [Event Posteriors] \n
line 5: [Observed] \n
```
All of the transitions are set to default naive parameters. For more information on learning new transition parameters, see the manual.
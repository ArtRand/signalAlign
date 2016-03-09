#ifndef EMISSIONS_MATRIX_H
#define EMISSIONS_MATRIX_H

#define MATRIX_SIZE 16777216 // assuming 4096 kmers

void emissions_kmer_setMatchProbsToDefaults(double *emissionMatchProbs);

void emissions_kmer_setGapProbsToDefaults(double *emissionGapProbs);

#endif

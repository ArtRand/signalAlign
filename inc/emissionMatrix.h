#ifndef EMISSIONS_MATRIX_H
#define EMISSIONS_MATRIX_H

#define KMER_LENGTH 6
#define NUM_OF_KMERS 4096    // for alphabet 'AGCT', may not need 'N' in there
#define MATRIX_SIZE 16777216 // assuming 4096 kmers

void emissions_kmer_setMatchProbsToDefaults(double *emissionMatchProbs);

void emissions_kmer_setGapProbsToDefaults(double *emissionGapProbs);

#endif

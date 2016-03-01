// Compare nanopore HDP distributions

#include <stdio.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include "sonLib.h"
#include "nanopore_hdp.h"
#include "hdp_math_utils.h"

#define FILENAME_BUFFER_LEN 100

char* kmer_from_index(int64_t index, const char* alphabet, int64_t alphabet_size, int64_t kmer_length) {
    char* kmer = (char*) malloc(sizeof(char) * (kmer_length + 1));
    int64_t index_remainder = index;
    kmer[kmer_length] = '\0';
    for (int64_t i = kmer_length - 1; i >= 0; i--) {
        kmer[i] = alphabet[index_remainder % alphabet_size];
        index_remainder /= alphabet_size;
    }
    return kmer;
}

void write_kmer_distr(NanoporeHDP* nhdp, char* kmer, double* eval_grid, int64_t grid_length, char *workingDirectory) {
    //char filename[FILENAME_BUFFER_LEN];
    char *filename = stString_print("%s/%s_distr.txt", workingDirectory, kmer);
    //sprintf(filename, "%s/%s_distr.txt", workingDirectory, kmer);
    FILE* out = fopen(filename, "w");

    //fprintf(out, "%s\t", kmer);

    double density;
    for (int64_t i = 0; i < grid_length - 1; i++) {
        density = get_nanopore_kmer_density(nhdp, kmer, eval_grid + i);
        fprintf(out, "%.17lg\n", density);
    }

    density = get_nanopore_kmer_density(nhdp, kmer, eval_grid + (grid_length - 1));
    fprintf(out, "%.17lg\n", density);
    
    fclose(out);
}

void write_all_kmer_distrs(NanoporeHDP* nhdp, double* eval_grid, int64_t grid_length, char *workingDirectory) {
    char *x_filename = stString_print("%s/x_vals.txt", workingDirectory);
    //char *distributions = stString_print("%s/distributions.tsv", workingDirectory);
    //sprintf(x_filename, "%s/x_vals.txt", workingDirectory);
    //sprintf(distributions, "%s/distributions.tsv", workingDirectory);
    FILE* x_vals = fopen(x_filename, "w");
    //FILE *distributionFileHandle = fopen(distributions, "a");

    for (int64_t i = 0; i < grid_length - 1; i++) {
        fprintf(x_vals, "%.17lg\n", eval_grid[i]);
    }

    fprintf(x_vals, "%.17lg", eval_grid[grid_length - 1]);

    fclose(x_vals);
    
    int64_t alphabet_size = get_nanopore_hdp_alphabet_size(nhdp);
    int64_t kmer_length = get_nanopore_hdp_kmer_length(nhdp);
    char* alphabet = get_nanopore_hdp_alphabet(nhdp);
    
    int64_t num_kmers = power(alphabet_size, kmer_length);
    
    for (int64_t kmer_index = 0; kmer_index < num_kmers; kmer_index++) {
        char* kmer = kmer_from_index(kmer_index, alphabet, alphabet_size, kmer_length);
        
        write_kmer_distr(nhdp, kmer, eval_grid, grid_length, workingDirectory);
        
        free(kmer);
    }

    free(alphabet);
    //fclose(workingDirectory);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "USAGE_NEW: compareDistributions [NanoporeHDP_file] [distribution_directory]\n");
        exit(EXIT_FAILURE);
    }

    char *modelFile = argv[1];
    fprintf(stderr, "[compareDistributions] NOTICE: Loading NanoporeHDP from %s\n", modelFile);

    char *workingDirectory = argv[2];
    //FILE *fH = fopen(destinationFile, "w");
    fprintf(stderr, "[compareDistributions] NOTICE: Putting distributions in %s\n", workingDirectory);

    NanoporeHDP *nHdp= deserialize_nhdp(modelFile);

    double gridStart = 30.0;
    double gridStop = 90.0;
    int64_t gridLength = 600;

    double *evaluationGrid = linspace(gridStart, gridStop, gridLength);

    write_all_kmer_distrs(nHdp, evaluationGrid, gridLength, workingDirectory);

    free(evaluationGrid);
    destroy_nanopore_hdp(nHdp);

    return 0;
}
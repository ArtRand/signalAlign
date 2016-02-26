//
//  nanopore_hdp.h
//  
//
//  Created by Jordan Eizenga on 1/8/16.
//
//

#ifndef nanopore_hdp_h
#define nanopore_hdp_h

#include <stdbool.h>
#include <inttypes.h>
#include "hdp.h"

typedef struct _nanoporeHDP {
    HierarchicalDirichletProcess* hdp;
    char* alphabet;
    int64_t alphabet_size;
    int64_t kmer_length;
    stSet* distr_metric_memos;
} NanoporeHDP;

typedef enum _nanoporeHdpType {
    singleLevelFixed = 0,
    singleLevelPrior = 1,
    multisetFixed = 2,
    multisetPrior = 3,
} NanoporeHdpType;

typedef struct _nanoporeDistributionMetricMemo {
    NanoporeHDP* nhdp;
    DistributionMetricMemo* memo;
} NanoporeDistributionMetricMemo;

void destroy_nanopore_hdp(NanoporeHDP* nhdp);

int64_t get_nanopore_hdp_kmer_length(NanoporeHDP* nhdp);
int64_t get_nanopore_hdp_alphabet_size(NanoporeHDP* nhdp);
char* get_nanopore_hdp_alphabet(NanoporeHDP* nhdp);

void execute_nhdp_gibbs_sampling(NanoporeHDP* nhdp, int64_t num_samples, int64_t burn_in, int64_t thinning, bool verbose);

void execute_nhdp_gibbs_sampling_with_snapshots(NanoporeHDP* nhdp, int64_t num_samples, int64_t burn_in, int64_t thinning,
                                                void (*snapshot_func)(HierarchicalDirichletProcess*, void*),
                                                void* snapshot_func_args, bool verbose);

void finalize_nhdp_distributions(NanoporeHDP* nhdp);

double get_nanopore_kmer_density(NanoporeHDP* nhdp, void *kmer, void *event);

void update_nhdp_from_alignment(NanoporeHDP* nhdp, const char* alignment_filepath, bool has_header);

// filter for only observations containing "strand_filter" in the strand column
void update_nhdp_from_alignment_with_filter(NanoporeHDP* nhdp, const char* alignment_filepath,
                                            bool has_header, const char* strand_filter);

// computing metrics on distributions

double get_kmer_distr_distance(NanoporeDistributionMetricMemo* memo, char* kmer_1, char* kmer_2);

NanoporeDistributionMetricMemo* new_nhdp_kl_divergence_memo(NanoporeHDP* nhdp);
NanoporeDistributionMetricMemo* new_nhdp_hellinger_distance_memo(NanoporeHDP* nhdp);
NanoporeDistributionMetricMemo* new_nhdp_l2_distance_memo(NanoporeHDP* nhdp);
NanoporeDistributionMetricMemo* new_nhdp_shannon_jensen_distance_memo(NanoporeHDP* nhdp);
// note: the lifetime of a NanoporeDistributionMetricMemo is tied to the lifetime of the
// NanoporeHDP that generated it

double compare_nhdp_distrs_kl_divergence(NanoporeHDP* nhdp_1, char* kmer_1,
                                         NanoporeHDP* nhdp_2, char* kmer_2);

double compare_nhdp_distrs_l2_distance(NanoporeHDP* nhdp_1, char* kmer_1,
                                       NanoporeHDP* nhdp_2, char* kmer_2) ;

double compare_nhdp_distrs_shannon_jensen_distance(NanoporeHDP* nhdp_1, char* kmer_1,
                                                   NanoporeHDP* nhdp_2, char* kmer_2);

double compare_nhdp_distrs_hellinger_distance(NanoporeHDP* nhdp_1, char* kmer_1,
                                              NanoporeHDP* nhdp_2, char* kmer_2);


// single level HDP
NanoporeHDP* flat_hdp_model(const char* alphabet, int64_t alphabet_size, int64_t kmer_length, double base_gamma,
                            double leaf_gamma, double sampling_grid_start, double sampling_grid_stop,
                            int64_t sampling_grid_length, const char* model_filepath);
NanoporeHDP* flat_hdp_model_2(const char* alphabet, int64_t alphabet_size, int64_t kmer_length,
                              double base_gamma_alpha, double base_gamma_beta, double leaf_gamma_alpha,
                              double leaf_gamma_beta, double sampling_grid_start, double sampling_grid_stop,
                              int64_t sampling_grid_length, const char* model_filepath);

// second level of HDP based on number of purines/pyrimidines
NanoporeHDP* purine_composition_hdp_model(char* purine_alphabet, int64_t num_purines,
                                          char* pyrimidine_alphabet, int64_t num_pyrimidines,
                                          int64_t kmer_length, double base_gamma, double middle_gamma,
                                          double leaf_gamma, double sampling_grid_start, double sampling_grid_stop,
                                          int64_t sampling_grid_length, const char* model_filepath);
NanoporeHDP* purine_composition_hdp_model_2(char* purine_alphabet, int64_t num_purines,
                                            char* pyrimidine_alphabet, int64_t num_pyrimidines,
                                            int64_t kmer_length, double base_gamma_alpha, double base_gamma_beta,
                                            double middle_gamma_alpha, double middle_gamma_beta,
                                            double leaf_gamma_alpha, double leaf_gamma_beta, double sampling_grid_start,
                                            double sampling_grid_stop, int64_t sampling_grid_length,
                                            const char* model_filepath);

// second level of HDP based on multiset of nucleotides
NanoporeHDP* multiset_hdp_model(const char* alphabet, int64_t alphabet_size, int64_t kmer_length,
                                double base_gamma, double middle_gamma, double leaf_gamma,
                                double sampling_grid_start, double sampling_grid_stop, int64_t sampling_grid_length,
                                const char* model_filepath);
NanoporeHDP* multiset_hdp_model_2(const char* alphabet, int64_t alphabet_size, int64_t kmer_length,
                                  double base_gamma_alpha, double base_gamma_beta, double middle_gamma_alpha,
                                  double middle_gamma_beta, double leaf_gamma_alpha, double leaf_gamma_beta,
                                  double sampling_grid_start, double sampling_grid_stop, int64_t sampling_grid_length,
                                  const char* model_filepath);

// second level of HDP based on middle 2 nucleotides
NanoporeHDP* middle_2_nts_hdp_model(const char* alphabet, int64_t alphabet_size, int64_t kmer_length,
                                    double base_gamma, double middle_gamma, double leaf_gamma,
                                    double sampling_grid_start, double sampling_grid_stop, int64_t sampling_grid_length,
                                    const char* model_filepath);
NanoporeHDP* middle_2_nts_hdp_model_2(const char* alphabet, int64_t alphabet_size, int64_t kmer_length,
                                      double base_gamma_alpha, double base_gamma_beta, double middle_gamma_alpha,
                                      double middle_gamma_beta, double leaf_gamma_alpha, double leaf_gamma_beta,
                                      double sampling_grid_start, double sampling_grid_stop,
                                      int64_t sampling_grid_length, const char* model_filepath);


void serialize_nhdp(NanoporeHDP* nhdp, const char* filepath);
NanoporeHDP* deserialize_nhdp(const char* filepath);

void nanoporeHdp_buildNanoporeHdpFromAlignment(NanoporeHdpType type,
                                               const char *templateModelFile, const char* complementModelFile,
                                               const char *alignments,
                                               const char *templateHDP, const char *complementHDP);

// n^k
int64_t power(int64_t n, int64_t k);
//  ((n k))
int64_t multiset_number(int64_t n, int64_t k);
// get word as int array by lexicographic order index
int64_t* get_word(int64_t word_id, int64_t alphabet_size, int64_t word_length);
// get multiset of a word as sorted int array by lexicographic order indx
int64_t* get_word_multiset(int64_t word_id, int64_t alphabet_size, int64_t word_length);
// get multiset lexicographic index from a multiset as sorted int array
int64_t multiset_id(int64_t* multiset, int64_t length, int64_t alphabet_size);
// get lexicographic index of multiset from lexicographic index of word
int64_t word_id_to_multiset_id(int64_t word_id, int64_t alphabet_size, int64_t word_length);

int64_t word_id(int64_t* word, int64_t alphabet_size, int64_t word_length);
int64_t* kmer_to_word(char* kmer, char* alphabet, int64_t alphabet_size, int64_t kmer_length);
int64_t kmer_id(char* kmer, char* alphabet, int64_t alphabet_size, int64_t kmer_length);
int64_t standard_kmer_id(char* kmer, int64_t kmer_length);


#endif /* nanopore_hdp_h */

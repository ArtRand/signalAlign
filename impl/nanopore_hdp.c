//
//  nanopore_hdp.c
//  
//
//  Created by Jordan Eizenga on 1/8/16.
//
//

// in 0-based index
#define ALIGNMENT_KMER_COL 9
#define ALIGNMENT_STRAND_COL 4
#define ALIGNMENT_SIGNAL_COL 13
#define ASSIGNMENT_KMER_COL 0
#define ASSIGNMENT_STRAND_COL 1
#define ASSIGNMENT_SIGNAL_COL 2
// number of expected column in the two kinds of input tables
#define NUM_ALIGNMENT_COLS 15
#define NUM_ASSIGNMENT_COLS 4

#define MODEL_ROW_HEADER_LENGTH 0
#define MODEL_MEAN_ENTRY 0
#define MODEL_NOISE_ENTRY 1
#define MODEL_ENTRY_LENGTH 5

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "pairwiseAligner.h"
#include "hdp_math_utils.h"


NanoporeHDP* package_nanopore_hdp(HierarchicalDirichletProcess* hdp, const char* alphabet, int64_t alphabet_size,
                                  int64_t kmer_length) {
    NanoporeHDP* nhdp = (NanoporeHDP*) malloc(sizeof(NanoporeHDP));
    
    // copy and sort alphabet
    char* internal_alphabet = (char*) malloc(sizeof(char) * (alphabet_size + 1));
    for (int64_t i = 0; i < alphabet_size; i++) {
        internal_alphabet[i] = alphabet[i];
    }
    
    int64_t min_idx;
    char temp;
    for (int64_t i = 0; i < alphabet_size; i++) {
        min_idx = i;
        for (int64_t j = i + 1; j < alphabet_size; j++) {
            if (internal_alphabet[j] < internal_alphabet[min_idx]) {
                min_idx = j;
            }
        }
        temp = internal_alphabet[i];
        internal_alphabet[i] = internal_alphabet[min_idx];
        internal_alphabet[min_idx] = temp;
    }
    
    for (int64_t i = 1; i < alphabet_size; i++) {
        if (alphabet[i - 1] == alphabet[i]) {
            fprintf(stderr, "Characters of alphabet must be distinct.\n");
            exit(EXIT_FAILURE);
        }
    }
    
    internal_alphabet[alphabet_size] = '\0';
    
    nhdp->hdp = hdp;
    nhdp->alphabet = internal_alphabet;
    nhdp->alphabet_size = alphabet_size;
    nhdp->kmer_length = kmer_length;
    
    // note: destroying the HDP housed in the NHDP will destroy the DistributionMetricMemo
    nhdp->distr_metric_memos = stSet_construct2(&free);
    
    return nhdp;
}

void destroy_nanopore_hdp(NanoporeHDP* nhdp) {
    destroy_hier_dir_proc(nhdp->hdp);
    stSet_destruct(nhdp->distr_metric_memos);
    free(nhdp->alphabet);
    free(nhdp);
}

int64_t get_nanopore_hdp_kmer_length(NanoporeHDP* nhdp) {
    return nhdp->kmer_length;
}


int64_t get_nanopore_hdp_alphabet_size(NanoporeHDP* nhdp) {
    return nhdp->alphabet_size;
}

char* get_nanopore_hdp_alphabet(NanoporeHDP* nhdp) {
    char* alphabet = nhdp->alphabet;
    int64_t alphabet_size = nhdp->alphabet_size;
    char* copy = (char*) malloc(sizeof(char) * (alphabet_size + 1));
    for (int64_t i = 0; i < alphabet_size; i++) {
        copy[i] = alphabet[i];
    }
    copy[alphabet_size] = '\0';
    return copy;
}


// wrappers
void execute_nhdp_gibbs_sampling(NanoporeHDP* nhdp, int64_t num_samples, int64_t burn_in,
                                 int64_t thinning, bool verbose) {
    execute_gibbs_sampling(nhdp->hdp, num_samples, burn_in, thinning, verbose);
}

void execute_nhdp_gibbs_sampling_with_snapshots(NanoporeHDP* nhdp,
                                                int64_t num_samples, int64_t burn_in, int64_t thinning,
                                                void (*snapshot_func)(HierarchicalDirichletProcess*, void*),
                                                void* snapshot_func_args, bool verbose) {
    execute_gibbs_sampling_with_snapshots(nhdp->hdp, num_samples, burn_in, thinning, snapshot_func, snapshot_func_args,
                                          verbose);
}
void finalize_nhdp_distributions(NanoporeHDP* nhdp) {
    finalize_distributions(nhdp->hdp);
}

void normal_inverse_gamma_params_from_minION(const char* model_filepath, double* mu_out, double* nu_out,
                                             double* alpha_out, double* beta_out) {
    // model format:
    // stateNumber \t alphabetSize \t alphabet \t kmerSize
    // [level_mean, level_stdv, noise_mean, noise_stdv, noise_lambda]
    FILE* model_file = fopen(model_filepath, "r");
    
    char* line = stFile_getLineFromFile(model_file);
    stList* tokens = stString_split(line);

    if (stList_length(tokens) != 4) {
        st_errAbort("normal_inverse_gamma_params_from_minION: Model format has changed invalid model"
                            "found here %s\n", model_filepath);
    }
    free(line);
    stList_destruct(tokens);

    // ignore transitions line
    line = stFile_getLineFromFile(model_file);
    tokens = stString_split(line);
    if (stList_length(tokens) != 10) {
        st_errnoAbort("More than 3-state hmm transitions parameters found\n");
    }

    line = stFile_getLineFromFile(model_file);
    tokens = stString_split(line);
    
    int64_t table_length = (stList_length(tokens) - MODEL_ROW_HEADER_LENGTH) / MODEL_ENTRY_LENGTH;
    double* means = (double*) malloc(sizeof(double) * table_length);
    double* precisions = (double*) malloc(sizeof(double) * table_length);
    
    int64_t mean_offset = MODEL_ROW_HEADER_LENGTH + MODEL_MEAN_ENTRY;  // 1
    int64_t noise_offset = MODEL_ROW_HEADER_LENGTH + MODEL_NOISE_ENTRY;  // 2
    char* mean_str;
    char* noise_str;
    double noise;
    for (int i = 0; i < table_length; i++) {
        mean_str = (char*) stList_get(tokens, mean_offset + i * MODEL_ENTRY_LENGTH);
        sscanf(mean_str, "%lf", &(means[i]));
        
        noise_str = (char*) stList_get(tokens, noise_offset + i * MODEL_ENTRY_LENGTH);
        sscanf(noise_str, "%lf", &noise);
        precisions[i] = 1.0 / (noise * noise);
    }
    
    free(line);
    stList_destruct(tokens);
    mle_normal_inverse_gamma_params(means, precisions, table_length, mu_out, nu_out, alpha_out, beta_out);
    free(means);
    free(precisions);
    
    fclose(model_file);
}

// fixed concentration parameters 'gamma' for each depth
HierarchicalDirichletProcess* minION_hdp(int64_t num_dps, int64_t depth, double* gamma, double sampling_grid_start,
                                         double sampling_grid_stop, int64_t sampling_grid_length,
                                         const char* model_filepath) {
    
    double mu, nu, alpha, beta;
    normal_inverse_gamma_params_from_minION(model_filepath, &mu, &nu, &alpha, &beta);
    return new_hier_dir_proc(num_dps, depth, gamma, sampling_grid_start, sampling_grid_stop,
                             sampling_grid_length, mu, nu, alpha, beta);
}

// Gamma distribution prior on the concentration parameters 'gamma'
// must designate vector of 'alpha' and 'beta' parameters of distribution for each depth
HierarchicalDirichletProcess* minION_hdp_2(int64_t num_dps, int64_t depth, double* gamma_alpha,
                                           double* gamma_beta, double sampling_grid_start,
                                           double sampling_grid_stop, int64_t sampling_grid_length,
                                           const char* model_filepath) {
    
    double mu, nu, alpha, beta;
    normal_inverse_gamma_params_from_minION(model_filepath, &mu, &nu, &alpha, &beta);
    return new_hier_dir_proc_2(num_dps, depth, gamma_alpha, gamma_beta, sampling_grid_start,
                               sampling_grid_stop, sampling_grid_length, mu, nu, alpha, beta);
}

void update_nhdp_from_alignment(NanoporeHDP* nhdp, const char* alignment_filepath, bool has_header) {
    update_nhdp_from_alignment_with_filter(nhdp, alignment_filepath, has_header, NULL);
}

void update_nhdp_from_alignment_with_filter(NanoporeHDP* nhdp, const char* alignment_filepath,
                                            bool has_header, const char* strand_filter) {
    
    stList* signal_list = stList_construct3(0, &free);
    stList* dp_id_list = stList_construct3(0, &free);
    
    FILE* align_file = fopen(alignment_filepath, "r");
    if (align_file == NULL) {
        fprintf(stderr, "Alignment %s file does not exist.\n", alignment_filepath);
        exit(EXIT_FAILURE);
    }
    
    stList* tokens;
    int64_t line_length;
    char* kmer;
    char* strand;
    char* signal_str;
    int64_t* dp_id_ptr;
    double* signal_ptr;
    bool warned = false;
    int proceed = 0;
    
    char* line = stFile_getLineFromFile(align_file);
    if (has_header) {
        line = stFile_getLineFromFile(align_file);
    }
    while (line != NULL) {
        tokens = stString_split(line);
        line_length = stList_length(tokens);
        
        if (!warned) {
            if ((line_length != NUM_ALIGNMENT_COLS) || (line_length != NUM_ASSIGNMENT_COLS)) {
                fprintf(stderr, "Input format has changed from design period, HDP may receive incorrect data.\n");
                warned = true;
                continue;
            }
        }
        
        bool using_alignment;
        if (line_length == NUM_ALIGNMENT_COLS) {
            using_alignment = true;
        } else {
            using_alignment = false;
        }

        int strand_col = using_alignment ? ALIGNMENT_STRAND_COL : ASSIGNMENT_STRAND_COL;
        int signal_col = using_alignment ? ALIGNMENT_SIGNAL_COL : ASSIGNMENT_SIGNAL_COL;
        int kmer_col   = using_alignment ? ALIGNMENT_KMER_COL : ASSIGNMENT_KMER_COL;
        
        strand = (char*) stList_get(tokens, strand_col);
        
        if (strand_filter != NULL) {
            proceed = strcmp(strand, strand_filter);
        }
        
        if (proceed == 0) {
            signal_str = (char*) stList_get(tokens, signal_col);
            kmer = (char*) stList_get(tokens, kmer_col);
            
            signal_ptr = (double*) malloc(sizeof(double));
            dp_id_ptr = (int64_t*) malloc(sizeof(int64_t));
            sscanf(signal_str, "%lf", signal_ptr);
            *dp_id_ptr = kmer_id(kmer, nhdp->alphabet, nhdp->alphabet_size, nhdp->kmer_length);
            
            stList_append(signal_list, signal_ptr);
            stList_append(dp_id_list, dp_id_ptr);
        }
        
        stList_destruct(tokens);
        free(line);
        line = stFile_getLineFromFile(align_file);
    }
    
    fclose(align_file);
    
    int64_t data_length;
    
    double* signal = stList_toDoublePtr(signal_list, &data_length);
    int64_t* dp_ids = stList_toIntPtr(dp_id_list, &data_length);
    
    stList_destruct(signal_list);
    stList_destruct(dp_id_list);
    
    reset_hdp_data(nhdp->hdp);
    pass_data_to_hdp(nhdp->hdp, signal, dp_ids, data_length);
}



// n^k
int64_t power(int64_t n, int64_t k) {
    int64_t num = 1;
    
    for (int64_t i = 0; i < k; i++) {
        num *= n;
    }
    
    return num;
}

//  ((n k))
int64_t multiset_number(int64_t n, int64_t k) {
    int64_t num = 1;
    for (int64_t m = n + k - 1; m >= n; m--) {
        num *= m;
    }
    for (int64_t m = k; m >= 2; m--) {
        num /= m;
    }
    return num;
}

int64_t* get_word(int64_t word_id, int64_t alphabet_size, int64_t word_length) {
    int64_t* word = (int64_t*) malloc(sizeof(int64_t) * word_length);
    int64_t id_remainder = word_id;
    for (int64_t i = 0; i < word_length; i++) {
        word[word_length - i - 1] = id_remainder % alphabet_size;
        id_remainder /= alphabet_size;
    }
    return word;
}

int64_t* get_word_multiset(int64_t word_id, int64_t alphabet_size, int64_t word_length) {
    int64_t* multiset = get_word(word_id, alphabet_size, word_length);
    
    // selection sort 'cause whatever
    int64_t min_idx;
    int64_t temp;
    for (int64_t i = 0; i < word_length; i++) {
        min_idx = i;
        for (int64_t j = i + 1; j < word_length; j++) {
            if (multiset[j] < multiset[min_idx]) {
                min_idx = j;
            }
        }
        temp = multiset[i];
        multiset[i] = multiset[min_idx];
        multiset[min_idx] = temp;
    }
    
    return multiset;
}

int64_t multiset_id_internal(int64_t* tail, int64_t tail_length, int64_t alphabet_min, int64_t alphabet_size) {
    int64_t head = tail[0];
    if (tail_length == 1) {
        return head - alphabet_min;
    }
    int64_t step = 0;
    for (int64_t i = alphabet_min; i < alphabet_size; i++) {
        if (head > i) {
            step += multiset_number(alphabet_size - i, tail_length - 1);
        }
        else {
            return step + multiset_id_internal(&(tail[1]), tail_length - 1, i, alphabet_size);
        }
    }
    fprintf(stderr, "Character outside alphabet included in multiset\n");
    exit(EXIT_FAILURE);
}

int64_t multiset_id(int64_t* multiset, int64_t length, int64_t alphabet_size) {
    return multiset_id_internal(multiset, length, 0, alphabet_size);
}

int64_t word_id_to_multiset_id(int64_t word_id, int64_t alphabet_size, int64_t word_length) {
    int64_t* multiset = get_word_multiset(word_id, alphabet_size, word_length);
    int64_t id = multiset_id(multiset, word_length, alphabet_size);
    free(multiset);
    return id;
}

int64_t word_id(int64_t* word, int64_t alphabet_size, int64_t word_length) {
    int64_t id = 0;
    int64_t step = 1;
    for (int64_t i = word_length - 1; i >= 0; i--) {
        id += step * word[i];
        step *= alphabet_size;
    }
    return id;
}

int64_t* kmer_to_word(char* kmer, char* alphabet, int64_t alphabet_size, int64_t kmer_length) {
    int64_t* word = (int64_t*) malloc(sizeof(int64_t) * kmer_length);
    for (int64_t i = 0; i < kmer_length; i++) {
        int64_t j = 0;
        while (kmer[i] != alphabet[j]) {
            j++;
            if (j == alphabet_size) {
                fprintf(stderr, "[signalAlign] - ERROR: K-mer contains character outside alphabet. "
                        "Got offending kmer is: %s. alphabet is %s kmer length %"PRId64"\n",
                        kmer, alphabet, kmer_length);
                exit(EXIT_FAILURE);
            }
        }
        word[i] = j;
    }
    return word;
}

int64_t kmer_id(char* kmer, char* alphabet, int64_t alphabet_size, int64_t kmer_length) {
    int64_t* word = kmer_to_word(kmer, alphabet, alphabet_size, kmer_length);
    int64_t id = word_id(word, alphabet_size, kmer_length);
    free(word);
    return id;
}

int64_t standard_kmer_id(char* kmer, int64_t kmer_length) {
    return kmer_id(kmer, "ACGT", 4, kmer_length);
}

int64_t nhdp_kmer_id(NanoporeHDP* nhdp, char* kmer) {
    return kmer_id(kmer, nhdp->alphabet, nhdp->alphabet_size, nhdp->kmer_length);
}

double get_nanopore_kmer_density(NanoporeHDP* nhdp, void *kmer, void *x) {
    if (kmer == NULL) {
        return LOG_ZERO;
    } else {
        double u = *(double *)x;
        //return dir_proc_density(nhdp->hdp, *(double *) x, nhdp_kmer_id(nhdp, (char *)kmer));
        return dir_proc_density(nhdp->hdp, u, nhdp_kmer_id(nhdp, (char *)kmer));
    }

}

double get_kmer_distr_distance(NanoporeDistributionMetricMemo* memo, char* kmer_1, char* kmer_2) {
    NanoporeHDP* nhdp = memo->nhdp;
    return get_dir_proc_distance(memo->memo, nhdp_kmer_id(nhdp, kmer_1), nhdp_kmer_id(nhdp, kmer_2));
}

NanoporeDistributionMetricMemo* package_nanopore_metric_memo(NanoporeHDP* nhdp, DistributionMetricMemo* memo) {
    NanoporeDistributionMetricMemo* nanopore_memo = (NanoporeDistributionMetricMemo*) malloc(sizeof(NanoporeDistributionMetricMemo));
    nanopore_memo->nhdp = nhdp;
    nanopore_memo->memo = memo;
    return nanopore_memo;
}

NanoporeDistributionMetricMemo* new_nhdp_kl_divergence_memo(NanoporeHDP* nhdp) {
    return package_nanopore_metric_memo(nhdp, new_kl_divergence_memo(nhdp->hdp));
}

NanoporeDistributionMetricMemo* new_nhdp_hellinger_distance_memo(NanoporeHDP* nhdp) {
    return package_nanopore_metric_memo(nhdp, new_hellinger_distance_memo(nhdp->hdp));
}

NanoporeDistributionMetricMemo* new_nhdp_l2_distance_memo(NanoporeHDP* nhdp) {
    return package_nanopore_metric_memo(nhdp, new_l2_distance_memo(nhdp->hdp));
}

NanoporeDistributionMetricMemo* new_nhdp_shannon_jensen_distance_memo(NanoporeHDP* nhdp) {
    return package_nanopore_metric_memo(nhdp, new_shannon_jensen_distance_memo(nhdp->hdp));
}

double compare_nhdp_distrs_kl_divergence(NanoporeHDP* nhdp_1, char* kmer_1,
                                         NanoporeHDP* nhdp_2, char* kmer_2) {
    return compare_hdp_distrs_kl_divergence(nhdp_1->hdp, nhdp_kmer_id(nhdp_1, kmer_1),
                                            nhdp_2->hdp, nhdp_kmer_id(nhdp_2, kmer_2));
}

double compare_nhdp_distrs_l2_distance(NanoporeHDP* nhdp_1, char* kmer_1,
                                       NanoporeHDP* nhdp_2, char* kmer_2) {
    return compare_hdp_distrs_l2_distance(nhdp_1->hdp, nhdp_kmer_id(nhdp_1, kmer_1),
                                          nhdp_2->hdp, nhdp_kmer_id(nhdp_2, kmer_2));
}


double compare_nhdp_distrs_shannon_jensen_distance(NanoporeHDP* nhdp_1, char* kmer_1,
                                                   NanoporeHDP* nhdp_2, char* kmer_2) {
    return compare_hdp_distrs_shannon_jensen_distance(nhdp_1->hdp, nhdp_kmer_id(nhdp_1, kmer_1),
                                                      nhdp_2->hdp, nhdp_kmer_id(nhdp_2, kmer_2));
}


double compare_nhdp_distrs_hellinger_distance(NanoporeHDP* nhdp_1, char* kmer_1,
                                              NanoporeHDP* nhdp_2, char* kmer_2) {
    return compare_hdp_distrs_hellinger_distance(nhdp_1->hdp, nhdp_kmer_id(nhdp_1, kmer_1),
                                                 nhdp_2->hdp, nhdp_kmer_id(nhdp_2, kmer_2));
}

double kmer_distr_expected_val(NanoporeHDP* nhdp, char* kmer) {
    return dir_proc_expected_val(nhdp->hdp, nhdp_kmer_id(nhdp, kmer));
}

double kmer_distr_variance(NanoporeHDP* nhdp, char* kmer) {
    return dir_proc_variance(nhdp->hdp, nhdp_kmer_id(nhdp, kmer));
}

int64_t flat_hdp_num_dps(int64_t alphabet_size, int64_t kmer_length) {
    int64_t num_leaves = power(alphabet_size, kmer_length);
    return num_leaves + 1;
}

void flat_hdp_model_internal(HierarchicalDirichletProcess* hdp, int64_t alphabet_size, int64_t kmer_length) {
    int64_t last_dp_id = power(alphabet_size, kmer_length);
    
    for (int64_t id = 0; id < last_dp_id; id++) {
        set_dir_proc_parent(hdp, id, last_dp_id);
    }
}

NanoporeHDP* flat_hdp_model(const char* alphabet, int64_t alphabet_size, int64_t kmer_length,
                            double base_gamma, double leaf_gamma,
                            double sampling_grid_start, double sampling_grid_stop, int64_t sampling_grid_length,
                            const char* model_filepath) {
    
    double* gamma_params = (double*) malloc(sizeof(double) * 2);
    gamma_params[0] = base_gamma;
    gamma_params[1] = leaf_gamma;
    
    int64_t num_dps = flat_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp(num_dps, 2, gamma_params, sampling_grid_start,
                                                   sampling_grid_stop, sampling_grid_length,
                                                   model_filepath);
    
    flat_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    finalize_hdp_structure(hdp);
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    return nhdp;
}

NanoporeHDP* flat_hdp_model_2(const char* alphabet, int64_t alphabet_size, int64_t kmer_length,
                              double base_gamma_alpha, double base_gamma_beta,
                              double leaf_gamma_alpha, double leaf_gamma_beta,
                              double sampling_grid_start, double sampling_grid_stop,
                              int64_t sampling_grid_length, const char* model_filepath) {
    
    double* gamma_alpha = (double*) malloc(sizeof(double) * 2);
    gamma_alpha[0] = base_gamma_alpha;
    gamma_alpha[1] = leaf_gamma_alpha;
    
    double* gamma_beta = (double*) malloc(sizeof(double) * 2);
    gamma_beta[0] = base_gamma_beta;
    gamma_beta[1] = leaf_gamma_beta;
    
    int64_t num_dps = flat_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp_2(num_dps, 2, gamma_alpha, gamma_beta, sampling_grid_start,
                                                     sampling_grid_stop, sampling_grid_length,
                                                     model_filepath);
    
    flat_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    finalize_hdp_structure(hdp);
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    return nhdp;
}

int64_t multiset_hdp_num_dps(int64_t alphabet_size, int64_t kmer_length) {
    int64_t num_leaves = power(alphabet_size, kmer_length);
    int64_t num_middle_dps = multiset_number(alphabet_size, kmer_length);
    return num_leaves + num_middle_dps + 1;
}

void multiset_hdp_model_internal(HierarchicalDirichletProcess* hdp, int64_t alphabet_size, int64_t kmer_length) {
    int64_t num_leaves = power(alphabet_size, kmer_length);
    int64_t num_middle_dps = multiset_number(alphabet_size, kmer_length);
    
    // set kmer parents to multisets
    int64_t multiset_id;
    for (int64_t kmer_id = 0; kmer_id < num_leaves; kmer_id++) {
        multiset_id = word_id_to_multiset_id(kmer_id, alphabet_size, kmer_length);
        set_dir_proc_parent(hdp, kmer_id, num_leaves + multiset_id);
    }
    
    // set multiset parents to base dp
    int64_t last_dp_id = num_leaves + num_middle_dps;
    for (int64_t middle_dp_id = num_leaves; middle_dp_id < last_dp_id; middle_dp_id++) {
        set_dir_proc_parent(hdp, middle_dp_id, last_dp_id);
    }
}

NanoporeHDP* multiset_hdp_model(const char* alphabet, int64_t alphabet_size, int64_t kmer_length,
                                double base_gamma, double middle_gamma, double leaf_gamma,
                                double sampling_grid_start, double sampling_grid_stop, int64_t sampling_grid_length,
                                const char* model_filepath) {
    
    double* gamma_params = (double*) malloc(sizeof(double) * 3);
    gamma_params[0] = base_gamma;
    gamma_params[1] = middle_gamma;
    gamma_params[2] = leaf_gamma;
    
    int64_t num_dps = multiset_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp(num_dps, 3, gamma_params, sampling_grid_start,
                                                   sampling_grid_stop, sampling_grid_length,
                                                   model_filepath);
    
    multiset_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    finalize_hdp_structure(hdp);
    
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    return nhdp;
}

NanoporeHDP* multiset_hdp_model_2(const char* alphabet, int64_t alphabet_size, int64_t kmer_length,
                                  double base_gamma_alpha, double base_gamma_beta,
                                  double middle_gamma_alpha, double middle_gamma_beta,
                                  double leaf_gamma_alpha, double leaf_gamma_beta,
                                  double sampling_grid_start, double sampling_grid_stop, int64_t sampling_grid_length,
                                  const char* model_filepath) {
    
    double* gamma_alpha = (double*) malloc(sizeof(double) * 3);
    gamma_alpha[0] = base_gamma_alpha;
    gamma_alpha[1] = middle_gamma_alpha;
    gamma_alpha[2] = leaf_gamma_alpha;
    
    
    double* gamma_beta = (double*) malloc(sizeof(double) * 3);
    gamma_beta[0] = base_gamma_beta;
    gamma_beta[1] = middle_gamma_beta;
    gamma_beta[2] = leaf_gamma_beta;
    
    int64_t num_dps = multiset_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp_2(num_dps, 3, gamma_alpha, gamma_beta, sampling_grid_start,
                                                     sampling_grid_stop, sampling_grid_length,
                                                     model_filepath);
    
    multiset_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    finalize_hdp_structure(hdp);
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    return nhdp;
}

int64_t middle_2_nts_hdp_num_dps(int64_t alphabet_size, int64_t kmer_length) {
    if (kmer_length <= 2) {
        fprintf(stderr, "k-mer is not long enough for middle 2 nucleotides HDP\n");
        exit(EXIT_FAILURE);
    }
    return power(alphabet_size, kmer_length) + power(alphabet_size, 2) + 1;
}

int64_t kmer_id_to_middle_nts_id(int64_t kmer_id, int64_t alphabet_size, int64_t kmer_length) {
    int64_t* kmer = get_word(kmer_id, alphabet_size, kmer_length);
    int64_t id = alphabet_size * kmer[kmer_length / 2 - 1] + kmer[kmer_length / 2];
    free(kmer);
    return id;
}

void middle_2_nts_hdp_model_internal(HierarchicalDirichletProcess* hdp, int64_t alphabet_size, int64_t kmer_length) {
    
    int64_t num_leaves = power(alphabet_size, kmer_length);
    int64_t num_middle_dps = power(alphabet_size, 2);
    
    int64_t middle_dp_id;
    for (int64_t kmer_id = 0; kmer_id < num_leaves; kmer_id++) {
        middle_dp_id = kmer_id_to_middle_nts_id(kmer_id, alphabet_size, kmer_length);
        set_dir_proc_parent(hdp, kmer_id, middle_dp_id + num_leaves);
    }
    
    int64_t last_dp_id = num_leaves + num_middle_dps;
    for (int64_t id = num_leaves; id < last_dp_id; id++) {
        set_dir_proc_parent(hdp, id, last_dp_id);
    }
}

NanoporeHDP* middle_2_nts_hdp_model(const char* alphabet, int64_t alphabet_size, int64_t kmer_length,
                                    double base_gamma, double middle_gamma, double leaf_gamma,
                                    double sampling_grid_start, double sampling_grid_stop, int64_t sampling_grid_length,
                                    const char* model_filepath) {
    if (kmer_length % 2 != 0) {
        fprintf(stderr, "Warning: middle two nucleotides of odd length kmer is ambiguous. Resolving arbitrarily.\n");
    }
    
    double* gamma_params = (double*) malloc(sizeof(double) * 3);
    gamma_params[0] = base_gamma;
    gamma_params[1] = middle_gamma;
    gamma_params[2] = leaf_gamma;
    
    int64_t num_dps = middle_2_nts_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp(num_dps, 3, gamma_params, sampling_grid_start,
                                                   sampling_grid_stop, sampling_grid_length,
                                                   model_filepath);
    
    middle_2_nts_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    finalize_hdp_structure(hdp);
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    return nhdp;
}

int64_t word_id_to_group_multiset_id(int64_t word_id, int64_t* char_groups, int64_t alphabet_size,
                                     int64_t word_length, int64_t num_groups) {
    int64_t* word = get_word(word_id, alphabet_size, word_length);

    for (int64_t i = 0; i < word_length; i++) {
        word[i] = char_groups[word[i]];
    }

    int64_t min_idx;
    int64_t temp;
    for (int64_t i = 0; i < word_length; i++) {
        min_idx = i;
        for (int64_t j = i + 1; j < word_length; j++) {
            if (word[j] < word[min_idx]) {
                min_idx = j;
            }
        }
        temp = word[i];
        word[i] = word[min_idx];
        word[min_idx] = temp;
    }

    int64_t id = multiset_id(word, word_length, num_groups);
    free(word);
    return id;
}

int64_t group_multiset_hdp_num_dps(int64_t alphabet_size, int64_t* char_groups, int64_t kmer_length) {
    int64_t num_groups = 0;
    for (int64_t i = 0; i < alphabet_size; i++) {
        if (char_groups[i] + 1 > num_groups) {
            num_groups = char_groups[i] + 1;
        }
    }

    int64_t num_leaves = power(alphabet_size, kmer_length);
    int64_t num_middle_dps = multiset_number(num_groups, kmer_length);
    return num_leaves + num_middle_dps + 1;
}

void group_multiset_hdp_model_internal(HierarchicalDirichletProcess* hdp, int64_t* char_groups,
                                       int64_t alphabet_size, int64_t kmer_length) {

    int64_t num_groups = 0;
    for (int64_t i = 0; i < alphabet_size; i++) {
        if (char_groups[i] + 1 > num_groups) {
            num_groups = char_groups[i] + 1;
        }
    }

    int64_t num_leaves = power(alphabet_size, kmer_length);
    int64_t num_middle_dps = multiset_number(num_groups, kmer_length);

    // set kmer parents to multisets
    int64_t multiset_id;
    for (int64_t kmer_id = 0; kmer_id < num_leaves; kmer_id++) {
        multiset_id = word_id_to_group_multiset_id(kmer_id, char_groups, alphabet_size, kmer_length, num_groups);
        set_dir_proc_parent(hdp, kmer_id, num_leaves + multiset_id);
    }

    // set multiset parents to base dp
    int64_t last_dp_id = num_leaves + num_middle_dps;
    for (int64_t middle_dp_id = num_leaves; middle_dp_id < last_dp_id; middle_dp_id++) {
        set_dir_proc_parent(hdp, middle_dp_id, last_dp_id);
    }
}

void confirm_valid_groupings(int64_t* char_groups, int64_t alphabet_size) {

    for (int64_t i = 0; i < alphabet_size; i++) {
        if (char_groups[i] < 0) {
            fprintf(stderr, "Group numbers must be non-negative.\n");
            exit(EXIT_FAILURE);
        }
    }

    int64_t num_groups = 0;
    for (int64_t i = 0; i < alphabet_size; i++) {
        if (char_groups[i] + 1 > num_groups) {
            num_groups = char_groups[i] + 1;
        }
    }

    for (int64_t i = 0; i < num_groups; i++) {
        bool found_group = false;
        for (int64_t j = 0; j < alphabet_size; j++) {
            if (char_groups[j] == i) {
                found_group = true;
                break;
            }
        }
        if (!found_group) {
            fprintf(stderr, "Groups must be consecutively numbered starting with 0.\n");
            exit(EXIT_FAILURE);
        }
    }
}

int64_t* alphabet_sort_groups(const char* alphabet, int64_t* char_groups, int64_t alphabet_size) {
    char* aux_alphabet = (char*) malloc(sizeof(char) * alphabet_size);
    int64_t* sorted_char_groups = (int64_t*) malloc(sizeof(int64_t) * alphabet_size);

    for (int64_t i = 0; i < alphabet_size; i++) {
        aux_alphabet[i] = alphabet[i];
        sorted_char_groups[i] = char_groups[i];
    }

    int64_t temp_group;
    char temp_char;
    int64_t min_idx;
    for (int64_t i = 0; i < alphabet_size; i++) {
        min_idx = i;
        for (int64_t j = i + 1; j < alphabet_size; j++) {
            if (aux_alphabet[j] < aux_alphabet[min_idx]) {
                min_idx = j;
            }
        }
        temp_char = aux_alphabet[i];
        aux_alphabet[i] = aux_alphabet[min_idx];
        aux_alphabet[min_idx] = temp_char;

        temp_group = sorted_char_groups[i];
        sorted_char_groups[i] = sorted_char_groups[min_idx];
        sorted_char_groups[min_idx] = temp_group;
    }

    free(aux_alphabet);
    return sorted_char_groups;
}

// assumes char_groups are 0-based and consecutively numbered
NanoporeHDP* group_multiset_hdp_model(const char* alphabet, int64_t* char_groups, int64_t alphabet_size, int64_t kmer_length,
                                      double base_gamma, double middle_gamma, double leaf_gamma,
                                      double sampling_grid_start, double sampling_grid_stop, int64_t sampling_grid_length,
                                      const char* model_filepath) {

    confirm_valid_groupings(char_groups, alphabet_size);

    double* gamma_params = (double*) malloc(sizeof(double) * 3);
    gamma_params[0] = base_gamma;
    gamma_params[1] = middle_gamma;
    gamma_params[2] = leaf_gamma;

    int64_t num_dps = group_multiset_hdp_num_dps(alphabet_size, char_groups, kmer_length);

    HierarchicalDirichletProcess* hdp = minION_hdp(num_dps, 3, gamma_params, sampling_grid_start,
                                                   sampling_grid_stop, sampling_grid_length,
                                                   model_filepath);

    int64_t* sorted_char_groups = alphabet_sort_groups(alphabet, char_groups, alphabet_size);
    group_multiset_hdp_model_internal(hdp, sorted_char_groups, alphabet_size, kmer_length);
    free(sorted_char_groups);

    finalize_hdp_structure(hdp);

    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);

    return nhdp;
}

// assumes char_groups are 0-based and consecutively numbered
NanoporeHDP* group_multiset_hdp_model_2(const char* alphabet, int64_t* char_groups, int64_t alphabet_size, int64_t kmer_length,
                                        double base_gamma_alpha, double base_gamma_beta, double middle_gamma_alpha,
                                        double middle_gamma_beta, double leaf_gamma_alpha, double leaf_gamma_beta,
                                        double sampling_grid_start, double sampling_grid_stop, int64_t sampling_grid_length,
                                        const char* model_filepath) {

    confirm_valid_groupings(char_groups, alphabet_size);

    double *gamma_alpha = (double *) malloc(sizeof(double) * 3);
    gamma_alpha[0] = base_gamma_alpha;
    gamma_alpha[1] = middle_gamma_alpha;
    gamma_alpha[2] = leaf_gamma_alpha;


    double *gamma_beta = (double *) malloc(sizeof(double) * 3);
    gamma_beta[0] = base_gamma_beta;
    gamma_beta[1] = middle_gamma_beta;
    gamma_beta[2] = leaf_gamma_beta;

    int64_t num_dps = group_multiset_hdp_num_dps(alphabet_size, char_groups, kmer_length);

    HierarchicalDirichletProcess *hdp = minION_hdp_2(num_dps, 3, gamma_alpha, gamma_beta, sampling_grid_start,
                                                     sampling_grid_stop, sampling_grid_length,
                                                     model_filepath);

    int64_t *sorted_char_groups = alphabet_sort_groups(alphabet, char_groups, alphabet_size);
    group_multiset_hdp_model_internal(hdp, sorted_char_groups, alphabet_size, kmer_length);
    free(sorted_char_groups);

    finalize_hdp_structure(hdp);

    NanoporeHDP *nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);

    return nhdp;
}

NanoporeHDP* middle_2_nts_hdp_model_2(const char* alphabet, int64_t alphabet_size, int64_t kmer_length,
                                      double base_gamma_alpha, double base_gamma_beta, double middle_gamma_alpha,
                                      double middle_gamma_beta, double leaf_gamma_alpha, double leaf_gamma_beta,
                                      double sampling_grid_start, double sampling_grid_stop,
                                      int64_t sampling_grid_length, const char* model_filepath) {
    if (kmer_length % 2 != 0) {
        fprintf(stderr, "Warning: middle 2 nucleotides of odd length kmer is ambiguous. Resolving arbitrarily.\n");
    }
    
    double* gamma_alpha = (double*) malloc(sizeof(double) * 3);
    gamma_alpha[0] = base_gamma_alpha;
    gamma_alpha[1] = middle_gamma_alpha;
    gamma_alpha[2] = leaf_gamma_alpha;
    
    double* gamma_beta = (double*) malloc(sizeof(double) * 3);
    gamma_beta[0] = base_gamma_beta;
    gamma_beta[1] = middle_gamma_beta;
    gamma_beta[2] = leaf_gamma_beta;
    
    int64_t num_dps = middle_2_nts_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp_2(num_dps, 3, gamma_alpha, gamma_beta, sampling_grid_start,
                                                     sampling_grid_stop, sampling_grid_length,
                                                     model_filepath);
    
    middle_2_nts_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    finalize_hdp_structure(hdp);
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    return nhdp;
}

int64_t purine_composition_hdp_num_dps(int64_t num_purines, int64_t num_pyrimidines, int64_t kmer_length) {
    int64_t num_leaves = power(num_purines + num_pyrimidines, kmer_length);
    int64_t num_middle_dps = kmer_length + 1;
    return num_leaves + num_middle_dps + 1;
}

void purine_composition_hdp_model_internal(HierarchicalDirichletProcess* hdp, bool* purine_alphabet,
                                           int64_t alphabet_size, int64_t kmer_length) {
    int64_t num_leaves = power(alphabet_size, kmer_length);
    int64_t num_middle_dps = kmer_length + 1;
    
    // set kmer parents to purine multisets
    int64_t num_purines;
    int64_t* word;
    for (int64_t kmer_id = 0; kmer_id < num_leaves; kmer_id++) {
        word = get_word(kmer_id, alphabet_size, kmer_length);
        num_purines = 0;
        for (int64_t i = 0; i < kmer_length; i++) {
            if (purine_alphabet[word[i]]) {
                num_purines++;
            }
        }
        free(word);
        set_dir_proc_parent(hdp, kmer_id, num_leaves + num_purines);
    }
    
    // set purine set parents to base dp
    int64_t last_dp_id = num_leaves + num_middle_dps;
    for (int64_t middle_dp_id = num_leaves; middle_dp_id < last_dp_id; middle_dp_id++) {
        set_dir_proc_parent(hdp, middle_dp_id, last_dp_id);
    }
}

NanoporeHDP* purine_composition_hdp_model(char* purine_alphabet, int64_t num_purines,
                                          char* pyrimidine_alphabet, int64_t num_pyrimidines,
                                          int64_t kmer_length, double base_gamma, double middle_gamma,
                                          double leaf_gamma, double sampling_grid_start, double sampling_grid_stop,
                                          int64_t sampling_grid_length, const char* model_filepath) {
    
    double* gamma_params = (double*) malloc(sizeof(double) * 3);
    gamma_params[0] = base_gamma;
    gamma_params[1] = middle_gamma;
    gamma_params[2] = leaf_gamma;
    
    int64_t num_dps = purine_composition_hdp_num_dps(num_purines, num_pyrimidines, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp(num_dps, 3, gamma_params, sampling_grid_start,
                                                   sampling_grid_stop, sampling_grid_length,
                                                   model_filepath);
    
    int64_t alphabet_size = num_purines + num_pyrimidines;
    char* alphabet = (char*) malloc(sizeof(char) * alphabet_size);
    for (int64_t i = 0; i < num_purines; i++) {
        alphabet[i] = purine_alphabet[i];
    }
    for (int64_t i = 0; i < num_pyrimidines; i++) {
        alphabet[i + num_purines] = pyrimidine_alphabet[i];
    }
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    // get back the alphabet in the internal ordering
    free(alphabet);
    alphabet = get_nanopore_hdp_alphabet(nhdp);
    bool* purines = (bool*) malloc(sizeof(bool) * alphabet_size);
    for (int64_t i = 0; i < num_purines; i++) {
        purines[i] = false;
        for (int64_t j = 0; j < num_purines; j++) {
            if (alphabet[i] == purine_alphabet[j]) {
                purines[i] = true;
                break;
            }
        }
    }
    free(alphabet);
    
    purine_composition_hdp_model_internal(hdp, purines, alphabet_size, kmer_length);
    free(purines);
    
    finalize_hdp_structure(hdp);
    
    return nhdp;
}

NanoporeHDP* purine_composition_hdp_model_2(char* purine_alphabet, int64_t num_purines,
                                            char* pyrimidine_alphabet, int64_t num_pyrimidines,
                                            int64_t kmer_length, double base_gamma_alpha, double base_gamma_beta,
                                            double middle_gamma_alpha, double middle_gamma_beta,
                                            double leaf_gamma_alpha, double leaf_gamma_beta, double sampling_grid_start,
                                            double sampling_grid_stop, int64_t sampling_grid_length,
                                            const char* model_filepath) {
    
    double* gamma_alpha = (double*) malloc(sizeof(double) * 3);
    gamma_alpha[0] = base_gamma_alpha;
    gamma_alpha[1] = middle_gamma_alpha;
    gamma_alpha[2] = leaf_gamma_alpha;
    
    
    double* gamma_beta = (double*) malloc(sizeof(double) * 3);
    gamma_beta[0] = base_gamma_beta;
    gamma_beta[1] = middle_gamma_beta;
    gamma_beta[2] = leaf_gamma_beta;
    
    int64_t num_dps = purine_composition_hdp_num_dps(num_purines, num_pyrimidines, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp_2(num_dps, 3, gamma_alpha, gamma_beta, sampling_grid_start,
                                                     sampling_grid_stop, sampling_grid_length,
                                                     model_filepath);
    
    int64_t alphabet_size = num_purines + num_pyrimidines;
    char* alphabet = (char*) malloc(sizeof(char) * alphabet_size);
    for (int64_t i = 0; i < num_purines; i++) {
        alphabet[i] = purine_alphabet[i];
    }
    for (int64_t i = 0; i < num_pyrimidines; i++) {
        alphabet[i + num_purines] = pyrimidine_alphabet[i];
    }
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    // get back the alphabet in the internal ordering
    free(alphabet);
    alphabet = get_nanopore_hdp_alphabet(nhdp);
    bool* purines = (bool*) malloc(sizeof(bool) * alphabet_size);
    for (int64_t i = 0; i < alphabet_size; i++) {
        purines[i] = false;
        for (int64_t j = 0; j < num_purines; j++) {
            if (alphabet[i] == purine_alphabet[j]) {
                purines[i] = true;
                break;
            }
        }
    }
    
    free(alphabet);
    
    purine_composition_hdp_model_internal(hdp, purines, alphabet_size, kmer_length);
    free(purines);
    
    finalize_hdp_structure(hdp);
    
    return nhdp;
}

void serialize_nhdp(NanoporeHDP* nhdp, const char* filepath) {
    FILE* out = fopen(filepath, "w");
    
    fprintf(out, "%"PRId64"\n", nhdp->alphabet_size);
    fprintf(out, "%s\n", nhdp->alphabet);
    fprintf(out, "%"PRId64"\n", nhdp->kmer_length);
    serialize_hdp(nhdp->hdp, out);
    
    fclose(out);
}

NanoporeHDP* deserialize_nhdp(const char* filepath) {
    FILE* in = fopen(filepath, "r");
    
    char* line = stFile_getLineFromFile(in);
    int64_t alphabet_size;
    sscanf(line, "%"SCNd64, &alphabet_size);
    free(line);
    
    line = stFile_getLineFromFile(in);
    char* alphabet = (char*) malloc(sizeof(char) * alphabet_size);
    sscanf(line, "%s", alphabet);
    free(line);
    
    line = stFile_getLineFromFile(in);
    int64_t kmer_length;
    sscanf(line, "%"SCNd64, &kmer_length);
    free(line);
    
    HierarchicalDirichletProcess* hdp = deserialize_hdp(in);
    
    fclose(in);
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    free(alphabet);
    
    return nhdp;
}

static void nanoporeHdp_checkThreeLevelPriorParameters(double baseGammaAlpha, double baseGammaBeta,
                                                       double middleGammaAlpha, double middleGammaBeta,
                                                       double leafGammaAlpha, double leafGammaBeta) {
    if ((baseGammaAlpha == NULL_HYPERPARAMETER) || (baseGammaBeta == NULL_HYPERPARAMETER) ||
        (middleGammaAlpha == NULL_HYPERPARAMETER) || (middleGammaBeta == NULL_HYPERPARAMETER) ||
        (leafGammaAlpha == NULL_HYPERPARAMETER) || (leafGammaBeta == NULL_HYPERPARAMETER)) {
        st_errAbort("loadNanoporeHdpFromScratch: You need to provide a alphas and betas for the base, middle, "
                            "and the leaf distributions for the prior for this NanoporeHdp");
    }
}

static void nanoporeHdp_checkThreeLevelFixedParameters(double baseGamma, double middleGamma, double leafGamma) {
    if ((baseGamma == NULL_HYPERPARAMETER) || (leafGamma == NULL_HYPERPARAMETER) ||
        (middleGamma == NULL_HYPERPARAMETER)) {
        st_errAbort("loadNanoporeHdpFromScratch: You need to provide a base gamma, middle gamma, and leaf gamma "
                            "for this NanoporeHdpType\n");
    }
}

static void nanoporeHdp_checkTwoLevelPriorParameters(double baseGammaAlpha, double baseGammaBeta,
                                                     double leafGammaAlpha, double leafGammaBeta) {
    if ((baseGammaAlpha == NULL_HYPERPARAMETER) || (baseGammaBeta == NULL_HYPERPARAMETER) ||
        (leafGammaAlpha == NULL_HYPERPARAMETER) || (leafGammaBeta == NULL_HYPERPARAMETER)) {
        st_errAbort("loadNanoporeHdpFromScratch: You need to provide a alphas and betas for the base and the leaf"
                            "distributions for the prior for this NanoporeHdp");
    }
}


static NanoporeHDP *loadNanoporeHdpFromScratch(NanoporeHdpType nHdpType, const char *modelFile, int64_t kmerLength,
                                               double baseGamma, double middleGamma, double leafGamma,
                                               double baseGammaAlpha, double baseGammaBeta,
                                               double middleGammaAlpha, double middleGammaBeta,
                                               double leafGammaAlpha, double leafGammaBeta,
                                               double samplingGridStart, double samplingGridEnd,
                                               int64_t samplingGridLength) {
    if (nHdpType == singleLevelFixed) {
        if ((baseGamma == NULL_HYPERPARAMETER) || (leafGamma == NULL_HYPERPARAMETER)) {
            st_errAbort("loadNanoporeHdpFromScratch: You need to provide a base gamma and leaf gamma "
                                "for this NanoporeHdpType\n");
        }

        NanoporeHDP *nHdp = flat_hdp_model(METHYL_HYDROXY_CYTOSINE_ALPHA, SYMBOL_NUMBER_EPIGENETIC_C, kmerLength,
                                           baseGamma, leafGamma,
                                           samplingGridStart, samplingGridEnd, samplingGridLength, modelFile);

        return nHdp;
    }
    if (nHdpType == singleLevelPrior) {
        nanoporeHdp_checkTwoLevelPriorParameters(baseGammaAlpha, baseGammaBeta, leafGammaAlpha, leafGammaBeta);

        NanoporeHDP *nHdp = flat_hdp_model_2(METHYL_HYDROXY_CYTOSINE_ALPHA, SYMBOL_NUMBER_EPIGENETIC_C, kmerLength,
                                             baseGammaAlpha, baseGammaBeta, leafGammaAlpha, leafGammaBeta,
                                             samplingGridStart, samplingGridEnd, samplingGridLength,
                                             modelFile);
        return nHdp;
    }
    if (nHdpType == singleLevelPrior2) {
        nanoporeHdp_checkTwoLevelPriorParameters(baseGammaAlpha, baseGammaBeta, leafGammaAlpha, leafGammaBeta);
        NanoporeHDP *nHdp = flat_hdp_model_2(METHYL_CYTOSINE_ALPHA, SYMBOL_NUMBER, kmerLength,
                                             baseGammaAlpha, baseGammaBeta, leafGammaAlpha, leafGammaBeta,
                                             samplingGridStart, samplingGridEnd, samplingGridLength,
                                             modelFile);
        return nHdp;
    }
    if (nHdpType == singleLevelPriorEcoli) {
        nanoporeHdp_checkTwoLevelPriorParameters(baseGammaAlpha, baseGammaBeta, leafGammaAlpha, leafGammaBeta);
        NanoporeHDP *nHdp = flat_hdp_model_2(METHYL_CYTOSINE_ADENOSINE_ALPHA, SYMBOL_NUMBER_METHYL_CA, kmerLength,
                                             baseGammaAlpha, baseGammaBeta, leafGammaAlpha, leafGammaBeta,
                                             samplingGridStart, samplingGridEnd, samplingGridLength,
                                             modelFile);
        return nHdp;
    }
    if (nHdpType == multisetFixed) {
        nanoporeHdp_checkThreeLevelFixedParameters(baseGamma, middleGamma, leafGamma);

        NanoporeHDP *nHdp = multiset_hdp_model(METHYL_HYDROXY_CYTOSINE_ALPHA, SYMBOL_NUMBER_EPIGENETIC_C, kmerLength,
                                               baseGamma, middleGamma, leafGamma,
                                               samplingGridStart, samplingGridEnd, samplingGridLength,
                                               modelFile);
        return nHdp;
    }
    if (nHdpType == multisetPrior) {
        nanoporeHdp_checkThreeLevelPriorParameters(baseGammaAlpha, baseGammaBeta,
                                                   middleGammaAlpha, middleGammaBeta,
                                                   leafGammaAlpha, leafGammaBeta);

        NanoporeHDP *nHdp = multiset_hdp_model_2(METHYL_HYDROXY_CYTOSINE_ALPHA, SYMBOL_NUMBER_EPIGENETIC_C, kmerLength,
                                                 baseGammaAlpha, baseGammaBeta,
                                                 middleGammaAlpha, middleGammaBeta,
                                                 leafGammaAlpha, leafGammaBeta,
                                                 samplingGridStart, samplingGridEnd, samplingGridLength,
                                                 modelFile);
        return nHdp;
    }
    if (nHdpType == multisetPrior2) {
        nanoporeHdp_checkThreeLevelPriorParameters(baseGammaAlpha, baseGammaBeta,
                                                   middleGammaAlpha, middleGammaBeta,
                                                   leafGammaAlpha, leafGammaBeta);

        NanoporeHDP *nHdp = multiset_hdp_model_2(METHYL_CYTOSINE_ALPHA, SYMBOL_NUMBER, kmerLength,
                                                 baseGammaAlpha, baseGammaBeta,
                                                 middleGammaAlpha, middleGammaBeta,
                                                 leafGammaAlpha, leafGammaBeta,
                                                 samplingGridStart, samplingGridEnd, samplingGridLength,
                                                 modelFile);
        return nHdp;
    }
    if (nHdpType == multisetPriorEcoli) {
        nanoporeHdp_checkThreeLevelPriorParameters(baseGammaAlpha, baseGammaBeta,
                                                   middleGammaAlpha, middleGammaBeta,
                                                   leafGammaAlpha, leafGammaBeta);
        NanoporeHDP *nHdp = multiset_hdp_model_2(METHYL_CYTOSINE_ADENOSINE_ALPHA, SYMBOL_NUMBER_METHYL_CA, kmerLength,
                                                 baseGammaAlpha, baseGammaBeta,
                                                 middleGammaAlpha, middleGammaBeta,
                                                 leafGammaAlpha, leafGammaBeta,
                                                 samplingGridStart, samplingGridEnd, samplingGridLength,
                                                 modelFile);
        return nHdp;
    }
    if (nHdpType == compFixed) {
        nanoporeHdp_checkThreeLevelFixedParameters(baseGamma, middleGamma, leafGamma);

        NanoporeHDP *nHdp = purine_composition_hdp_model(PURINES, 2, PYRIMIDINES, 4, kmerLength,
                                                         baseGamma, middleGamma, leafGamma,
                                                         samplingGridStart, samplingGridEnd,
                                                         samplingGridLength, modelFile);
        return nHdp;
    }
    if (nHdpType == compPrior) {
        nanoporeHdp_checkThreeLevelPriorParameters(baseGammaAlpha, baseGammaBeta, middleGammaAlpha,
                                                   middleGammaBeta, leafGammaAlpha, leafGammaBeta);

        NanoporeHDP *nHdp = purine_composition_hdp_model_2(PURINES, 2, PYRIMIDINES, 4, kmerLength,
                                                           baseGammaAlpha, baseGammaBeta,
                                                           middleGammaAlpha, middleGammaBeta,
                                                           leafGammaAlpha, leafGammaBeta,
                                                           samplingGridStart, samplingGridEnd,
                                                           samplingGridLength, modelFile);
        return nHdp;

    }
    if (nHdpType == middleNtsFixed) {
        nanoporeHdp_checkThreeLevelFixedParameters(baseGamma, middleGamma, leafGamma);

        NanoporeHDP *nHdp = middle_2_nts_hdp_model(METHYL_HYDROXY_CYTOSINE_ALPHA, SYMBOL_NUMBER_EPIGENETIC_C, kmerLength,
                                                   baseGamma, middleGamma, leafGamma,
                                                   samplingGridStart, samplingGridEnd, samplingGridLength,
                                                   modelFile);
        return nHdp;

    }
    if (nHdpType == middleNtsPrior) {
        nanoporeHdp_checkThreeLevelPriorParameters(baseGammaAlpha, baseGammaBeta, middleGammaAlpha,
                                                   middleGammaBeta, leafGammaAlpha, leafGammaBeta);

        NanoporeHDP *nHdp = middle_2_nts_hdp_model_2(METHYL_HYDROXY_CYTOSINE_ALPHA, SYMBOL_NUMBER_EPIGENETIC_C, kmerLength,
                                                     baseGammaAlpha, baseGammaBeta,
                                                     middleGammaAlpha, middleGammaBeta,
                                                     leafGammaAlpha, leafGammaBeta,
                                                     samplingGridStart, samplingGridEnd, samplingGridLength,
                                                     modelFile);
        return nHdp;
    }
    if (nHdpType == groupMultisetFixed) {
        nanoporeHdp_checkThreeLevelFixedParameters(baseGamma, middleGamma, leafGamma);
        // ACEGOT
        // {0, 1, 1, 2, 1, 3}
        int64_t groups[6] = {0, 1, 1, 2, 1, 3};

        NanoporeHDP *nHdp = group_multiset_hdp_model(METHYL_HYDROXY_CYTOSINE_ALPHA, groups, SYMBOL_NUMBER_EPIGENETIC_C, kmerLength,
                                                     baseGamma, middleGamma, leafGamma,
                                                     samplingGridStart, samplingGridEnd, samplingGridLength,
                                                     modelFile);
        return nHdp;
    }
    if (nHdpType == groupMultisetPrior) {
        nanoporeHdp_checkThreeLevelPriorParameters(baseGammaAlpha, baseGammaBeta, middleGammaAlpha,
                                                   middleGammaBeta, leafGammaAlpha, leafGammaBeta);
        // ACEGOT
        // {0, 1, 1, 2, 1, 3}
        int64_t groups[6] = {0, 1, 1, 2, 1, 3};

        NanoporeHDP *nHdp = group_multiset_hdp_model_2(METHYL_HYDROXY_CYTOSINE_ALPHA, groups,
                                                       SYMBOL_NUMBER_EPIGENETIC_C, kmerLength,
                                                       baseGammaAlpha, baseGammaBeta,
                                                       middleGammaAlpha, middleGammaBeta,
                                                       leafGammaAlpha, leafGammaBeta,
                                                       samplingGridStart, samplingGridEnd, samplingGridLength,
                                                       modelFile);
        return nHdp;
    }
    else {
        fprintf(stderr, "loadNanoporeHdpFromScratch: - error making HDP from scratch\n");
        exit(EXIT_FAILURE);
    }
}

void nanoporeHdp_buildNanoporeHdpFromAlignment(NanoporeHdpType type, int64_t kmerLength,
                                               const char *templateModelFile, const char* complementModelFile,
                                               const char *alignments,
                                               const char *templateHDP, const char *complementHDP,
                                               int64_t nbSamples, int64_t burnIn, int64_t thinning, bool verbose,
                                               double baseGamma, double middleGamma, double leafGamma,
                                               double baseGammaAlpha, double baseGammaBeta,
                                               double middleGammaAlpha, double middleGammaBeta,
                                               double leafGammaAlpha, double leafGammaBeta,
                                               double samplingGridStart, double samplingGridEnd,
                                               int64_t samplingGridLength) {
    fprintf(stderr, "Building Nanopore HDP\n");
#pragma omp parallel sections
 {
    {
        fprintf(stderr, "Updating Template HDP from alignments...\n");
        NanoporeHDP *nHdpT = loadNanoporeHdpFromScratch(type, templateModelFile, kmerLength,
                                                        baseGamma, middleGamma, leafGamma,
                                                        baseGammaAlpha, baseGammaBeta,
                                                        middleGammaAlpha, middleGammaBeta,
                                                        leafGammaAlpha, leafGammaBeta,
                                                        samplingGridStart, samplingGridEnd, samplingGridLength);
        update_nhdp_from_alignment_with_filter(nHdpT, alignments, FALSE, "t");

        fprintf(stderr, "Running Gibbs for template doing %"PRId64"samples, %"PRId64"burn in, %"PRId64"thinning.\n",
                nbSamples, burnIn, thinning);

        execute_nhdp_gibbs_sampling(nHdpT, nbSamples, burnIn, thinning, verbose);
        finalize_nhdp_distributions(nHdpT);

        fprintf(stderr, "Serializing template to %s...\n", templateHDP);
        serialize_nhdp(nHdpT, templateHDP);
        destroy_nanopore_hdp(nHdpT);
    }
#pragma omp section
    {
        fprintf(stderr, "Updating Complement HDP from alignments...\n");
        NanoporeHDP *nHdpC = loadNanoporeHdpFromScratch(type, complementModelFile, kmerLength,
                                                        baseGamma, middleGamma, leafGamma,
                                                        baseGammaAlpha, baseGammaBeta,
                                                        middleGammaAlpha, middleGammaBeta,
                                                        leafGammaAlpha, leafGammaBeta,
                                                        samplingGridStart, samplingGridEnd, samplingGridLength);
        update_nhdp_from_alignment_with_filter(nHdpC, alignments, FALSE, "c");

        fprintf(stderr, "Running Gibbs for complement doing %"PRId64"samples, %"PRId64"burn in, %"PRId64"thinning.\n",
                nbSamples, burnIn, thinning);
        execute_nhdp_gibbs_sampling(nHdpC, nbSamples, burnIn, thinning, verbose);
        finalize_nhdp_distributions(nHdpC);

        fprintf(stderr, "Serializing complement to %s...\n", complementHDP);
        serialize_nhdp(nHdpC, complementHDP);
        destroy_nanopore_hdp(nHdpC);
    }
}
}

void nanoporeHdp_buildOneDHdpFromAlignment(NanoporeHdpType type, int64_t kmerLength,
                                           const char *templateModelFile,
                                           const char *alignments,
                                           const char *templateHDP,
                                           int64_t nbSamples, int64_t burnIn, int64_t thinning, bool verbose,
                                           double baseGamma, double middleGamma, double leafGamma,
                                           double baseGammaAlpha, double baseGammaBeta,
                                           double middleGammaAlpha, double middleGammaBeta,
                                           double leafGammaAlpha, double leafGammaBeta,
                                           double samplingGridStart, double samplingGridEnd,
                                           int64_t samplingGridLength) {
    fprintf(stderr, "Updating Template HDP from alignments...\n");
    NanoporeHDP *nHdpT = loadNanoporeHdpFromScratch(type, templateModelFile, kmerLength,
                                                    baseGamma, middleGamma, leafGamma,
                                                    baseGammaAlpha, baseGammaBeta,
                                                    middleGammaAlpha, middleGammaBeta,
                                                    leafGammaAlpha, leafGammaBeta,
                                                    samplingGridStart, samplingGridEnd, samplingGridLength);
    update_nhdp_from_alignment_with_filter(nHdpT, alignments, FALSE, "t");

    fprintf(stderr, "Running Gibbs for template doing %"PRId64"samples, %"PRId64"burn in, %"PRId64"thinning.\n",
            nbSamples, burnIn, thinning);

    execute_nhdp_gibbs_sampling(nHdpT, nbSamples, burnIn, thinning, verbose);
    finalize_nhdp_distributions(nHdpT);

    fprintf(stderr, "Serializing template to %s...\n", templateHDP);
    serialize_nhdp(nHdpT, templateHDP);
    destroy_nanopore_hdp(nHdpT);

}


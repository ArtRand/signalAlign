#ifndef HDP_H_INCLUDED
#define HDP_H_INCLUDED

#include <inttypes.h>
#include <stdbool.h>

typedef struct HierarchicalDirichletProcess HierarchicalDirichletProcess;
typedef struct DistributionMetricMemo DistributionMetricMemo;

// constructors and destructor

HierarchicalDirichletProcess* new_hier_dir_proc(int64_t num_dps, int64_t depth, double* gamma, double sampling_grid_start,
                                                double sampling_grid_stop, int64_t sampling_grid_length, double mu,
                                                double nu, double alpha, double beta);

HierarchicalDirichletProcess* new_hier_dir_proc_2(int64_t num_dps, int64_t depth, double* gamma_alpha, double* gamma_beta,
                                                  double sampling_grid_start, double sampling_grid_stop,
                                                  int64_t sampling_grid_length, double mu, double nu, double alpha,
                                                  double beta);

void destroy_hier_dir_proc(HierarchicalDirichletProcess* hdp);

// topology

void set_dir_proc_parent(HierarchicalDirichletProcess* hdp, int64_t child_id, int64_t parent_id);

void finalize_hdp_structure(HierarchicalDirichletProcess* hdp);

// data management

void pass_data_to_hdp(HierarchicalDirichletProcess* hdp, double* data, int64_t* dp_id, int64_t length);

void reset_hdp_data(HierarchicalDirichletProcess* hdp);

// Gibbs sampling

void execute_gibbs_sampling(HierarchicalDirichletProcess* hdp, int64_t num_samples, int64_t burn_in,
                            int64_t thinning, bool verbose);

void execute_gibbs_sampling_with_snapshots(HierarchicalDirichletProcess* hdp,
                                           int64_t num_samples, int64_t burn_in, int64_t thinning,
                                           void (*snapshot_func)(HierarchicalDirichletProcess*, void*),
                                           void* snapshot_func_args, bool verbose);

void finalize_distributions(HierarchicalDirichletProcess* hdp);

// querying the HDP

double dir_proc_density(HierarchicalDirichletProcess* hdp, double x, int64_t dp_id);

void take_snapshot(HierarchicalDirichletProcess* hdp, int64_t** num_dp_fctrs_out, int64_t* num_dps_out,
                   double** gamma_params_out, int64_t* num_gamma_params_out, double* log_likelihood_out,
                   double* log_density_out);

// get methods

bool is_structure_finalized(HierarchicalDirichletProcess* hdp);
bool is_gamma_random(HierarchicalDirichletProcess* hdp);
bool is_sampling_finalized(HierarchicalDirichletProcess* hdp);
int64_t get_num_dir_proc(HierarchicalDirichletProcess* hdp);
int64_t get_depth(HierarchicalDirichletProcess* hdp);
int64_t get_num_data(HierarchicalDirichletProcess* hdp);
double* get_data_copy(HierarchicalDirichletProcess* hdp);
int64_t* get_data_pt_dp_ids_copy(HierarchicalDirichletProcess* hdp);
double* get_gamma_params_copy(HierarchicalDirichletProcess* hdp);
double get_mu(HierarchicalDirichletProcess* hdp);
double get_nu(HierarchicalDirichletProcess* hdp);
double get_alpha(HierarchicalDirichletProcess* hdp);
double get_beta(HierarchicalDirichletProcess* hdp);
int64_t get_grid_length(HierarchicalDirichletProcess* hdp);
double* get_sampling_grid_copy(HierarchicalDirichletProcess* hdp);
double* get_gamma_alpha_params_copy(HierarchicalDirichletProcess* hdp);
double* get_gamma_beta_params_copy(HierarchicalDirichletProcess* hdp);
int64_t get_dir_proc_num_factors(HierarchicalDirichletProcess* hdp, int64_t dp_id);
int64_t get_dir_proc_parent_id(HierarchicalDirichletProcess* hdp, int64_t dp_id);

// computing distance between DP distributions

double get_dir_proc_distance(DistributionMetricMemo* memo, int64_t dp_id_1, int64_t dp_id_2);

DistributionMetricMemo* new_kl_divergence_memo(HierarchicalDirichletProcess* hdp);
DistributionMetricMemo* new_hellinger_distance_memo(HierarchicalDirichletProcess* hdp);
DistributionMetricMemo* new_l2_distance_memo(HierarchicalDirichletProcess* hdp);
DistributionMetricMemo* new_shannon_jensen_distance_memo(HierarchicalDirichletProcess* hdp);
// note: the lifetime of a DistributionMetricMemo is tied to the lifetime of the
// HierarchicalDirichletProcess that generated it

// computing distances between HDPs

double compare_hdp_distrs_kl_divergence(HierarchicalDirichletProcess* hdp_1, int64_t dp_id_1,
                                        HierarchicalDirichletProcess* hdp_2, int64_t dp_id_2);

double compare_hdp_distrs_l2_distance(HierarchicalDirichletProcess* hdp_1, int64_t dp_id_1,
                                      HierarchicalDirichletProcess* hdp_2, int64_t dp_id_2);

double compare_hdp_distrs_shannon_jensen_distance(HierarchicalDirichletProcess* hdp_1, int64_t dp_id_1,
                                                  HierarchicalDirichletProcess* hdp_2, int64_t dp_id_2);

double compare_hdp_distrs_hellinger_distance(HierarchicalDirichletProcess* hdp_1, int64_t dp_id_1,
                                             HierarchicalDirichletProcess* hdp_2, int64_t dp_id_2);

// serialization

// note: only allowed for HDPs with finalized structure
void serialize_hdp(HierarchicalDirichletProcess* hdp, FILE* out);
HierarchicalDirichletProcess* deserialize_hdp(FILE* in);

#endif // HDP_H_INCLUDED

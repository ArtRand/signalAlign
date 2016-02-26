#ifndef HDP_MATH_UTILS_H_INCLUDED
#define HDP_MATH_UTILS_H_INCLUDED

#include <stdbool.h>
#include <inttypes.h>
#include "sonLib.h"

void parallel_cdf(double* cdf, double* probs, int64_t length, int64_t chunk_size);
double parallel_max(double* x, int64_t length);
void parallel_add(double add_val, double* x, int64_t length);
void parallel_exp(double* x, int64_t length);

typedef struct SumOfLogsMemo SumOfLogsMemo;

SumOfLogsMemo* new_log_sum_memo();
void destroy_log_sum_memo(SumOfLogsMemo* memo);

// returns log(Gamma(n / 2)) in amortized constant time with low risk of overflow
double log_gamma_half(int64_t n, SumOfLogsMemo* sum_of_logs_memo);
double sum_of_logs(SumOfLogsMemo* memo, int64_t n);

// returns log(x + y) without leaving log transformed space
double add_logs(double log_x, double log_y);
// quick-select on array copy (does not alter original array)
double median(double* arr, int64_t length);
double max(double* arr, int64_t length);

// returns the index of the first element of arr greater or equal to x, assuming arr is sorted
// returns final index if x is greater than all elements of arr
int64_t bisect_left(double x, double* arr, int64_t length);

double* spline_knot_slopes(double* x, double* y, int64_t length);
double spline_interp(double query_x, double* x, double* y, double* slope, int64_t length);
double grid_spline_interp(double query_x, double* x, double* y, double* slope, int64_t length);

double* linspace(double start, double stop, int64_t length);

double rand_uniform(double a);
double rand_beta(double a, double b);
bool rand_bernoulli(double p);

// explained further in Jordan's math notebook section "Cached variables for improved performance"
double log_posterior_conditional_term(double nu_post, double two_alpha_post, double beta_post);//,
//SumOfLogsMemo* memo);


void normal_inverse_gamma_params(double* x, int64_t length, double* mu_out, double* nu_out,
                                 double* alpha_out, double* beta_out);
void mle_normal_inverse_gamma_params(double* mus, double* taus, int64_t length, double* mu_0_out,
                                     double* nu_out, double* alpha_out, double* beta_out);

int64_t* stList_toIntPtr(stList* list, int64_t* length_out);
double* stList_toDoublePtr(stList* list, int64_t* length_out);

#endif // HDP_MATH_UTILS_H_INCLUDED

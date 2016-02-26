#include <math.h>
#include <tgmath.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <stdbool.h>
#include <inttypes.h>
#include "hdp_math_utils.h"
#include "sonLib.h"

#define LOG_ROOT_PI 0.572364942924700087071713
#define LOG_4 1.386294361119890618834464

#ifndef M_PI
#define M_PI 3.14159265358979323846264338
#endif

#ifndef EULER_MASCHERONI
#define EULER_MASCHERONI 0.57721566490153286060651209008240243
#endif

#ifndef MACHEP
#define MACHEP 1.11022302462515654042E-16
#endif

#ifndef MINUS_INF
#define MINUS_INF -0.5 * DBL_MAX
#endif

void parallel_cdf(double* cdf, double* probs, int64_t length, int64_t chunk_size) {
    
    if (2 * chunk_size >= length) {
        double cumul = 0.0;
        for (int64_t i = 0; i < length; i++) {
            cumul += probs[i];
            cdf[i] = cumul;
        }
        return;
    }
    
    int64_t num_chunks = (length - 1) / chunk_size + 1;
    
#pragma omp parallel for shared(cdf,probs)
    for (int64_t i = 0; i < num_chunks; i++) {
        int64_t start = i * chunk_size;
        int64_t stop = start + chunk_size;
        if (stop > length) {
            stop = length;
        }
        
        double partial_cumul = 0.0;
        for (int64_t j = start; j < stop; j++) {
            partial_cumul += probs[j];
            cdf[j] = partial_cumul;
        }
    }
    
    double* partial_sums = (double*) malloc(sizeof(double) * num_chunks);
    double partial_sums_cumul = 0.0;
    for (int64_t i = chunk_size - 1; i < length; i += chunk_size) {
        partial_sums_cumul += cdf[i];
        partial_sums[i / chunk_size] = partial_sums_cumul;
    }
    
#pragma omp parallel for shared(cdf,partial_sums)
    for (int64_t i = chunk_size; i < length; i++) {
        cdf[i] += partial_sums[i / chunk_size - 1];
    }
    
    free(partial_sums);
}

double parallel_max(double* x, int64_t length) {
    double max_val = MINUS_INF;
#pragma omp parallel shared(max_val)
    {
        double local_max = MINUS_INF;
#pragma omp for nowait
        for (int64_t i = 0; i < length; i++) {
            if (x[i] > local_max) {
                local_max = x[i];
            }
        }
#pragma omp critical
        {
            if (local_max > max_val) {
                max_val = local_max;
            }
        }
    }
    return max_val;
}

void parallel_add(double add_val, double* x, int64_t length) {
#pragma omp parallel for
    for (int64_t i = 0; i < length; i++) {
        x[i] += add_val;
    }
}

void parallel_exp(double* x, int64_t length) {
#pragma omp parallel for
    for (int64_t i = 0; i < length; i++) {
        x[i] = exp(x[i]);
    }
}

typedef struct LogGammaHalfMemo LogGammaHalfMemo;

struct LogGammaHalfMemo {
    double alpha;
    double* zero_offset_memo;
    int64_t zero_offset_final_entry;
    int64_t zero_offset_length;
    double* half_offset_memo;
    int64_t half_offset_final_entry;
    int64_t half_offset_length;
};

LogGammaHalfMemo* new_log_gamma_memo(double alpha) {
    LogGammaHalfMemo* memo = (LogGammaHalfMemo*) malloc(sizeof(LogGammaHalfMemo));
    memo->alpha = alpha;
    double* zero_base_case = (double*) malloc(sizeof(double));
    zero_base_case[0] = lgamma(alpha);
    memo->zero_offset_final_entry = 0;
    memo->zero_offset_memo = zero_base_case;
    memo->zero_offset_length = 1;
    
    double* half_base_case = (double*) malloc(sizeof(double));
    half_base_case[0] = lgamma(alpha + .5);
    memo->half_offset_final_entry = 0;
    memo->half_offset_memo = half_base_case;
    memo->half_offset_length = 1;
    
    return memo;
}

void destroy_log_gamma_memo(LogGammaHalfMemo* memo) {
    free(memo->half_offset_memo);
    free(memo->zero_offset_memo);
    free(memo);
}

void extend_gamma_zero_offset_memo(LogGammaHalfMemo* memo) {
    int64_t final_entry = memo->half_offset_final_entry + 1;
    memo->zero_offset_final_entry = final_entry;
    double* current_array = memo->zero_offset_memo;
    
    int64_t current_length = memo->zero_offset_length;
    if (current_length == final_entry) {
        
        int64_t new_array_length = current_length * 2;
        double* new_array = (double*) malloc(sizeof(double) * new_array_length);
        
        for (int64_t i = 0; i < current_length; i++) {
            new_array[i] = current_array[i];
        }
        
        memo->zero_offset_length = new_array_length;
        memo->zero_offset_memo = new_array;
        free(current_array);
        current_array = new_array;
    }
    double log_term = log(memo->alpha - 1.0 + (double) final_entry);
    current_array[final_entry] = current_array[final_entry - 1] + log_term;
}

void extend_gamma_half_offset_memo(LogGammaHalfMemo* memo) {
    int64_t final_entry = memo->half_offset_final_entry + 1;
    memo->half_offset_final_entry = final_entry;
    double* current_array = memo->half_offset_memo;
    
    int64_t current_length = memo->half_offset_length;
    if (current_length == final_entry) {
        
        int64_t new_array_length = current_length * 2;
        double* new_array = (double*) malloc(sizeof(double) * new_array_length);
        
        for (int64_t i = 0; i < current_length; i++) {
            new_array[i] = current_array[i];
        }
        
        memo->half_offset_length = new_array_length;
        memo->half_offset_memo = new_array;
        free(current_array);
        current_array = new_array;
    }
    double log_term = log(memo->alpha -.5 + (double) final_entry);
    current_array[final_entry] = current_array[final_entry - 1] + log_term;
}

// returns log(Gamma(memo->alpha + n / 2))
double offset_log_gamma_half(int64_t n, LogGammaHalfMemo* memo) {
    int64_t idx = n / 2;
    if (n % 2 == 0) {
        while (memo->zero_offset_final_entry < idx) {
            extend_gamma_zero_offset_memo(memo);
        }
        return memo->zero_offset_memo[idx];
    }
    else {
        while (memo->half_offset_final_entry < idx) {
            extend_gamma_half_offset_memo(memo);
        }
        return memo->half_offset_memo[idx];
    }
}

struct SumOfLogsMemo {
    double* memo_array;
    int64_t final_entry;
    int64_t array_length;
};

SumOfLogsMemo* new_log_sum_memo() {
    SumOfLogsMemo* memo = (SumOfLogsMemo*) malloc(sizeof(SumOfLogsMemo));
    double* base_case = (double*) malloc(sizeof(double));
    base_case[0] = 0.0;
    memo->memo_array = base_case;
    memo->final_entry = 1;
    memo->array_length = 1;
    return memo;
}

void destroy_log_sum_memo(SumOfLogsMemo* memo) {
    free(memo->memo_array);
    free(memo);
}

void extend_log_sum_memo(SumOfLogsMemo* memo) {
    int64_t final_entry = memo->final_entry;
    int64_t current_length = memo->array_length;
    if (current_length == final_entry) {
        double* current_array = memo->memo_array;
        
        int64_t new_array_length = current_length * 2;
        double* new_array = (double*) malloc(sizeof(double) * new_array_length);
        
        for (int64_t i = 0; i < current_length; i++) {
            new_array[i] = current_array[i];
        }
        
        memo->array_length = new_array_length;
        memo->memo_array = new_array;
        free(current_array);
    }
    double log_term = log((double) final_entry + 1);
    memo->memo_array[final_entry] = memo->memo_array[final_entry - 1] + log_term;
    (memo->final_entry)++;
}

double sum_of_logs(SumOfLogsMemo* memo, int64_t n) {
    while (n > memo->final_entry) {
        extend_log_sum_memo(memo);
    }
    return memo->memo_array[n - 1];
}

// returns log(Gamma(n / 2)) in amortized constant time with low risk of overflow
double log_gamma_half(int64_t n, SumOfLogsMemo* sum_of_logs_memo) {
    if (n <= 2) {
        fprintf(stderr, "log_gamma_half only supports n > 2\n");
        exit(EXIT_FAILURE);
    }
    if (n % 2 == 0) {
        return sum_of_logs(sum_of_logs_memo, n / 2 - 1);
    }
    else {
        return LOG_ROOT_PI - (n / 2) * LOG_4 + sum_of_logs(sum_of_logs_memo, n - 1)
        - sum_of_logs(sum_of_logs_memo, n / 2);
    }
}

// returns log(x + y) without leaving log transformed space
double add_logs(double log_x, double log_y) {
    if (log_x > log_y) {
        return log_x + log(1.0 + exp(log_y - log_x));
    }
    else {
        return log_y + log(1.0 + exp(log_x - log_y));
    }
}

// quick-select algorithm on array copy (does not alter original array)
double quickselect(double* arr, int64_t length, int64_t target_idx) {
    if (target_idx < 0 || target_idx >= length) {
        fprintf(stderr, "Order statistic outside of array bounds\n");
        exit(EXIT_FAILURE);
    }
    
    double* arr_copy = (double*) malloc(sizeof(double) * length);
    for (int64_t i = 0; i < length; i++ ) {
        arr_copy[i] = arr[i];
    }
    
    int64_t low = 0;
    int64_t hi = length - 1;
    int64_t mid;
    int64_t median;
    double temp;
    
    while (true) {
        // median of three technique
        mid = (hi + low) / 2;
        if (arr_copy[hi] > arr_copy[mid]) {
            if (arr_copy[hi] > arr_copy[low]) {
                if (arr_copy[mid] > arr_copy[low]) {
                    median = mid;
                }
                else {
                    median = low;
                }
            }
            else {
                median = hi;
            }
        }
        else {
            if (arr_copy[hi] > arr_copy[low]) {
                median = hi;
            }
            else {
                if (arr_copy[mid] > arr_copy[low]) {
                    median = low;
                }
                else {
                    median = mid;
                }
            }
        }
        
        // remove pivot
        temp = arr_copy[median];
        arr_copy[median] = arr_copy[hi];
        arr_copy[hi] = temp;
        
        // partition array
        int64_t pivot = low;
        for (int64_t i = low; i < hi; i++) {
            if (arr_copy[i] < arr_copy[hi]) {
                temp = arr_copy[i];
                arr_copy[i] = arr_copy[pivot];
                arr_copy[pivot] = temp;
                pivot++;
            }
        }
        
        temp = arr_copy[pivot];
        arr_copy[pivot] = arr_copy[hi];
        arr_copy[hi] = temp;
        
        if (pivot == target_idx) {
            return arr_copy[pivot];
        }
        else if (pivot < target_idx) {
            low = pivot + 1;
        }
        else {
            hi = pivot - 1;
        }
    }
}

double median(double* arr, int64_t length) {
    return quickselect(arr, length, length / 2);
}

double max(double* arr, int64_t length) {
    double curr_max = arr[0];
    for (int64_t i = 1; i < length; i++) {
        if (arr[i] > curr_max) {
            curr_max = arr[i];
        }
    }
    return curr_max;
}

// returns the index of the first element of arr greater or equal to x, assuming arr is sorted
// returns final index if x is greater than all elements of arr
int64_t bisect_left(double x, double* arr, int64_t length) {
    if (x <= arr[0]) {
        return 0;
    }
    int64_t low = 0;
    int64_t hi = length - 1;
    int64_t mid;
    double arr_mid;
    while (hi > low + 1) {
        mid = (hi + low) / 2;
        
        arr_mid = arr[mid];
        if (x <= arr_mid) {
            hi = mid;
        }
        else {
            low = mid;
        }
    }
    return hi;
}

void spline_knot_slopes_internal(double* x, double* y, double* k, int64_t idx, double center_coef_prev,
                                 double right_coef_prev, double rhs_prev, int64_t final_idx) {
    
    if (idx == final_idx) {
        double left_coef = 1.0 / (x[idx] - x[idx - 1]);
        double center_coef = 2.0 * left_coef;
        double rhs = 3.0 * (y[idx] - y[idx - 1]) * left_coef * left_coef;
        // Cramer's rule
        k[idx] = (rhs * center_coef_prev - rhs_prev * left_coef) /
        (center_coef * center_coef_prev - right_coef_prev * left_coef);
        return;
    }
    
    double left_coef = 1.0 / (x[idx] - x[idx - 1]);
    double right_coef = 1.0 / (x[idx + 1] - x[idx]);
    double center_coef = 2.0 * (left_coef + right_coef);
    
    double rhs = 3.0 * ((y[idx] - y[idx - 1]) * left_coef * left_coef +
                        (y[idx + 1] - y[idx]) * right_coef * right_coef);
    
    center_coef -= left_coef * right_coef_prev / center_coef_prev;
    rhs -= left_coef * rhs_prev / center_coef_prev;
    
    spline_knot_slopes_internal(x, y, k, idx + 1, center_coef, right_coef, rhs, final_idx);
    
    k[idx] = (rhs - right_coef * k[idx + 1]) / center_coef;
}

double* spline_knot_slopes(double* x, double* y, int64_t length) {
    double* k = (double*) malloc(sizeof(double) * length);
    
    double right_coef = 1.0 / (x[1] - x[0]);
    double center_coef = 2.0 * right_coef;
    double rhs = 3.0 * (y[1] - y[0]) * right_coef * right_coef;
    
    spline_knot_slopes_internal(x, y, k, 1, center_coef, right_coef, rhs, length - 1);
    
    k[0] = (rhs - right_coef * k[1]) / center_coef;
    
    return k;
}

double spline_interp(double query_x, double* x, double* y, double* slope, int64_t length) {
    if (query_x <= x[0]) {
        return y[0] - slope[0] * (x[0] - query_x);
    }
    else if (query_x >= x[length - 1]) {
        int64_t n = length - 1;
        return y[n] + slope[n] * (query_x - x[n]);
    }
    else {
        int64_t idx_right = bisect_left(query_x, x, length);
        int64_t idx_left = idx_right - 1;
        
        double dx = x[idx_right] - x[idx_left];
        double dy = y[idx_right] - y[idx_left];
        
        double a = slope[idx_left] * dx - dy;
        double b = dy - slope[idx_right] * dx;
        
        double t_left = (query_x - x[idx_left]) / dx;
        double t_right = 1.0 - t_left;
        
        return t_right * y[idx_left] + t_left * y[idx_right] +
        t_left * t_right * (a * t_right + b * t_left);
    }
}

// assumes even spacing of x points
double grid_spline_interp(double query_x, double* x, double* y, double* slope, int64_t length) {
    if (query_x <= x[0]) {
        return y[0] - slope[0] * (x[0] - query_x);
    }
    else if (query_x >= x[length - 1]) {
        int64_t n = length - 1;
        return y[n] + slope[n] * (query_x - x[n]);
    }
    else {
        double dx = x[1] - x[0];
        int64_t idx_left = (int64_t) ((query_x - x[0]) / dx);
        int64_t idx_right = idx_left + 1;
        
        double dy = y[idx_right] - y[idx_left];
        
        double a = slope[idx_left] * dx - dy;
        double b = dy - slope[idx_right] * dx;
        
        double t_left = (query_x - x[idx_left]) / dx;
        double t_right = 1.0 - t_left;
        
        return t_right * y[idx_left] + t_left * y[idx_right]
        + t_left * t_right * (a * t_right + b * t_left);
    }
}

double* linspace(double start, double stop, int64_t length) {
    if (start >= stop) {
        fprintf(stderr, "linspace requires stop > start\n");
        exit(EXIT_FAILURE);
    }
    double* lin = (double*) malloc(sizeof(double) * length);
    int64_t n = length - 1;
    double dx = (stop - start) / ((double) n);
    for (int64_t i = 0; i < n; i++) {
        lin[i] = start +  i * dx;
    }
    lin[n] = stop;
    return lin;
}

double rand_standard_uniform() {
    return ((double) rand()) / ((double) RAND_MAX);
}

double rand_uniform(double a) {
    return ((double) rand()) / ((double) RAND_MAX / a);
}

bool rand_bernoulli(double p) {
    return (rand_standard_uniform() < p);
}

double rand_exponential(double lambda) {
    double draw;
    do {
        draw = rand_standard_uniform();
    } while (draw == 1.0);
    return -log(1.0 - draw) / lambda;
}

double log_posterior_conditional_term(double nu_post, double two_alpha_post,
                                      double beta_post) {//, SumOfLogsMemo* memo) {
    
    //    return log_gamma_half((int64_t) two_alpha_post, memo)
    //           - .5 * (log(nu_post) + two_alpha_post * log(beta_post));
    return lgamma( 0.5 * two_alpha_post) - .5 * (log(nu_post) + two_alpha_post * log(beta_post));
}

void normal_inverse_gamma_params(double* x, int64_t length, double* mu_out, double* nu_out,
                                 double* alpha_out, double* beta_out) {
    double mean = 0.0;
    for (int64_t i = 0; i < length; i++) {
        mean += x[i];
    }
    mean /= (double) length;
    
    double dev;
    double sum_sq_devs = 0.0;
    for (int64_t i = 0; i < length; i++) {
        dev = x[i] - mean;
        sum_sq_devs += dev * dev;
    }
    
    *mu_out = mean;
    *nu_out = (double) length;
    *alpha_out = ((double) length - 1.0) / 2.0;
    *beta_out = .5 * sum_sq_devs;
}

static double A_digamma[] = {
    8.33333333333333333333E-2,
    -2.10927960927960927961E-2,
    7.57575757575757575758E-3,
    -4.16666666666666666667E-3,
    3.96825396825396825397E-3,
    -8.33333333333333333333E-3,
    8.33333333333333333333E-2
};

// modified from Scipy source: https://github.com/scipy/scipy/blob/master/scipy/special/cephes/psi.c
static double polevl(double x, double coef[], int N)
{
    double ans;
    int i;
    double *p;
    
    p = coef;
    ans = *p++;
    i = N;
    
    do
        ans = ans * x + *p++;
    while (--i);
    
    return (ans);
}

double digamma(double x) {
    double p, q, nz, s, w, y, z;
    int i, n, negative;
    
    negative = 0;
    nz = 0.0;
    
    if (x <= 0.0) {
        negative = 1;
        q = x;
        p = floor(q);
        if (p == q) {
            fprintf(stderr, "Digamma evaluated at singularity.\n");
            exit(EXIT_FAILURE);
        }
        /* Remove the zeros of tan(NPY_PI x)
         * by subtracting the nearest integer from x
         */
        nz = q - p;
        if (nz != 0.5) {
            if (nz > 0.5) {
                p += 1.0;
                nz = q - p;
            }
            nz = M_PI / tan(M_PI * nz);
        }
        else {
            nz = 0.0;
        }
        x = 1.0 - x;
    }
    
    /* check for positive integer up to 10 */
    if ((x <= 10.0) && (x == floor(x))) {
        y = 0.0;
        n = x;
        for (i = 1; i < n; i++) {
            w = i;
            y += 1.0 / w;
        }
        y -= EULER_MASCHERONI;
        goto digamma_done;
    }
    
    s = x;
    w = 0.0;
    while (s < 10.0) {
        w += 1.0 / s;
        s += 1.0;
    }
    
    if (s < 1.0e17) {
        z = 1.0 / (s * s);
        y = z * polevl(z, A_digamma, 6);
    }
    else
        y = 0.0;
    
    y = log(s) - (0.5 / s) - y - w;
    
digamma_done:
    
    if (negative) {
        y -= nz;
    }
    
    return y;
}

// modified from SciPy source: https://github.com/scipy/scipy/blob/master/scipy/special/cephes/zeta.c
static double A_zeta[] = {
    12.0,
    -720.0,
    30240.0,
    -1209600.0,
    47900160.0,
    -1.8924375803183791606e9,	/*1.307674368e12/691 */
    7.47242496e10,
    -2.950130727918164224e12,	/*1.067062284288e16/3617 */
    1.1646782814350067249e14,	/*5.109094217170944e18/43867 */
    -4.5979787224074726105e15,	/*8.028576626982912e20/174611 */
    1.8152105401943546773e17,	/*1.5511210043330985984e23/854513 */
    -7.1661652561756670113e18	/*1.6938241367317436694528e27/236364091 */
};

double hurwitz_zeta(double x, double q)
{
    int i;
    double a, b, k, s, t, w;
    
    if (x == 1.0)
        goto retinf;
    
    if (x < 1.0) {
    domerr:
        fprintf(stderr, "Domain error in zeta function.\n");
        exit(EXIT_FAILURE);
    }
    
    if (q <= 0.0) {
        if (q == floor(q)) {
        retinf:
            fprintf(stderr, "Evaluted zeta function at singularity.\n");
            exit(EXIT_FAILURE);
        }
        if (x != floor(x))
            goto domerr;	/* because q^-x not defined */
    }
    
    /* Asymptotic expansion
     * http://dlmf.nist.gov/25.11#E43
     */
    if (q > 1e8) {
        return (1/(x - 1) + 1/(2*q)) * pow(q, 1 - x);
    }
    
    /* Euler-Maclaurin summation formula */
    
    /* Permit negative q but continue sum until n+q > +9 .
     * This case should be handled by a reflection formula.
     * If q<0 and x is an integer, there is a relation to
     * the polyGamma function.
     */
    s = pow(q, -x);
    a = q;
    i = 0;
    b = 0.0;
    while ((i < 9) || (a <= 9.0)) {
        i += 1;
        a += 1.0;
        b = pow(a, -x);
        s += b;
        if (fabs(b / s) < MACHEP)
            goto zeta_done;
    }
    
    w = a;
    s += b * w / (x - 1.0);
    s -= 0.5 * b;
    a = 1.0;
    k = 0.0;
    for (i = 0; i < 12; i++) {
        a *= x + k;
        b /= w;
        t = a * b / A_zeta[i];
        s = s + t;
        t = fabs(t / s);
        if (t < MACHEP)
            goto zeta_done;
        k += 1.0;
        a *= x + k;
        b /= w;
        k += 1.0;
    }
zeta_done:
    return (s);
}

double trigamma(double x) {
    return hurwitz_zeta(2.0, x);
}

double newton_approx_alpha(int64_t length, double sum_log_tau, double sum_tau) {
    double constant = sum_log_tau / length - log( sum_tau / length);
    
    double alpha = (double) 1.0;
    
    double f_alpha;
    double df_alpha;
    double alpha_prime;
    while (true) {
        f_alpha = log(alpha) - digamma(alpha) + constant;
        df_alpha = 1.0 / alpha - trigamma(alpha);
        
        if (df_alpha == 0.0 || df_alpha != df_alpha) {
            fprintf(stderr, "MLE estimation of alpha numerically unstable at designated starting value.\n");
            exit(EXIT_FAILURE);
        }
        
        alpha_prime = alpha - f_alpha / df_alpha;
        if (fabs(alpha - alpha_prime) < MACHEP) {
            return alpha_prime;
        }
        alpha = alpha_prime;
    }
}


void mle_normal_inverse_gamma_params(double* mus, double* taus, int64_t length, double* mu_0_out,
                                     double* nu_out, double* alpha_out, double* beta_out) {
    double sum_tau = 0.0;
    double sum_log_tau = 0.0;
    for (int64_t i = 0; i < length; i++) {
        sum_tau += taus[i];
        sum_log_tau += log(taus[i]);
    }
    
    double mu_0 = 0.0;
    for (int64_t i = 0; i < length; i++) {
        mu_0 += mus[i] * taus[i];
    }
    
    mu_0 /= sum_tau;
    
    double sum_weighted_sq_devs = 0.0;
    double dev;
    for (int64_t i = 0; i < length; i++) {
        dev = mus[i] - mu_0;
        sum_weighted_sq_devs += taus[i] * dev * dev;
    }
    
    double nu = ((double) length) / sum_weighted_sq_devs;
    
    double alpha = newton_approx_alpha(length, sum_log_tau, sum_tau);
    
    double beta = length * alpha / sum_tau;
    
    *mu_0_out = mu_0;
    *nu_out = nu;
    *alpha_out = alpha;
    *beta_out = beta;
}

int64_t* stList_toIntPtr(stList* list, int64_t* length_out) {
    int64_t length = (int64_t) stList_length(list);
    int64_t* int_arr = (int64_t*) malloc(sizeof(int64_t) * length);
    int64_t* entry;
    for (int64_t i = 0; i < length; i++) {
        entry = (int64_t*) stList_get(list, i);
        int_arr[i] = *entry;
    }
    *length_out = length;
    return int_arr;
}

double* stList_toDoublePtr(stList* list, int64_t* length_out) {
    int64_t length  = stList_length(list);
    double* double_arr = (double*) malloc(sizeof(double) * length);
    double* entry;
    for (int64_t i = 0; i < length; i++) {
        entry = (double*) stList_get(list, i);
        double_arr[i] = *entry;
    }
    *length_out = length;
    return double_arr;
}

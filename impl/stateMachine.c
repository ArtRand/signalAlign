/*
 * stateMachine.c
 *
 *  Created on: 1 Aug 2014
 *      Author: benedictpaten
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include "stateMachine.h"
#include "bioioC.h"
#include "pairwiseAligner.h"
#include "discreteHmm.h"

//////////////////////////////////////////////////////////////////////////////
// StateMachine Emission functions for discrete alignments (symbols and kmers)
//////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////// STATIC FUNCTIONS ////////////////////////////////////////////////////////
static void state_check(StateMachine *sM, State s) {
    assert(s >= 0 && s < sM->stateNumber);
}


void stateMachine_index_check(int64_t c) {
    assert(c >= 0 && c < NUM_OF_KMERS);
    if ((c < 0) || (c > NUM_OF_KMERS)) {
        st_errAbort("stateMachine_index_check: Got invalid kmerIndex %lld\n", c);
    }
}

static inline void emissions_discrete_initializeEmissionsMatrices(StateMachine *sM) {
    sM->EMISSION_GAP_X_PROBS = st_malloc(sM->parameterSetSize*sizeof(double));
    sM->EMISSION_GAP_Y_MATRIX = st_malloc(sM->parameterSetSize * sizeof(double));
    sM->EMISSION_MATCH_MATRIX = st_malloc(sM->parameterSetSize * sM->parameterSetSize * sizeof(double));
}
/*
void emissions_symbol_setEmissionsToDefaults(StateMachine *sM) {
    // initialize
    emissions_discrete_initializeEmissionsMatrices(sM);

    // Set Match probs to default values
    const double EMISSION_MATCH=-2.1149196655034745; //log(0.12064298095701059);
    const double EMISSION_TRANSVERSION=-4.5691014376830479; //log(0.010367271172731285);
    const double EMISSION_TRANSITION=-3.9833860032220842; //log(0.01862247669752685);

    const double M[SYMBOL_NUMBER_NO_N*SYMBOL_NUMBER_NO_N] = {
            EMISSION_MATCH, EMISSION_TRANSVERSION, EMISSION_TRANSITION, EMISSION_TRANSVERSION,
            EMISSION_TRANSVERSION, EMISSION_MATCH, EMISSION_TRANSVERSION, EMISSION_TRANSITION,
            EMISSION_TRANSITION, EMISSION_TRANSVERSION, EMISSION_MATCH, EMISSION_TRANSVERSION,
            EMISSION_TRANSVERSION, EMISSION_TRANSITION, EMISSION_TRANSVERSION, EMISSION_MATCH };

    memcpy(sM->EMISSION_MATCH_MATRIX, M, sizeof(double)*SYMBOL_NUMBER_NO_N*SYMBOL_NUMBER_NO_N);

    // Set Gap probs to default values
    const double EMISSION_GAP = -1.6094379124341003; //log(0.2)
    const double G[4] = { EMISSION_GAP, EMISSION_GAP, EMISSION_GAP, EMISSION_GAP };
    memcpy(sM->EMISSION_GAP_X_PROBS, G, sizeof(double)*SYMBOL_NUMBER_NO_N);
    memcpy(sM->EMISSION_GAP_Y_MATRIX, G, sizeof(double)*SYMBOL_NUMBER_NO_N);
}
*/
static inline void emissions_discrete_initMatchProbsToZero(double *emissionMatchProbs, int64_t symbolSetSize) {
    memset(emissionMatchProbs, 0, symbolSetSize*symbolSetSize*sizeof(double));
}

static inline void emissions_discrete_initGapProbsToZero(double *emissionGapProbs, int64_t symbolSetSize) {
    memset(emissionGapProbs, 0, symbolSetSize*sizeof(double));
}

///////////////////////////////////////////// CORE FUNCTIONS ////////////////////////////////////////////////////////

void emissions_discrete_initEmissionsToZero(StateMachine *sM) {
    // initialize
    emissions_discrete_initializeEmissionsMatrices(sM);
    // set match matrix to zeros
    emissions_discrete_initMatchProbsToZero(sM->EMISSION_MATCH_MATRIX, sM->parameterSetSize);
    // set gap matrix to zeros
    emissions_discrete_initGapProbsToZero(sM->EMISSION_GAP_X_PROBS, sM->parameterSetSize);
    emissions_discrete_initGapProbsToZero(sM->EMISSION_GAP_Y_MATRIX, sM->parameterSetSize);
}

int64_t emissions_discrete_getBaseIndex(void *base) {
    char b = *(char*) base;
    switch (b) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'E':
            return 2;
        case 'G':
            return 3;
        case 'O':
            return 4;
        case 'T':
            return 5;
        default: // N or n. Hack to make sure that we get a zero model mean
            return NUM_OF_KMERS + 1;
    }
}

int64_t emissions_discrete_getKmerIndexFromKmer(void *kmer) {
    int64_t kmerLen = strlen((char*) kmer);
    if (kmerLen == 0) {
        return NUM_OF_KMERS + 1;
    }
    int64_t axisLength = NUM_OF_KMERS;
    int64_t l = axisLength / SYMBOL_NUMBER_NO_N;
    int64_t i = 0;
    int64_t x = 0;
    while(l > 1) {
        x += l * emissions_discrete_getBaseIndex((char *)kmer + i);
        i += 1;
        l = l / SYMBOL_NUMBER_NO_N;
    }
    int64_t last = KMER_LENGTH - 1;
    x += emissions_discrete_getBaseIndex((char *)kmer + last);

    return x;
}

int64_t emissions_discrete_getKmerIndexFromPtr(void *kmer) {
    // make temp kmer meant to work with getKmer
    char *kmer_i = malloc((KMER_LENGTH) * sizeof(char));
    for (int64_t x = 0; x < KMER_LENGTH; x++) {
        kmer_i[x] = *((char *)kmer+x);
    }
    kmer_i[KMER_LENGTH] = '\0';
    int64_t i = emissions_discrete_getKmerIndexFromKmer(kmer_i);
    //index_check(i);
    free(kmer_i);
    return i;
}

double emissions_symbol_getGapProb(const double *emissionGapProbs, void *base) {
    int64_t i = emissions_discrete_getBaseIndex(base);
    if(i == 4) {
        return -1.386294361; //log(0.25)
    }
    return emissionGapProbs[i];
}

double emissions_symbol_getMatchProb(const double *emissionMatchProbs, void *x, void *y) {
    int64_t iX = emissions_discrete_getBaseIndex(x);
    int64_t iY = emissions_discrete_getBaseIndex(y);
    if(iX == 4 || iY == 4) {
        return -2.772588722; //log(0.25**2)
    }
    return emissionMatchProbs[iX * SYMBOL_NUMBER_NO_N + iY];
}

double emissions_kmer_getGapProb(StateMachine *sM, const double *emissionGapProbs, void *x_i) {
    // even though this shouldn't happen, check for null x_i and return logzero
    if (x_i == NULL) {
        return LOG_ZERO;
    }
    // make temp x_i meant to work with getKmer
    char *kmer_i = malloc((sM->kmerLength) * sizeof(char));
    for (int64_t x = 0; x < sM->kmerLength; x++) {
        kmer_i[x] = *((char *) x_i + x);
    }
    kmer_i[sM->kmerLength] = '\0';

    int64_t i = kmer_id(kmer_i, sM->alphabet, sM->alphabetSize, sM->kmerLength);
    //index_check(i);
    free(kmer_i);
    double p = i > NUM_OF_KMERS ? LOG_ZERO : emissionGapProbs[i];
    return p;
}

double emissions_kmer_getMatchProb(const double *emissionMatchProbs, void *x, void *y) {
    int64_t iX = emissions_discrete_getKmerIndexFromKmer(x);
    int64_t iY = emissions_discrete_getKmerIndexFromKmer(y);
    int64_t tableIndex = iX * NUM_OF_KMERS + iY;
    return emissionMatchProbs[tableIndex];
}
/////////////////////////////////////////
// functions for signal/kmer alignment //
/////////////////////////////////////////

/////////////////////////////////////////// STATIC FUNCTIONS ////////////////////////////////////////////////////////

static inline void emissions_signal_initializeEmissionsMatrices(StateMachine *sM, int64_t nbSkipParams) {
    // changed to 30 for skip prob bins
    // the kmer/gap and skip (f(|ui-1 - ui|)) probs have smaller tables, either 30 for the skip or parameterSetSize
    // for the naive kmer skip one
    // see note at emissions_signal_initEmissionsToZero about this change
    sM->EMISSION_GAP_X_PROBS = st_malloc(nbSkipParams * sizeof(double));

    // both the Iy and M - type states use the event/kmer match model so the matrices need to be the same size
    sM->EMISSION_GAP_Y_MATRIX = st_malloc((sM->parameterSetSize * MODEL_PARAMS) * sizeof(double));
    sM->EMISSION_MATCH_MATRIX = st_malloc((sM->parameterSetSize * MODEL_PARAMS) * sizeof(double));
}

static inline void emissions_signal_initMatchMatrixToZero(double *matchModel, int64_t parameterSetSize) {
    memset(matchModel, 0, ((parameterSetSize * MODEL_PARAMS)) * sizeof(double));
}

static inline void emissions_signal_initKmerSkipTableToZero(double *skipModel, int64_t parameterSetSize) {
    memset(skipModel, 0, parameterSetSize * sizeof(double));
}

static inline double emissions_signal_getModelLevelMean(const double *eventModel, int64_t kmerIndex) {
    return kmerIndex > NUM_OF_KMERS ? 0.0 : eventModel[(kmerIndex * MODEL_PARAMS)];
}

static inline double emissions_signal_getModelLevelSd(const double *eventModel, int64_t kmerIndex) {
    return kmerIndex > NUM_OF_KMERS ? 0.0 : eventModel[(kmerIndex * MODEL_PARAMS + 1)];
}

static inline double emissions_signal_getModelFluctuationMean(const double *eventModel, int64_t kmerIndex) {
    return kmerIndex > NUM_OF_KMERS ? 0.0 : eventModel[(kmerIndex * MODEL_PARAMS + 2)];
}

static inline double emissions_signal_getModelFluctuationSd(const double *eventModel, int64_t kmerIndex) {
    return kmerIndex > NUM_OF_KMERS ? 0.0 : eventModel[(kmerIndex * MODEL_PARAMS + 3)];
}

static inline double emissions_signal_getModelFluctuationLambda(const double *eventModel, int64_t kmerIndex) {
    return kmerIndex > NUM_OF_KMERS ? 0.0 : eventModel[(kmerIndex * MODEL_PARAMS + 4)];
}

static inline double emissions_signal_logInvGaussPdf(double eventNoise, double modelNoiseMean,
                                                     double modelNoiseLambda) {
    double l_twoPi = 1.8378770664093453;  // log(2*pi)
    double l_eventNoise = log(eventNoise);
    double a = (eventNoise - modelNoiseMean) / modelNoiseMean;
    double l_modelNoseLambda = log(modelNoiseLambda);

    // returns Log-space
    return (l_modelNoseLambda - l_twoPi - 3 * l_eventNoise - modelNoiseLambda * a * a / eventNoise) / 2;
}

static inline double emissions_signal_logGaussPdf(double x, double mu, double sigma) {
    if (sigma == 0.0) {
        return LOG_ZERO;
    }
    double log_inv_sqrt_2pi = -0.91893853320467267;
    double l_sigma = log(sigma);
    double a = (x - mu) / sigma;

    // returns Log-space
    return log_inv_sqrt_2pi - l_sigma + (-0.5 * a * a);
}

static double emissions_signal_poissonPosteriorProb(int64_t n, double duration) {
    assert(n <= 5);

    // Experimented with different values of c,
    //double c = 0.00570570570571; // max of PDF
    double c = 0.00332005312085; // mode of all test durations
    //double c = 0.0045; // guess of somewhere in between

    // Experimenting with changing the rate parameter, started with 2, but then the p(c|N=2) == p(c|N=1) which
    // doesn't make sense. When 0 < beta < 1 then the p(c|N=0) > p(c|N=1). At 1.25, it seems to have the correct
    // curve.
    //double l_beta = 0.22314355131420976; // log(1.25)
    double l_beta = 0.1397619423751586; // log(1.15)
    double lambda = duration / c;
    double l_factorials[6] = {0.0, 0.0, 0.69314718056, 1.79175946923, 3.17805383035, 4.78749174278};

    //result = ((n+1)*np.log(2)) + (n*np.log(lam)) - np.log(factorial(n)) - (2*lam)
    double a = (n+1) * l_beta;
    double b = n * log(lambda);
    double d = 2 * lambda;

    // returns log-space
    //double prob = a + b - l_factorials[n] - d;
    //return prob;
    return a + b - l_factorials[n] - d;
}

///////////////////////////////////////////// CORE FUNCTIONS ////////////////////////////////////////////////////////
double emissions_signal_descaleEventMean_JordanStyle(double scaledEvent, double levelMean,
                                                     double scale, double shift, double var) {
    // (x  + var * level_mean - scale * level_mean - shift) / var
    return (scaledEvent + var * levelMean - scale * levelMean - shift) / var;
}

void emissions_signal_initEmissionsToZero(StateMachine *sM, int64_t nbSkipParams) {
    // initialize
    emissions_signal_initializeEmissionsMatrices(sM, nbSkipParams);

    // set kmer skip matrix to zeros
    emissions_signal_initKmerSkipTableToZero(sM->EMISSION_GAP_X_PROBS, nbSkipParams);

    // set extra event matrix to zeros
    emissions_signal_initMatchMatrixToZero(sM->EMISSION_GAP_Y_MATRIX, sM->parameterSetSize);

    // set match matrix to zeros
    emissions_signal_initMatchMatrixToZero(sM->EMISSION_MATCH_MATRIX, sM->parameterSetSize);
}

int64_t emissions_signal_getKmerSkipBin(double *matchModel, void *kmers) {
    char *kmer_im1 = malloc((KMER_LENGTH) * sizeof(char));
    for (int64_t x = 0; x < KMER_LENGTH; x++) {
        kmer_im1[x] = *((char *)kmers+x);
    }
    kmer_im1[KMER_LENGTH] = '\0';

    // make kmer_i
    char *kmer_i = malloc((KMER_LENGTH) * sizeof(char));
    for (int64_t x = 0; x < KMER_LENGTH; x++) {
        kmer_i[x] = *((char *)kmers+(x+1));
    }
    kmer_i[KMER_LENGTH] = '\0';

    // get indices
    int64_t k_i= emissions_discrete_getKmerIndexFromKmer(kmer_i);
    int64_t k_im1 = emissions_discrete_getKmerIndexFromKmer(kmer_im1);

    // get the expected mean current for each one
    double u_ki = emissions_signal_getModelLevelMean(matchModel, k_i);
    double u_kim1 = emissions_signal_getModelLevelMean(matchModel, k_im1);

    // find the difference
    double d = fabs(u_ki - u_kim1);

    // get the 'bin' for skip prob, clamp to the last bin
    int64_t bin = (int64_t)(d / 0.5); // 0.5 pA bins right now
    bin = bin >= 30 ? 29 : bin;
    free(kmer_im1);
    free(kmer_i);
    return bin;
}

double emissions_signal_getBetaOrAlphaSkipProb(StateMachine *sM, void *kmers, bool getAlpha) {
    // downcast
    //StateMachine3Vanilla *sM3v = (StateMachine3Vanilla *) sM;
    // get the skip bin
    int64_t bin = emissions_signal_getKmerSkipBin(sM->EMISSION_MATCH_MATRIX, kmers);
    return getAlpha ? sM->EMISSION_GAP_X_PROBS[bin+30] : sM->EMISSION_GAP_X_PROBS[bin];
}

double emissions_signal_logGaussMatchProb(const double *eventModel, void *kmer, void *event) {
    // make temp kmer meant to work with getKmer2
    char *kmer_i = malloc((KMER_LENGTH) * sizeof(char));
    for (int64_t x = 0; x < KMER_LENGTH; x++) {
        kmer_i[x] = *((char *)kmer+(x+1));
    }
    kmer_i[KMER_LENGTH] = '\0';

    // get event mean, and kmer index
    double eventMean = *(double *) event;
    int64_t kmerIndex = emissions_discrete_getKmerIndexFromKmer(kmer_i);
    double l_inv_sqrt_2pi = log(0.3989422804014327); // constant
    double modelMean = emissions_signal_getModelLevelMean(eventModel, kmerIndex);
    double modelStdDev = emissions_signal_getModelLevelSd(eventModel, kmerIndex);
    double l_modelSD = log(modelStdDev);
    double a = (eventMean - modelMean) / modelStdDev;

    // clean up
    free(kmer_i);
    /// / debugging
    //double prob = l_inv_sqrt_2pi - l_modelSD + (-0.5f * a * a);
    //st_uglyf("MATCHING--kmer:%s (index: %lld), event mean: %f, levelMean: %f, prob: %f\n", kmer_i, kmerIndex, eventMean, modelMean, prob);
    // returns log space
    return l_inv_sqrt_2pi - l_modelSD + (-0.5f * a * a);
}

double emissions_signal_getEventMatchProbWithTwoDists(const double *eventModel, void *kmer, void *event) {
    // make temp kmer meant to work with getKmer2
    char *kmer_i = malloc((KMER_LENGTH) * sizeof(char));
    for (int64_t x = 0; x < KMER_LENGTH; x++) {
        kmer_i[x] = *((char *)kmer+(x+1));
    }
    kmer_i[KMER_LENGTH] = '\0';

    // get event mean, and noise
    double eventMean = *(double *) event;
    double eventNoise = *(double *) ((char *)event + sizeof(double));

    // get the kmer index
    int64_t kmerIndex = emissions_discrete_getKmerIndexFromKmer(kmer_i);

    // first calculate the prob of the level mean
    double expectedLevelMean = emissions_signal_getModelLevelMean(eventModel, kmerIndex);
    double expectedLevelSd = emissions_signal_getModelLevelSd(eventModel, kmerIndex);
    double levelProb = emissions_signal_logGaussPdf(eventMean, expectedLevelMean, expectedLevelSd);

    // now calculate the prob of the noise mean
    double expectedNoiseMean = emissions_signal_getModelFluctuationMean(eventModel, kmerIndex);
    double modelNoiseLambda = emissions_signal_getModelFluctuationLambda(eventModel, kmerIndex);
    double noiseProb = emissions_signal_logInvGaussPdf(eventNoise, expectedNoiseMean, modelNoiseLambda);

    // clean up
    free(kmer_i);

    return levelProb + noiseProb;
}

double emissions_signal_multipleKmerMatchProb(const double *eventModel, void *kmers, void *event, int64_t n) {
    // this is meant to work with getKmer2
    double p = 0.0;
    for (int64_t i = 0; i < n; i++) {
        // check if we're going to run off the end of the sequence
        int64_t l = (KMER_LENGTH * n);
        //char lastBase = *(char *) (kmers + l);
        char lastBase = *((char *)kmers + l); // new way
        // if we're not, logAdd up all the probs for the next n kmers
        if (isupper(lastBase)) {
            //char *x_i = &kmers[i];
            char *x_i = ((char *)kmers + i); // new way
            p = logAdd(p, emissions_signal_getEventMatchProbWithTwoDists(eventModel, x_i, event));
        } else { // otherwise return zero prob of emitting this many kmers from this event
            return LOG_ZERO;
        }
    }

    return p - log(n);
}

double emissions_signal_getDurationProb(void *event, int64_t n) {
    double duration = *(double *) ((char *)event + (2 * sizeof(double)));
    return emissions_signal_poissonPosteriorProb(n, duration);
}

double emissions_signal_getBivariateGaussPdfMatchProb(const double *eventModel, void *kmer, void *event) {
    // this is meant to work with getKmer2
    // wrangle event data
    double eventMean = *(double *) event;
    double eventNoise = *(double *) ((char*)event + sizeof(double)); // aaah pointers

    // correlation coefficient is the 0th member of the event model
    double p = eventModel[0];
    double pSq = p * p;

    // make temp kmer
    char *kmer_i = malloc((KMER_LENGTH) * sizeof(char));
    for (int64_t x = 0; x < KMER_LENGTH; x++) {
        kmer_i[x] = *((char *)kmer+(x+1));
    }
    kmer_i[KMER_LENGTH] = '\0';

    int64_t kmerIndex = emissions_discrete_getKmerIndexFromKmer(kmer_i);

    // get the µ and σ for the level and noise for the model
    double levelMean = emissions_signal_getModelLevelMean(eventModel, kmerIndex);
    double levelStdDev = emissions_signal_getModelLevelSd(eventModel, kmerIndex);
    double noiseMean = emissions_signal_getModelFluctuationMean(eventModel, kmerIndex);
    double noiseStdDev = emissions_signal_getModelFluctuationSd(eventModel, kmerIndex);

    // do calculation
    double log_inv_2pi = -1.8378770664093453;
    double expC = -1 / (2 * (1 - pSq));
    double xu = (eventMean - levelMean) / levelStdDev;
    double yu = (eventNoise - noiseMean) / noiseStdDev;
    double a = expC * ((xu * xu) + (yu * yu) - (2 * p * xu * yu));
    double c = log_inv_2pi - log(levelStdDev * noiseStdDev * sqrt(1 - pSq));

    // clean
    free(kmer_i);

    return c + a;
}

double emissions_signal_getHdpKmerDensity(StateMachine *sM, void *x_i, void *e_j, bool ignore) {
    (void) ignore;
    StateMachine3_HDP *self = (StateMachine3_HDP *)sM;
    if (x_i == NULL) {
        return LOG_ZERO;
    }

    // make temp x_i
    char *kmer_i = malloc((self->model.kmerLength) * sizeof(char));
    for (int64_t x = 0; x < self->model.kmerLength; x++) {
        kmer_i[x] = *((char *) x_i + x);
    }
    kmer_i[self->model.kmerLength] = '\0';

    // wrangle e_j data
    double eventMean = *(double *)e_j;
    int64_t kmerIndex = kmer_id(kmer_i, self->model.alphabet, self->model.alphabetSize, self->model.kmerLength);

    double levelMean = emissions_signal_getModelLevelMean(self->model.EMISSION_MATCH_MATRIX, kmerIndex);
    double *normedMean = (double *)st_malloc(sizeof(double));
    *normedMean = emissions_signal_descaleEventMean_JordanStyle(eventMean, levelMean,
                                                                self->model.scale, self->model.shift, self->model.var);

    double density = (1 / self->model.var) * get_nanopore_kmer_density(self->hdpModel, x_i, normedMean);
    free(normedMean);
    free(kmer_i);
    return log(density);
}


double emissions_signal_strawManGetKmerEventMatchProbWithDescaling(StateMachine *sM, void *x_i, void *e_j, bool match) {
    StateMachine3 *self = (StateMachine3 *)sM;  // downcast

    if (x_i == NULL) {
        return LOG_ZERO;
    }
    // this is meant to work with getKmer (NOT getKmer2)
    // wrangle e_j data
    double eventMean = *(double *) e_j;
    double eventNoise = *(double *) ((char *) e_j + sizeof(double)); // aaah pointers

    // make temp x_i
    char *kmer_i = malloc((self->model.kmerLength) * sizeof(char));
    for (int64_t x = 0; x < self->model.kmerLength; x++) {
        kmer_i[x] = *((char *) x_i + x);
    }
    kmer_i[self->model.kmerLength] = '\0';

    // get index
    int64_t kmerIndex = kmer_id(kmer_i, self->model.alphabet, self->model.alphabetSize, self->model.kmerLength);

    double *eventModel = match ? self->model.EMISSION_MATCH_MATRIX : self->model.EMISSION_GAP_Y_MATRIX;
    // get the µ and σ for the level and noise for the model
    double levelMean = emissions_signal_getModelLevelMean(eventModel, kmerIndex);
    double levelStdDev = emissions_signal_getModelLevelSd(eventModel, kmerIndex);
    eventMean = emissions_signal_descaleEventMean_JordanStyle(eventMean, levelMean, self->model.scale,
                                                              self->model.shift, self->model.var);

    double noiseMean = emissions_signal_getModelFluctuationMean(eventModel, kmerIndex);
    //double noiseStdDev = emissions_signal_getModelFluctuationSd(eventModel, kmerIndex);

    double modelNoiseLambda = emissions_signal_getModelFluctuationLambda(eventModel, kmerIndex);

    double l_probEventMean = emissions_signal_logGaussPdf(eventMean, levelMean, levelStdDev);
    // I tried using the inverse gaussian here
    //double l_probEventNoise = emissions_signal_logGaussPdf(eventNoise, noiseMean, noiseStdDev);
    double l_probEventNoise = emissions_signal_logInvGaussPdf(eventNoise, noiseMean, modelNoiseLambda);

    // clean
    free(kmer_i);

    // debugging
    //double prob = l_probEventMean + l_probEventNoise;
    //st_uglyf("MATCHING--x_i:%s (index: %lld), e_j mean: %f, \n modelMean: %f, modelLsd: %f probEvent: %f probNoise: %f, combined: %f\n",
    //         kmer_i, kmerIndex, eventMean, levelMean, levelStdDev, l_probEventMean, l_probEventNoise, prob);

    return l_probEventMean + l_probEventNoise;
}

double emissions_signal_strawManGetKmerEventMatchProb(StateMachine *sM, void *x_i, void *e_j, bool match) {
    StateMachine3 *self = (StateMachine3 *)sM;
    if (x_i == NULL) {
        return LOG_ZERO;
    }
    // this is meant to work with getKmer (NOT getKmer2)
    // wrangle e_j data
    double eventMean = *(double *) e_j;
    double eventNoise = *(double *) ((char *) e_j + sizeof(double)); // aaah pointers

    // make temp x_i
    char *kmer_i = malloc((self->model.kmerLength) * sizeof(char));
    for (int64_t x = 0; x < self->model.kmerLength; x++) {
        kmer_i[x] = *((char *) x_i + x);
    }
    kmer_i[self->model.kmerLength] = '\0';

    // get index
    int64_t kmerIndex = kmer_id(kmer_i, self->model.alphabet, self->model.alphabetSize, self->model.kmerLength);
    double *eventModel = match ? self->model.EMISSION_MATCH_MATRIX : self->model.EMISSION_GAP_Y_MATRIX;

    // get the µ and σ for the level and noise for the model
    double levelMean = emissions_signal_getModelLevelMean(eventModel, kmerIndex);
    double levelStdDev = emissions_signal_getModelLevelSd(eventModel, kmerIndex);
    double noiseMean = emissions_signal_getModelFluctuationMean(eventModel, kmerIndex);
    //double noiseStdDev = emissions_signal_getModelFluctuationSd(eventModel, kmerIndex);

    double modelNoiseLambda = emissions_signal_getModelFluctuationLambda(eventModel, kmerIndex);

    double l_probEventMean = emissions_signal_logGaussPdf(eventMean, levelMean, levelStdDev);
    //double l_probEventNoise = emissions_signal_logGaussPdf(eventNoise, noiseMean, noiseStdDev);
    double l_probEventNoise = emissions_signal_logInvGaussPdf(eventNoise, noiseMean, modelNoiseLambda);

    // clean
    free(kmer_i);

    // debugging
    //double prob = l_probEventMean + l_probEventNoise;
    //st_uglyf("MATCHING--x_i:%s (index: %lld), e_j mean: %f, \n modelMean: %f, modelLsd: %f probEvent: %f probNoise: %f, combined: %f\n",
    //         kmer_i, kmerIndex, eventMean, levelMean, levelStdDev, l_probEventMean, l_probEventNoise, prob);

    return l_probEventMean + l_probEventNoise;
}

void emissions_signal_scaleEmissions(StateMachine *sM, double scale, double shift, double var) {
    // model is arranged: level_mean, level_stdev, sd_mean, sd_stdev, sd_lambda per kmer
    // already been adjusted for correlation coeff.
    for (int64_t i = 0; i < (sM->parameterSetSize * MODEL_PARAMS); i += MODEL_PARAMS) {
        // Level adjustments
        // level_mean = mean * scale + shift
        sM->EMISSION_MATCH_MATRIX[i] = sM->EMISSION_MATCH_MATRIX[i] * scale + shift;
        // level_stdev = stdev * var
        sM->EMISSION_MATCH_MATRIX[i + 1] = sM->EMISSION_MATCH_MATRIX[i + 1] * var;

        // adjusting extra event also
        // level_mean = mean * scale + shift
        sM->EMISSION_GAP_Y_MATRIX[i] = sM->EMISSION_GAP_Y_MATRIX[i] * scale + shift;
        // level_stdev = stdev * var
        sM->EMISSION_GAP_Y_MATRIX[i + 1] = sM->EMISSION_GAP_Y_MATRIX[i + 1] * var;
    }
}

void emissions_signal_scaleNoise(StateMachine *sM, NanoporeReadAdjustmentParameters npp) {
    for (int64_t i = 0; i < (sM->parameterSetSize * MODEL_PARAMS); i += MODEL_PARAMS) {
        // Fluctuation (noise) adjustments
        // noise_mean *= scale_sd
        sM->EMISSION_MATCH_MATRIX[i + 2] = sM->EMISSION_MATCH_MATRIX[i + 2] * npp.scale_sd;
        // noise_lambda *= var_sd
        sM->EMISSION_MATCH_MATRIX[i + 4] = sM->EMISSION_MATCH_MATRIX[i + 4] * npp.var_sd;
        // noise_sd = sqrt(adjusted_noise_mean**3 / adjusted_noise_lambda);
        sM->EMISSION_MATCH_MATRIX[i + 3] = sqrt(pow(sM->EMISSION_MATCH_MATRIX[i + 2], 3.0)
                                                / sM->EMISSION_MATCH_MATRIX[i + 4]);

        // Fluctuation (noise) adjustments
        // noise_mean *= scale_sd
        sM->EMISSION_GAP_Y_MATRIX[i + 2] = sM->EMISSION_GAP_Y_MATRIX[i + 2] * npp.scale_sd;
        // noise_lambda *= var_sd
        sM->EMISSION_GAP_Y_MATRIX[i + 4] = sM->EMISSION_GAP_Y_MATRIX[i + 4] * npp.var_sd;
        // noise_sd = sqrt(adjusted_noise_mean**3 / adjusted_noise_lambda);
        sM->EMISSION_GAP_Y_MATRIX[i + 3] = sqrt(pow(sM->EMISSION_GAP_Y_MATRIX[i + 2], 3.0)
                                                / sM->EMISSION_MATCH_MATRIX[i + 4]);
    }
}

void emissions_signal_scaleModel(StateMachine *sM,
                                 double scale, double shift, double var,
                                 double scale_sd, double var_sd) {
    // model is arranged: level_mean, level_stdev, sd_mean, sd_stdev, sd_lambda per kmer
    // already been adjusted for correlation coeff.
    for (int64_t i = 0; i < (sM->parameterSetSize * MODEL_PARAMS); i += MODEL_PARAMS) {
        // Level adjustments
        // level_mean = mean * scale + shift
        sM->EMISSION_MATCH_MATRIX[i] = sM->EMISSION_MATCH_MATRIX[i] * scale + shift;
        // level_stdev = stdev * var
        sM->EMISSION_MATCH_MATRIX[i + 1] = sM->EMISSION_MATCH_MATRIX[i + 1] * var;

        // adjusting extra event also
        // level_mean = mean * scale + shift
        sM->EMISSION_GAP_Y_MATRIX[i] = sM->EMISSION_GAP_Y_MATRIX[i] * scale + shift;
        // level_stdev = stdev * var
        sM->EMISSION_GAP_Y_MATRIX[i + 1] = sM->EMISSION_GAP_Y_MATRIX[i + 1] * var;

        // Fluctuation (noise) adjustments
        // noise_mean *= scale_sd
        sM->EMISSION_MATCH_MATRIX[i + 2] = sM->EMISSION_MATCH_MATRIX[i + 2] * scale_sd;
        // noise_lambda *= var_sd
        sM->EMISSION_MATCH_MATRIX[i + 4] = sM->EMISSION_MATCH_MATRIX[i + 4] * var_sd;
        // noise_sd = sqrt(adjusted_noise_mean**3 / adjusted_noise_lambda);
        sM->EMISSION_MATCH_MATRIX[i + 3] = sqrt(pow(sM->EMISSION_MATCH_MATRIX[i + 2], 3.0)
                                                / sM->EMISSION_MATCH_MATRIX[i + 4]);

        // Fluctuation (noise) adjustments
        // noise_mean *= scale_sd
        sM->EMISSION_GAP_Y_MATRIX[i + 2] = sM->EMISSION_GAP_Y_MATRIX[i + 2] * scale_sd;
        // noise_lambda *= var_sd
        sM->EMISSION_GAP_Y_MATRIX[i + 4] = sM->EMISSION_GAP_Y_MATRIX[i + 4] * var_sd;
        // noise_sd = sqrt(adjusted_noise_mean**3 / adjusted_noise_lambda);
        sM->EMISSION_GAP_Y_MATRIX[i + 3] = sqrt(pow(sM->EMISSION_GAP_Y_MATRIX[i + 2], 3.0)
                                                / sM->EMISSION_MATCH_MATRIX[i + 4]);
    }
}

////////////////////////////
// EM emissions functions //
////////////////////////////
static void emissions_em_loadMatchProbs(double *emissionMatchProbs, Hmm *hmm, int64_t matchState) {
    //Load the matches
    HmmDiscrete *hmmD = (HmmDiscrete *)hmm;
    for(int64_t x = 0; x < hmm->parameterSetSize; x++) {
        for(int64_t y = 0; y < hmm->parameterSetSize; y++) {
            emissionMatchProbs[x * hmm->parameterSetSize + y] = log(hmmD->getEmissionExpFcn(hmm, matchState, x, y));
        }
    }
}

static void emissions_em_loadMatchProbsSymmetrically(double *emissionMatchProbs, Hmm *hmm, int64_t matchState) {
    //Load the matches
    HmmDiscrete *hmmD = (HmmDiscrete *)hmm;
    for(int64_t x = 0; x < hmm->parameterSetSize; x++) {
        emissionMatchProbs[x * hmm->parameterSetSize + x] = log(hmmD->getEmissionExpFcn(hmm, matchState, x, x));
        for(int64_t y=x+1; y<hmm->parameterSetSize; y++) {
            double d = log((hmmD->getEmissionExpFcn(hmm, matchState, x, y) +
                    hmmD->getEmissionExpFcn(hmm, matchState, y, x)) / 2.0);
            emissionMatchProbs[x * hmm->parameterSetSize + y] = d;
            emissionMatchProbs[y * hmm->parameterSetSize + x] = d;
        }
    }
}

static void emissions_em_collapseMatrixEmissions(Hmm *hmm, int64_t state, double *gapEmissions, bool collapseToX) {
    HmmDiscrete *hmmD = (HmmDiscrete *)hmm;
    for(int64_t x=0; x<hmm->parameterSetSize; x++) {
        for(int64_t y=0; y<hmm->parameterSetSize; y++) {
            gapEmissions[collapseToX ? x : y] += hmmD->getEmissionExpFcn(hmm, state, x, y);
        }
    }
}


static void emissions_em_loadGapProbs(double *emissionGapProbs, Hmm *hmm,
                                      int64_t *xGapStates, int64_t xGapStateNo,
                                      int64_t *yGapStates, int64_t yGapStateNo) {
    //Initialise to 0.0
    for(int64_t i=0; i < hmm->parameterSetSize; i++) {
        emissionGapProbs[i] = 0.0;
    }
    //Load the probs taking the average over all the gap states
    for(int64_t i=0; i < xGapStateNo; i++) {
        emissions_em_collapseMatrixEmissions(hmm, xGapStates[i], emissionGapProbs, 1);
    }
    for(int64_t i=0; i<yGapStateNo; i++) {
        emissions_em_collapseMatrixEmissions(hmm, yGapStates[i], emissionGapProbs, 0);
    }
    //Now normalise
    double total = 0.0;
    for(int64_t i=0; i < hmm->parameterSetSize; i++) {
        total += emissionGapProbs[i];
    }
    for(int64_t i=0; i< hmm->parameterSetSize; i++) {
        emissionGapProbs[i] = log(emissionGapProbs[i]/total);
    }
}


///////////////////////////////////
///////////////////////////////////
//Five state state-machine
///////////////////////////////////
///////////////////////////////////

/////////////////////////////////////////// STATIC FUNCTIONS ////////////////////////////////////////////////////////

static double stateMachine5_startStateProb(StateMachine *sM, int64_t state) {
    //Match state is like going to a match.
    state_check(sM, state);
    return state == match ? 0 : LOG_ZERO;
}

static double stateMachine5_raggedStartStateProb(StateMachine *sM, int64_t state) {
    state_check(sM, state);
    return (state == longGapX || state == longGapY) ? 0 : LOG_ZERO;
}

static double stateMachine5_endStateProb(StateMachine *sM, int64_t state) {
    //End state is like to going to a match
    StateMachine5 *sM5 = (StateMachine5 *) sM;
    state_check(sM, state);
    switch (state) {
    case match:
        return sM5->TRANSITION_MATCH_CONTINUE;
    case shortGapX:
        return sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X;
    case shortGapY:
        return sM5->TRANSITION_MATCH_FROM_SHORT_GAP_Y;
    case longGapX:
        return sM5->TRANSITION_MATCH_FROM_LONG_GAP_X;
    case longGapY:
        return sM5->TRANSITION_MATCH_FROM_LONG_GAP_Y;
    }
    return 0.0;
}

static double stateMachine5_raggedEndStateProb(StateMachine *sM, int64_t state) {
    StateMachine5 *sM5 = (StateMachine5 *) sM;
    state_check(sM, state);
    switch (state) {
    case match:
        return sM5->TRANSITION_GAP_LONG_OPEN_X;
    case shortGapX:
        return sM5->TRANSITION_GAP_LONG_OPEN_X;
    case shortGapY:
        return sM5->TRANSITION_GAP_LONG_OPEN_Y;
    case longGapX:
        return sM5->TRANSITION_GAP_LONG_EXTEND_X;
    case longGapY:
        return sM5->TRANSITION_GAP_LONG_EXTEND_Y;
    }
    return 0.0;
}

static void stateMachine5_cellCalculate(StateMachine *sM,
                                        void *current, void *lower, void *middle, void *upper,
                                        void *cX, void *cY,
                                        void (*doTransition)(double *, double *, // fromCells, toCells
                                                             int64_t, int64_t,   // from, to
                                                             double, double,     // emissionProb, transitionProb
                                                             void *),            // extraArgs
                                        void *extraArgs) {
    st_errAbort("5-state stateMachine not implemented\n");

    StateMachine5 *sM5 = (StateMachine5 *) sM;
    if (lower != NULL) {
        double eP = sM5->getXGapProbFcn(sM5->model.EMISSION_GAP_X_PROBS, cX);
        doTransition(lower, current, match, shortGapX, eP, sM5->TRANSITION_GAP_SHORT_OPEN_X, extraArgs);
        doTransition(lower, current, shortGapX, shortGapX, eP, sM5->TRANSITION_GAP_SHORT_EXTEND_X, extraArgs);
        // how come these are commented out?
        //doTransition(lower, current, shortGapY, shortGapX, eP, sM5->TRANSITION_GAP_SHORT_SWITCH_TO_X, extraArgs);
        doTransition(lower, current, match, longGapX, eP, sM5->TRANSITION_GAP_LONG_OPEN_X, extraArgs);
        doTransition(lower, current, longGapX, longGapX, eP, sM5->TRANSITION_GAP_LONG_EXTEND_X, extraArgs);
        //doTransition(lower, current, longGapY, longGapX, eP, sM5->TRANSITION_GAP_LONG_SWITCH_TO_X, extraArgs);
    }
    if (middle != NULL) {
        double eP = sM5->getMatchProbFcn(sM5->model.EMISSION_MATCH_MATRIX, cX, cY);
        doTransition(middle, current, match, match, eP, sM5->TRANSITION_MATCH_CONTINUE, extraArgs);
        doTransition(middle, current, shortGapX, match, eP, sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X, extraArgs);
        doTransition(middle, current, shortGapY, match, eP, sM5->TRANSITION_MATCH_FROM_SHORT_GAP_Y, extraArgs);
        doTransition(middle, current, longGapX, match, eP, sM5->TRANSITION_MATCH_FROM_LONG_GAP_X, extraArgs);
        doTransition(middle, current, longGapY, match, eP, sM5->TRANSITION_MATCH_FROM_LONG_GAP_Y, extraArgs);
    }
    if (upper != NULL) {
        double eP = sM5->getYGapProbFcn(sM5->model.EMISSION_GAP_Y_MATRIX, cY);
        doTransition(upper, current, match, shortGapY, eP, sM5->TRANSITION_GAP_SHORT_OPEN_Y, extraArgs);
        doTransition(upper, current, shortGapY, shortGapY, eP, sM5->TRANSITION_GAP_SHORT_EXTEND_Y, extraArgs);
        //doTransition(upper, current, shortGapX, shortGapY, eP, sM5->TRANSITION_GAP_SHORT_SWITCH_TO_Y, extraArgs);
        doTransition(upper, current, match, longGapY, eP, sM5->TRANSITION_GAP_LONG_OPEN_Y, extraArgs);
        doTransition(upper, current, longGapY, longGapY, eP, sM5->TRANSITION_GAP_LONG_EXTEND_Y, extraArgs);
        //doTransition(upper, current, longGapX, longGapY, eP, sM5->TRANSITION_GAP_LONG_SWITCH_TO_Y, extraArgs);
    }
}

///////////////////////////////////////////// CORE FUNCTIONS ////////////////////////////////////////////////////////

StateMachine *stateMachine5_construct(StateMachineType type, int64_t parameterSetSize,
                                      void (*setEmissionsDefaults)(StateMachine *sM),
                                      double (*gapXProbFcn)(const double *, void *),
                                      double (*gapYProbFcn)(const double *, void *),
                                      double (*matchProbFcn)(const double *, void *, void *),
                                      void (*cellCalcUpdateExpFcn)(double *fromCells, double *toCells,
                                                                   int64_t from, int64_t to,
                                                                   double eP, double tP, void *extraArgs)) {
    /*
     * Description of (potentially ambigious) arguments:
     * parameterSetSize = the number of kmers that we are using, of len(kmer) = 1, then the number is 4 (or 5 if we're
     * including N). It's 25 if len(kmer) = 2, it's 4096 in the 6-mer model.
     *
     */
    st_errAbort("5-state stateMachine not implemented");
    StateMachine5 *sM5 = st_malloc(sizeof(StateMachine5));
    if(type != fiveState && type != fiveStateAsymmetric) {
        st_errAbort("Wrong type for five state %i", type);
    }
    // setup transitions, specific to stateMachine5
    sM5->TRANSITION_MATCH_CONTINUE = -0.030064059121770816; //0.9703833696510062f
    sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X = -1.272871422049609; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    sM5->TRANSITION_MATCH_FROM_LONG_GAP_X = -5.673280173170473; //1.0 - gapExtend = 0.00343657420938
    sM5->TRANSITION_GAP_SHORT_OPEN_X = -4.34381910900448; //0.0129868352330243
    sM5->TRANSITION_GAP_SHORT_EXTEND_X = -0.3388262689231553; //0.7126062401851738f;
    sM5->TRANSITION_GAP_SHORT_SWITCH_TO_X = -4.910694825551255; //0.0073673675173412815f;
    sM5->TRANSITION_GAP_LONG_OPEN_X = -6.30810595366929; //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    sM5->TRANSITION_GAP_LONG_EXTEND_X = -0.003442492794189331; //0.99656342579062f;
    sM5->TRANSITION_GAP_LONG_SWITCH_TO_X = -6.30810595366929; //0.99656342579062f;
    // make it symmetric
    sM5->TRANSITION_MATCH_FROM_SHORT_GAP_Y = sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X;
    sM5->TRANSITION_MATCH_FROM_LONG_GAP_Y = sM5->TRANSITION_MATCH_FROM_LONG_GAP_X;
    sM5->TRANSITION_GAP_SHORT_OPEN_Y = sM5->TRANSITION_GAP_SHORT_OPEN_X;
    sM5->TRANSITION_GAP_SHORT_EXTEND_Y = sM5->TRANSITION_GAP_SHORT_EXTEND_X;
    sM5->TRANSITION_GAP_SHORT_SWITCH_TO_Y = sM5->TRANSITION_GAP_SHORT_SWITCH_TO_X;
    sM5->TRANSITION_GAP_LONG_OPEN_Y = sM5->TRANSITION_GAP_LONG_OPEN_X;
    sM5->TRANSITION_GAP_LONG_EXTEND_Y = sM5->TRANSITION_GAP_LONG_EXTEND_X;
    sM5->TRANSITION_GAP_LONG_SWITCH_TO_Y = sM5->TRANSITION_GAP_LONG_SWITCH_TO_X;
    // setup the parent class
    sM5->model.type = type;
    sM5->model.parameterSetSize = parameterSetSize;
    sM5->model.stateNumber = 5;
    sM5->model.matchState = match;
    sM5->model.startStateProb = stateMachine5_startStateProb;
    sM5->model.endStateProb = stateMachine5_endStateProb;
    sM5->model.raggedStartStateProb = stateMachine5_raggedStartStateProb;
    sM5->model.raggedEndStateProb = stateMachine5_raggedEndStateProb;
    sM5->model.cellCalculate = stateMachine5_cellCalculate;
    sM5->model.cellCalculateUpdateExpectations = cellCalcUpdateExpFcn;

    sM5->getXGapProbFcn = gapXProbFcn;
    sM5->getYGapProbFcn = gapYProbFcn;
    sM5->getMatchProbFcn = matchProbFcn;

    // set emissions to defaults (or zeros)
    setEmissionsDefaults((StateMachine *) sM5);

    return (StateMachine *) sM5;
}

static void switchDoubles(double *a, double *b) {
    double c = *a;
    *a = *b;
    *b = c;
}

///////////////////////////////
// EM - StateMachine5 functions/
///////////////////////////////

static void stateMachine5_loadAsymmetric(StateMachine5 *sM5, Hmm *hmm) {
    if (hmm->type != fiveStateAsymmetric) {
        st_errAbort("Wrong hmm type");
    }
    sM5->TRANSITION_MATCH_CONTINUE = log(hmm->getTransitionsExpFcn(hmm, match, match)); //0.9703833696510062f

    sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X = log(hmm->getTransitionsExpFcn(hmm, shortGapX, match));
    sM5->TRANSITION_MATCH_FROM_LONG_GAP_X = log(hmm->getTransitionsExpFcn(hmm, longGapX, match));
    sM5->TRANSITION_GAP_SHORT_OPEN_X = log(hmm->getTransitionsExpFcn(hmm, match, shortGapX));
    sM5->TRANSITION_GAP_SHORT_EXTEND_X = log(hmm->getTransitionsExpFcn(hmm, shortGapX, shortGapX));
    sM5->TRANSITION_GAP_SHORT_SWITCH_TO_X = log(hmm->getTransitionsExpFcn(hmm, shortGapY, shortGapX));
    sM5->TRANSITION_GAP_LONG_OPEN_X = log(hmm->getTransitionsExpFcn(hmm, match, longGapX));
    sM5->TRANSITION_GAP_LONG_EXTEND_X = log(hmm->getTransitionsExpFcn(hmm, longGapX, longGapX));
    sM5->TRANSITION_GAP_LONG_SWITCH_TO_X = log(hmm->getTransitionsExpFcn(hmm, longGapY, longGapX));

    if(sM5->TRANSITION_GAP_SHORT_EXTEND_X > sM5->TRANSITION_GAP_LONG_EXTEND_X) {
        // Switch the long and short gap parameters if one the "long states" have a smaller
        // extend probability than the "short states", as can randomly happen during EM training.
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_EXTEND_X), &(sM5->TRANSITION_GAP_LONG_EXTEND_X));
        switchDoubles(&(sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X), &(sM5->TRANSITION_MATCH_FROM_LONG_GAP_X));
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_OPEN_X), &(sM5->TRANSITION_GAP_LONG_OPEN_X));
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_SWITCH_TO_X), &(sM5->TRANSITION_GAP_LONG_SWITCH_TO_X));
    }

    sM5->TRANSITION_MATCH_FROM_SHORT_GAP_Y = log(hmm->getTransitionsExpFcn(hmm, shortGapY, match));
    sM5->TRANSITION_MATCH_FROM_LONG_GAP_Y = log(hmm->getTransitionsExpFcn(hmm, longGapY, match));
    sM5->TRANSITION_GAP_SHORT_OPEN_Y = log(hmm->getTransitionsExpFcn(hmm, match, shortGapY));
    sM5->TRANSITION_GAP_SHORT_EXTEND_Y = log(hmm->getTransitionsExpFcn(hmm, shortGapY, shortGapY));
    sM5->TRANSITION_GAP_SHORT_SWITCH_TO_Y = log(hmm->getTransitionsExpFcn(hmm, shortGapX, shortGapY));
    sM5->TRANSITION_GAP_LONG_OPEN_Y = log(hmm->getTransitionsExpFcn(hmm, match, longGapY));
    sM5->TRANSITION_GAP_LONG_EXTEND_Y = log(hmm->getTransitionsExpFcn(hmm, longGapY, longGapY));
    sM5->TRANSITION_GAP_LONG_SWITCH_TO_Y = log(hmm->getTransitionsExpFcn(hmm, longGapX, longGapY));

    if(sM5->TRANSITION_GAP_SHORT_EXTEND_Y > sM5->TRANSITION_GAP_LONG_EXTEND_Y) {
        // Switch the long and short gap parameters if one the "long states" have a smaller
        // extend probability than the "short states", as can randomly happen during EM training.
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_EXTEND_Y), &(sM5->TRANSITION_GAP_LONG_EXTEND_Y));
        switchDoubles(&(sM5->TRANSITION_MATCH_FROM_SHORT_GAP_Y), &(sM5->TRANSITION_MATCH_FROM_LONG_GAP_Y));
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_OPEN_Y), &(sM5->TRANSITION_GAP_LONG_OPEN_Y));
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_SWITCH_TO_Y), &(sM5->TRANSITION_GAP_LONG_SWITCH_TO_Y));
    }

    emissions_em_loadMatchProbs(sM5->model.EMISSION_MATCH_MATRIX, hmm, match);
    int64_t xGapStates[2] = { shortGapX, longGapX };
    int64_t yGapStates[2] = { shortGapY, longGapY };
    emissions_em_loadGapProbs(sM5->model.EMISSION_GAP_X_PROBS, hmm, xGapStates, 2, NULL, 0);
    emissions_em_loadGapProbs(sM5->model.EMISSION_GAP_Y_MATRIX, hmm, NULL, 0, yGapStates, 2);
}

static void stateMachine5_loadSymmetric(StateMachine5 *sM5, Hmm *hmm) {
    if (hmm->type != fiveState) {
        printf("Wrong hmm type");
        st_errAbort("Wrong hmm type");
    }

    sM5->TRANSITION_MATCH_CONTINUE = log(hmm->getTransitionsExpFcn(hmm, match, match));
    sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X = log(
            (hmm->getTransitionsExpFcn(hmm, shortGapX, match) +
             hmm->getTransitionsExpFcn(hmm, shortGapY, match)) / 2); //1.0 - gapExtend - gapSwitch = 0.280026392297485
    sM5->TRANSITION_MATCH_FROM_LONG_GAP_X = log(
            (hmm->getTransitionsExpFcn(hmm, longGapX, match) +
             hmm->getTransitionsExpFcn(hmm, longGapY, match)) / 2); //1.0 - gapExtend = 0.00343657420938
    sM5->TRANSITION_GAP_SHORT_OPEN_X = log(
            (hmm->getTransitionsExpFcn(hmm, match, shortGapX) +
             hmm->getTransitionsExpFcn(hmm, match, shortGapY)) / 2); //0.0129868352330243
    sM5->TRANSITION_GAP_SHORT_EXTEND_X = log(
            (hmm->getTransitionsExpFcn(hmm, shortGapX, shortGapX) +
             hmm->getTransitionsExpFcn(hmm, shortGapY, shortGapY)) / 2); //0.7126062401851738f;
    sM5->TRANSITION_GAP_SHORT_SWITCH_TO_X = log(
            (hmm->getTransitionsExpFcn(hmm, shortGapX, shortGapY) +
             hmm->getTransitionsExpFcn(hmm, shortGapY, shortGapX)) / 2); //0.0073673675173412815f;
    sM5->TRANSITION_GAP_LONG_OPEN_X = log(
            (hmm->getTransitionsExpFcn(hmm, match, longGapX) +
             hmm->getTransitionsExpFcn(hmm, match, longGapY)) / 2); //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    sM5->TRANSITION_GAP_LONG_EXTEND_X = log(
            (hmm->getTransitionsExpFcn(hmm, longGapX, longGapX) +
             hmm->getTransitionsExpFcn(hmm, longGapY, longGapY)) / 2);
    sM5->TRANSITION_GAP_LONG_SWITCH_TO_X = log(
            (hmm->getTransitionsExpFcn(hmm, longGapX, longGapY) +
             hmm->getTransitionsExpFcn(hmm, longGapY, longGapX)) / 2); //0.0073673675173412815f;

    if(sM5->TRANSITION_GAP_SHORT_EXTEND_X > sM5->TRANSITION_GAP_LONG_EXTEND_X) {
        //Switch the long and short gap parameters if one the "long states" have a smaller extend probability than the "short states", as can randomly happen during EM training.
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_EXTEND_X), &(sM5->TRANSITION_GAP_LONG_EXTEND_X));
        switchDoubles(&(sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X), &(sM5->TRANSITION_MATCH_FROM_LONG_GAP_X));
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_OPEN_X), &(sM5->TRANSITION_GAP_LONG_OPEN_X));
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_SWITCH_TO_X), &(sM5->TRANSITION_GAP_LONG_SWITCH_TO_X));
    }

    sM5->TRANSITION_MATCH_FROM_SHORT_GAP_Y = sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X;
    sM5->TRANSITION_MATCH_FROM_LONG_GAP_Y = sM5->TRANSITION_MATCH_FROM_LONG_GAP_X;
    sM5->TRANSITION_GAP_SHORT_OPEN_Y = sM5->TRANSITION_GAP_SHORT_OPEN_X;
    sM5->TRANSITION_GAP_SHORT_EXTEND_Y = sM5->TRANSITION_GAP_SHORT_EXTEND_X;
    sM5->TRANSITION_GAP_SHORT_SWITCH_TO_Y = sM5->TRANSITION_GAP_SHORT_SWITCH_TO_X;
    sM5->TRANSITION_GAP_LONG_OPEN_Y = sM5->TRANSITION_GAP_LONG_OPEN_X;
    sM5->TRANSITION_GAP_LONG_EXTEND_Y = sM5->TRANSITION_GAP_LONG_EXTEND_X;
    sM5->TRANSITION_GAP_LONG_SWITCH_TO_Y = sM5->TRANSITION_GAP_LONG_SWITCH_TO_X;

    emissions_em_loadMatchProbsSymmetrically(sM5->model.EMISSION_MATCH_MATRIX, hmm, match);
    int64_t xGapStates[2] = { shortGapX, longGapX };
    int64_t yGapStates[2] = { shortGapY, longGapY };
    emissions_em_loadGapProbs(sM5->model.EMISSION_GAP_X_PROBS, hmm, xGapStates, 2, yGapStates, 2);
    emissions_em_loadGapProbs(sM5->model.EMISSION_GAP_Y_MATRIX, hmm, xGapStates, 2, yGapStates, 2);
}


//////////////////////////////////////////////////////////////////////
//Three state state-machine [StateMachine3 and StateMachine3Vanilla]
//////////////////////////////////////////////////////////////////////

/////////////////////////////////////////// STATIC FUNCTIONS ////////////////////////////////////////////////////////

// Transitions //
typedef enum {
    match0 = 0, match1 = 1, match2 = 2, match3 = 3, match4 = 4, match5 = 5, gapX = 6
} SignalState;

static bool stateMachine3_checkStateNumber(int64_t stateNumber, StateMachineType type) {
    if (((type == threeState) || (type == threeStateHdp)) && (stateNumber == 3)) {
        return TRUE;
    }
    if ((type == fiveState) && (stateNumber == 5)) {
        return TRUE;
    }
    return FALSE;
}

static double stateMachine3_startStateProb(StateMachine *sM, int64_t state) {
    //Match state is like going to a match.
    state_check(sM, state);
    return state == match ? 0 : LOG_ZERO;
}

static double stateMachine3_raggedStartStateProb(StateMachine *sM, int64_t state) {
    state_check(sM, state);
    return (state == shortGapX || state == shortGapY) ? 0 : LOG_ZERO;
}

static double stateMachine3_endStateProb(StateMachine *sM, int64_t state) {
    //End state is like to going to a match
    StateMachine3 *sM3 = (StateMachine3 *) sM;
    state_check(sM, state);
    switch (state) {
    case match:
        return sM3->TRANSITION_MATCH_CONTINUE;
    case shortGapX:
        return sM3->TRANSITION_MATCH_FROM_GAP_X;
    case shortGapY:
        return sM3->TRANSITION_MATCH_FROM_GAP_Y;
    }
    return 0.0;
}

static double stateMachine3_raggedEndStateProb(StateMachine *sM, int64_t state) {
    //End state is like to going to a match
    StateMachine3 *sM3 = (StateMachine3 *) sM;
    state_check(sM, state);
    switch (state) {
    case match:
        return (sM3->TRANSITION_GAP_OPEN_X + sM3->TRANSITION_GAP_OPEN_Y) / 2.0;
    case shortGapX:
        return sM3->TRANSITION_GAP_EXTEND_X;
    case shortGapY:
        return sM3->TRANSITION_GAP_EXTEND_Y;
    }
    return 0.0;
}

void stateMachine3_setTransitionsToNucleotideDefaults(StateMachine *sM) {
    StateMachine3 *sM3 = (StateMachine3 *) sM;
    sM3->TRANSITION_MATCH_CONTINUE = -0.030064059121770816; //0.9703833696510062f
    sM3->TRANSITION_MATCH_FROM_GAP_X = -1.272871422049609; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    sM3->TRANSITION_MATCH_FROM_GAP_Y = -1.272871422049609; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    sM3->TRANSITION_GAP_OPEN_X = -4.21256642; //0.0129868352330243
    sM3->TRANSITION_GAP_OPEN_Y = -4.21256642; //0.0129868352330243
    sM3->TRANSITION_GAP_EXTEND_X = -0.3388262689231553; //0.7126062401851738f;
    sM3->TRANSITION_GAP_EXTEND_Y = -0.3388262689231553; //0.7126062401851738f;
    sM3->TRANSITION_GAP_SWITCH_TO_X = -4.910694825551255; //0.0073673675173412815f;
    sM3->TRANSITION_GAP_SWITCH_TO_Y = -4.910694825551255; //0.0073673675173412815f;
}


void stateMachine3_setTransitionsToNanoporeDefaults(StateMachine *sM) {
    StateMachine3 *sM3 = (StateMachine3 *) sM;
    sM3->TRANSITION_MATCH_CONTINUE = -0.23552123624314988; // log(step_prob) (0.79015888282447311)
    sM3->TRANSITION_MATCH_FROM_GAP_X = -0.21880828092192281; // log(1 - skip_prob) (1 - 0.19652425498269727)
    sM3->TRANSITION_MATCH_FROM_GAP_Y = -0.013406326748077823; // log(1 - stay_prob) (0.98668313780708949)
    sM3->TRANSITION_GAP_OPEN_X = -1.6269694202638481; // log(skip_prob) (0.19652425498269727)
    sM3->TRANSITION_GAP_OPEN_Y = -4.3187242127300092; // log(1 - (skip_prob + step_prob)) (0.013316862192829682)
    sM3->TRANSITION_GAP_EXTEND_X = -1.6269694202638481; // log(skip_prob) (0.19652425498269727)
    sM3->TRANSITION_GAP_EXTEND_Y = -4.3187242127239411; // log(stay_prob) 0.013316862192910478
    sM3->TRANSITION_GAP_SWITCH_TO_X = LOG_ZERO;
    sM3->TRANSITION_GAP_SWITCH_TO_Y = LOG_ZERO;
}

static void stateMachine3_loadTransitionsFromFile(StateMachine *sM, stList *transitions) {
    StateMachine3 *self = (StateMachine3 *)sM;
    int64_t j;
    double transition;

    // match->match
    j = sscanf(stList_get(transitions, 0), "%lf", &transition);
    if (j != 1) {
        st_errAbort("stateMachine3_loadFromFile: error parsing match->match transition\n");
    }
    self->TRANSITION_MATCH_CONTINUE = log(transition);
    // match->gapX
    j = sscanf(stList_get(transitions, 1), "%lf", &transition);
    if (j != 1) {
        st_errAbort("stateMachine3_loadFromFile: error parsing match->gapX transition\n");
    }
    self->TRANSITION_GAP_OPEN_X = log(transition);
    // match->gapY
    j = sscanf(stList_get(transitions, 2), "%lf", &transition);
    if (j != 1) {
        st_errAbort("stateMachine3_loadFromFile: error parsing match->gapY transition\n");
    }
    self->TRANSITION_GAP_OPEN_Y = log(transition);

    // gapX->match
    j = sscanf(stList_get(transitions, 3), "%lf", &transition);
    if (j != 1) {
        st_errAbort("stateMachine3_loadFromFile: error parsing gapX->match transition\n");
    }
    self->TRANSITION_MATCH_FROM_GAP_X = log(transition);
    // gapX->gapX
    j = sscanf(stList_get(transitions, 4), "%lf", &transition);
    if (j != 1) {
        st_errAbort("stateMachine3_loadFromFile: error parsing gapX->gapX transition\n");
    }
    self->TRANSITION_GAP_EXTEND_X = log(transition);
    // gapX->gapY skip, set to log zero

    // gapY->match
    j = sscanf(stList_get(transitions, 6), "%lf", &transition);
    if (j != 1) {
        st_errAbort("stateMachine3_loadFromFile: error parsing gapY->match transition\n");
    }
    self->TRANSITION_MATCH_FROM_GAP_Y = log(transition);
    // gapY->gapX
    j = sscanf(stList_get(transitions, 7), "%lf", &transition);
    if (j != 1) {
        st_errAbort("stateMachine3_loadFromFile: error parsing gapY->gapX transition\n");
    }
    self->TRANSITION_GAP_SWITCH_TO_Y = log(transition);
    // gapY->gapY
    j = sscanf(stList_get(transitions, 8), "%lf", &transition);
    if (j != 1) {
        st_errAbort("stateMachine3_loadFromFile: error parsing gapY->gapY transition\n");
    }
    self->TRANSITION_GAP_EXTEND_Y = log(transition);
}
/* Temp change for r9 experiment

void stateMachine3_setTransitionsToNanoporeDefaults(StateMachine *sM) {
    StateMachine3 *sM3 = (StateMachine3 *) sM;
    sM3->TRANSITION_MATCH_CONTINUE = -0.43078291609245423; // log(step_prob) (0.65)
    sM3->TRANSITION_MATCH_FROM_GAP_X = -0.41551544396166595; // log(1 - skip_prob) (1 - 0.34)
    sM3->TRANSITION_MATCH_FROM_GAP_Y = -0.010050335853501451; // log(1 - stay_prob) (1 - 0.01)
    sM3->TRANSITION_GAP_OPEN_X = -1.0788096613719298; // log(skip_prob) (0.34)
    sM3->TRANSITION_GAP_OPEN_Y = -4.6051701859880909; // log(1 - (skip_prob + step_prob)) (1 - (0.34 + 0.65))
    sM3->TRANSITION_GAP_EXTEND_X = -1.0788096613719298; // log(skip_prob) (0.34)
    sM3->TRANSITION_GAP_EXTEND_Y = -4.6051701859880909; // log(stay_prob) 0.013316862192910478
    sM3->TRANSITION_GAP_SWITCH_TO_X = LOG_ZERO;
    sM3->TRANSITION_GAP_SWITCH_TO_Y = LOG_ZERO;
}
*/

void stateMachine3_setModelToHdpExpectedValues(StateMachine *sM, NanoporeHDP *nhdp) {
    // get all the kmers that the HDP has data for
    if (strcmp(sM->alphabet, nhdp->alphabet) != 0) {
        st_errAbort("stateMachine3_setModelToHdpExpectedValues: "
                            "This Nanopore Hdp and stateMachine have different alphabets\n");
    }
    if (sM->alphabetSize != nhdp->alphabet_size) {
        st_errAbort("stateMachine3_setModelToHdpExpectedValues: StateMachine alphabet size and Nanopore HDP"
                            "alphabet size aren't the same");
    }
    if (sM->kmerLength != nhdp->kmer_length) {
        st_errAbort("stateMachine3_setModelToHdpExpectedValues: StateMachine kmer length is not the same as "
                            "NanoporeHdp kmer length");
    }

    stList *kmers = path_listPotentialKmers(nhdp->kmer_length, nhdp->alphabet_size, nhdp->alphabet);
    for (int64_t i = 0; i < stList_length(kmers); i++) {
        char *kmer = (char *)stList_get(kmers, i);
        int64_t kmerIndex = kmer_id(kmer, sM->alphabet, sM->alphabetSize, sM->kmerLength);
        int64_t meanIndex = kmerIndex * MODEL_PARAMS;
        int64_t sdIndex = kmerIndex * MODEL_PARAMS + 1;
        bool isObs = hdp_check_for_observed(nhdp->hdp, kmerIndex);
        if (isObs) {
            double newMean = dir_proc_expected_val(nhdp->hdp, kmerIndex);
            double newSd = sqrt(dir_proc_variance(nhdp->hdp, kmerIndex));
            sM->EMISSION_MATCH_MATRIX[meanIndex] = newMean;
            sM->EMISSION_MATCH_MATRIX[sdIndex] = newSd;
        }
    }
}

static void stateMachine3_cellCalculate(StateMachine *sM,
                                        void *current, void *lower, void *middle, void *upper,
                                        void *cX, void *cY,
                                        void (*doTransition)(double *, double *,
                                                             int64_t, int64_t,
                                                             double, double,
                                                             void *),
                                        void *extraArgs) {
    StateMachine3 *sM3 = (StateMachine3 *)sM;
    HDCell *hdCurrent = current == NULL ? NULL : (HDCell *)current;
    HDCell *hdLower = lower == NULL ? NULL : (HDCell *)lower;
    HDCell *hdMiddle = middle == NULL ? NULL : (HDCell *)middle;
    HDCell *hdUpper = upper == NULL ? NULL : (HDCell *)upper;

    if (hdLower != NULL) {
        for (int64_t p = 0; p < hdCurrent->numberOfPaths; p++) {
            Path *pathC = hdCell_getPath(hdCurrent, p);
            for (int64_t q = 0; q < hdLower->numberOfPaths; q++) {
                Path *pathL = hdCell_getPath(hdLower, q);
                if (path_checkLegal(pathL, pathC)) {
                    double *lowerCells = path_getCell(pathL);
                    double *currentCells = path_getCell(pathC);
                    double eP = sM3->getXGapProbFcn(sM, sM3->model.EMISSION_GAP_X_PROBS, pathC->kmer);
                    doTransition(lowerCells, currentCells, match, shortGapX, eP, sM3->TRANSITION_GAP_OPEN_X, extraArgs);
                    doTransition(lowerCells, currentCells, shortGapX, shortGapX, eP, sM3->TRANSITION_GAP_EXTEND_X, extraArgs);
                    doTransition(lowerCells, currentCells, shortGapY, shortGapX, eP, sM3->TRANSITION_GAP_SWITCH_TO_X, extraArgs);
                }
            }
        }
    }
    if (hdMiddle != NULL) {
        for (int64_t p = 0; p < hdCurrent->numberOfPaths; p++) {
            Path *pathC = hdCell_getPath(hdCurrent, p);
            for (int64_t q = 0; q < hdMiddle->numberOfPaths; q++) {
                Path *pathM = hdCell_getPath(hdMiddle, q);
                if (path_checkLegal(pathM, pathC)) {
                    double *middleCells = path_getCell(pathM);
                    double *currentCells = path_getCell(pathC);
                    double eP = sM3->getMatchProbFcn(sM, pathC->kmer, cY, TRUE);
                    doTransition(middleCells, currentCells, match, match, eP, sM3->TRANSITION_MATCH_CONTINUE, extraArgs);
                    doTransition(middleCells, currentCells, shortGapX, match, eP, sM3->TRANSITION_MATCH_FROM_GAP_X, extraArgs);
                    doTransition(middleCells, currentCells, shortGapY, match, eP, sM3->TRANSITION_MATCH_FROM_GAP_Y, extraArgs);
                }
            }
        }
    }
    if (hdUpper != NULL) {
        for (int64_t p = 0; p < hdCurrent->numberOfPaths; p++) {
            Path *pathC = hdCell_getPath(hdCurrent, p);
            for (int64_t q = 0; q < hdUpper->numberOfPaths; q++) {
                Path *pathU = hdCell_getPath(hdUpper, q);
                if (stString_eq(pathC->kmer, pathU->kmer)) {
                    double *upperCells = path_getCell(pathU);
                    double *currentCells = path_getCell(pathC);
                    double eP = sM3->getYGapProbFcn(sM, pathC->kmer, cY, FALSE);
                    doTransition(upperCells, currentCells, match, shortGapY, eP, sM3->TRANSITION_GAP_OPEN_Y, extraArgs);
                    doTransition(upperCells, currentCells, shortGapY, shortGapY, eP, sM3->TRANSITION_GAP_EXTEND_Y, extraArgs);
                    // shortGapX -> shortGapY not allowed, this would be going from a kmer skip to extra event?
                    //doTransition(upper, current, shortGapX, shortGapY, eP, sM3->TRANSITION_GAP_SWITCH_TO_Y, extraArgs);
                }
            }
        }
    }
}

static void stateMachine3HDP_cellCalculate(StateMachine *sM,
                                           void *current, void *lower, void *middle, void *upper,
                                           void *cX, void *cY,
                                           void (*doTransition)(double *, double *,
                                                                int64_t, int64_t,
                                                                double, double,
                                                                void *),
                                           void *extraArgs) {
    StateMachine3_HDP *sM3 = (StateMachine3_HDP *) sM;
    HDCell *hdCurrent = current == NULL ? NULL : (HDCell *)current;
    HDCell *hdLower = lower == NULL ? NULL : (HDCell *)lower;
    HDCell *hdMiddle = middle == NULL ? NULL : (HDCell *)middle;
    HDCell *hdUpper = upper == NULL ? NULL : (HDCell *)upper;

    if (hdLower != NULL) {
        for (int64_t p = 0; p < hdCurrent->numberOfPaths; p ++) {
            Path *pathC = hdCell_getPath(hdCurrent, p);
            for (int64_t q = 0; q < hdLower->numberOfPaths; q++) {
                Path *pathL = hdCell_getPath(hdLower, q);
                if (path_checkLegal(pathL, pathC)) {
                    //st_uglyf("SENTINAL - legal LOWER : pathC kmer %s\n", pathC->kmer);
                    double *lowerCells = path_getCell(pathL);
                    double *currentCells = path_getCell(pathC);
                    double eP = -2.3025850929940455; // log(0.1)
                    doTransition(lowerCells, currentCells, match, shortGapX, eP, sM3->TRANSITION_GAP_OPEN_X, extraArgs);
                    doTransition(lowerCells, currentCells, shortGapX, shortGapX, eP, sM3->TRANSITION_GAP_EXTEND_X, extraArgs);
                    doTransition(lowerCells, currentCells, shortGapY, shortGapX, eP, sM3->TRANSITION_GAP_SWITCH_TO_X, extraArgs);
                }
            }
        }
    }
    if (hdMiddle != NULL) {
        for (int64_t p = 0; p < hdCurrent->numberOfPaths; p++) {
            Path *pathC = hdCell_getPath(hdCurrent, p);
            for (int64_t q = 0; q < hdMiddle->numberOfPaths; q++) {
                Path *pathM = hdCell_getPath(hdMiddle, q);
                if (path_checkLegal(pathM, pathC)) {
                    //st_uglyf("SENTINAL - legal MIDDLE : pathC kmer %s\n", pathC->kmer);
                    double *middleCells = path_getCell(pathM);
                    double *curentCells = path_getCell(pathC);
                    double eP = sM3->getMatchProbFcn(sM, pathC->kmer, cY, TRUE);
                    doTransition(middleCells, curentCells, match, match, eP, sM3->TRANSITION_MATCH_CONTINUE, extraArgs);
                    doTransition(middleCells, curentCells, shortGapX, match, eP, sM3->TRANSITION_MATCH_FROM_GAP_X, extraArgs);
                    doTransition(middleCells, curentCells, shortGapY, match, eP, sM3->TRANSITION_MATCH_FROM_GAP_Y, extraArgs);
                }
            }
        }
    }
    if (hdUpper != NULL) {
        for (int64_t p = 0; p < hdCurrent->numberOfPaths; p++) {
            Path *pathC = hdCell_getPath(hdCurrent, p);
            for (int64_t q = 0; q < hdUpper->numberOfPaths; q++) {
                Path *pathU = hdCell_getPath(hdUpper, q);
                if (stString_eq(pathC->kmer, pathU->kmer)) {
                    //st_uglyf("SENTINAL - legal UPPER : pathC kmer %s\n", pathC->kmer);
                    double *upperCells = path_getCell(pathU);
                    double *currentCells = path_getCell(pathC);
                    double eP = sM3->getMatchProbFcn(sM, pathC->kmer, cY, FALSE);
                    doTransition(upperCells, currentCells, match, shortGapY, eP, sM3->TRANSITION_GAP_OPEN_Y, extraArgs);
                    doTransition(upperCells, currentCells, shortGapY, shortGapY, eP, sM3->TRANSITION_GAP_EXTEND_Y, extraArgs);
                    // shortGapX -> shortGapY not allowed, this would be going from a kmer skip to extra event?
                    //doTransition(upper, current, shortGapX, shortGapY, eP, sM3->TRANSITION_GAP_SWITCH_TO_Y, extraArgs);
                }
            }
        }
    }
}

///////////////////////////////////////////// CORE FUNCTIONS ////////////////////////////////////////////////////////
StateMachine *stateMachine3_loadFromFile(const char *modelFile, StateMachineType type,
                                         double (*gapXProbFcn)(StateMachine *, const double *, void *),
                                         double (*matchProbFcn)(StateMachine *, void *, void *, bool ),
                                         void (*loadTransitionsFcn)(StateMachine *, stList *),
                                         NanoporeHDP *nHdp) {
    /*
     *  the model file has the format:
     *  line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
     *  line 1: match->match \t match->gapX \t match->gapY \t
     *          gapX->match \t gapX->gapX \t gapX->gapY \t
     *          gapY->match \t gapY->gapX \t gapY->gapY \n
     *  line 1: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
     */
    if (!stFile_exists(modelFile)) {
        st_errAbort("stateMachine3_loadFromFile: Couldn't find .model file here %s\n", modelFile);
    }

    FILE *fH = fopen(modelFile, "r");

    // Line 0: parse the stateMachine stateNumber, alphabet, alphabet length, and kmer size parameteres
    char *string = stFile_getLineFromFile(fH);
    stList *tokens = stString_split(string);
    if (stList_length(tokens) != 4) {
        st_errAbort("stateMachine3_loadFromFile: Model file %s does not have the correct number of stateMachine "
                            "parameters should be 4 got %"PRId64"\n", modelFile, stList_length(tokens));
    }

    char *alphabet;
    int64_t stateNumber, alphabetSize, kmerLength, j, h;

    j = sscanf(stList_get(tokens, 0), "%"SCNd64, &stateNumber);
    if (j != 1) {
        st_errAbort("stateMachine3_loadFromFile: error parsing alphabet size\n");
    }
    j = sscanf(stList_get(tokens, 1), "%"SCNd64, &alphabetSize);
    if (j != 1) {
        st_errAbort("stateMachine3_loadFromFile: error parsing alphabet size\n");
    }
    alphabet = (char *)stList_get(tokens, 2);
    j = sscanf(stList_get(tokens, 3), "%"SCNd64, &kmerLength);
    if (j != 1) {
        st_errAbort("stateMachine3_loadFromFile: error parsing kmer length\n");
    }

    // check for correct number of states for given stateMachineType
    if (!stateMachine3_checkStateNumber(stateNumber, type)) {
        st_errAbort("stateMachine3_loadFromFile: Got invalid stateNumber for this stateMachine. Got stateNumber %s"
                            " and StateMachineType %i\n", stateNumber, type);
    };

    StateMachine *sM = stateMachine3_signalMachineBuilder(type, alphabet, alphabetSize, kmerLength, gapXProbFcn,
                                                          matchProbFcn, nHdp);
    free(string);
    stList_destruct(tokens);

    // line 1: Transitions line
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    if (stList_length(tokens) != (sM->stateNumber * sM->stateNumber + 1)) {
        st_errAbort("stateMachine3_loadFromFile: Got invalid number of tokens on transitions line, got %lld should be"
                            " %lld", stList_length(tokens), (sM->stateNumber * sM->stateNumber + 1));
    }

    // load the transitions
    loadTransitionsFcn(sM, tokens);
    free(string);
    stList_destruct(tokens);

    // Line 2: parse the match emissions line
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    // check to make sure that the model will fit in the stateMachine
    if (stList_length(tokens) != sM->parameterSetSize * MODEL_PARAMS) {
        st_errAbort("This stateMachine is not correct for signal model (match emissions) got %lld, should be %lld\n",
                    stList_length(tokens), sM->parameterSetSize * MODEL_PARAMS);
    }
    // load the model into the state machine emissions
    for (int64_t i = 0; i < sM->parameterSetSize * MODEL_PARAMS; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(sM->EMISSION_MATCH_MATRIX[i]));
        h = sscanf(stList_get(tokens, i), "%lf", &(sM->EMISSION_GAP_Y_MATRIX[i]));
        if ((j != 1) || (h != 1)) {
            st_errAbort("emissions_signal_loadPoreModel: error loading pore model (match emissions)\n");
        }
    }
    // clean up match emissions line
    free(string);
    stList_destruct(tokens);

    // increase the level_sd for the GapY state by 75% see Simpson et al.
    // start at 1 bc. the level_sd is the second param in MODEL PRAMS
    for (int64_t i = 1; i < (sM->parameterSetSize * MODEL_PARAMS); i += MODEL_PARAMS) {
        sM->EMISSION_GAP_Y_MATRIX[i] *= EXTRA_EVENT_NOISE_MULTIPLIER;
    }

    // close file
    fclose(fH);

    return sM;
}

StateMachine *stateMachine3_construct(StateMachineType type,
                                      const char *alphabet, int64_t alphabetSize, int64_t kmerLength,
                                      void (*setTransitionsToDefaults)(StateMachine *),
                                      void (*setEmissionsDefaults)(StateMachine *, int64_t),
                                      double (*gapXProbFcn)(StateMachine *, const double *, void *),
                                      double (*gapYProbFcn)(StateMachine *, void *, void *, bool ),
                                      double (*matchProbFcn)(StateMachine *, void *, void *, bool ),
                                      void (*cellCalcUpdateExpFcn)(double *fromCells, double *toCells,
                                                                   int64_t from, int64_t to,
                                                                   double eP, double tP, void *extraArgs)) {
    StateMachine3 *sM3 = st_malloc(sizeof(StateMachine3));
    if (type != threeState && type != threeStateAsymmetric) {
        st_errAbort("Tried to create a three state state-machine with the wrong type");
    }

    // setup the parent class
    sM3->model.type = type;
    sM3->model.parameterSetSize = intPow(alphabetSize, kmerLength);

    sM3->model.alphabetSize = alphabetSize;
    sM3->model.alphabet = sequence_prepareAlphabet(alphabet, alphabetSize);
    sM3->model.kmerLength = kmerLength;

    sM3->model.stateNumber = 3;
    sM3->model.matchState = match;
    sM3->model.startStateProb = stateMachine3_startStateProb;
    sM3->model.endStateProb = stateMachine3_endStateProb;
    sM3->model.raggedStartStateProb = stateMachine3_raggedStartStateProb;
    sM3->model.raggedEndStateProb = stateMachine3_raggedEndStateProb;
    sM3->model.cellCalculate = stateMachine3_cellCalculate;
    sM3->model.cellCalculateUpdateExpectations = cellCalcUpdateExpFcn;

    // setup functions
    sM3->getXGapProbFcn = gapXProbFcn;
    sM3->getYGapProbFcn = gapYProbFcn;
    sM3->getMatchProbFcn = matchProbFcn;

    // setup transitions
    setTransitionsToDefaults((StateMachine *) sM3);

    // set emissions to defaults or zeros
    setEmissionsDefaults((StateMachine *) sM3, sM3->model.parameterSetSize); // mallocs are in here

    // set gap probs
    for (int64_t i = 0; i < sM3->model.parameterSetSize; i++) {
        sM3->model.EMISSION_GAP_X_PROBS[i] = -2.3025850929940455; // log(0.1)
    }

    return (StateMachine *) sM3;
}

StateMachine *stateMachine3Hdp_construct(StateMachineType type,
                                         const char *alphabet, int64_t alphabetSize, int64_t kmerLength,
                                         void (*setTransitionsToDefaults)(StateMachine *),
                                         void (*setEmissionsDefaults)(StateMachine *, int64_t),
                                         NanoporeHDP *hdpModel,
                                         double (*matchProbFcn)(StateMachine *, void *, void *, bool),
                                         void (*cellCalcUpdateExpFcn)(double *, double *, int64_t, int64_t,
                                                                      double , double , void *)) {
    StateMachine3_HDP *sM3 = st_malloc(sizeof(StateMachine3_HDP));
    if (type != threeStateHdp) {
        st_errAbort("Tried to create a three state state-machine with the wrong type");
    }

    // setup the parent class
    sM3->model.type = type;
    //sM3->model.parameterSetSize = parameterSetSize;

    sM3->model.parameterSetSize = intPow(alphabetSize, kmerLength);
    sM3->model.alphabetSize = alphabetSize;
    sM3->model.alphabet = sequence_prepareAlphabet(alphabet, alphabetSize);
    sM3->model.kmerLength = kmerLength;

    sM3->model.stateNumber = 3;
    sM3->model.matchState = match;
    sM3->model.startStateProb = stateMachine3_startStateProb;
    sM3->model.endStateProb = stateMachine3_endStateProb;
    sM3->model.raggedStartStateProb = stateMachine3_raggedStartStateProb;
    sM3->model.raggedEndStateProb = stateMachine3_raggedEndStateProb;
    sM3->model.cellCalculate = stateMachine3HDP_cellCalculate;
    sM3->model.cellCalculateUpdateExpectations = cellCalcUpdateExpFcn;

    // setup functions
    sM3->getMatchProbFcn = matchProbFcn;

    // setup HDP
    sM3->hdpModel = hdpModel;

    // setup transitions
    setTransitionsToDefaults((StateMachine *) sM3);
    // set emissions to defaults or zeros
    setEmissionsDefaults((StateMachine *) sM3, sM3->model.parameterSetSize);

    // initialize kmer emission gap probs
    for (int64_t i = 0; i < sM3->model.parameterSetSize; i++) {
        sM3->model.EMISSION_GAP_X_PROBS[i] = -2.3025850929940455; // log(0.1)
    }
    return (StateMachine *) sM3;
}

StateMachine *stateMachine3_signalMachineBuilder(StateMachineType type, char *alphabet, int64_t alphabetSize,
                                                 int64_t kmerLength,
                                                 double (*gapXProbFcn)(StateMachine *, const double *, void *),
                                                 double (*matchProbFcn)(StateMachine *, void *, void *, bool),
                                                 NanoporeHDP *nHdp) {
    StateMachine *sM;
    switch (type) {
        case threeStateHdp:
            sM = stateMachine3Hdp_construct(threeStateHdp, alphabet, alphabetSize, kmerLength,
                                            stateMachine3_setTransitionsToNanoporeDefaults,
                                            emissions_signal_initEmissionsToZero,
                                            nHdp, matchProbFcn, cell_signal_updateExpectationsAndAssignments);
            return sM;
        default:
            sM = stateMachine3_construct(threeState, alphabet, alphabetSize, kmerLength,
                                         stateMachine3_setTransitionsToNanoporeDefaults,
                                         emissions_signal_initEmissionsToZero,
                                         gapXProbFcn, matchProbFcn, matchProbFcn,
                                         cell_signal_updateExpectations);
            return sM;
    }
}

/*
static void stateMachine3_loadAsymmetric(StateMachine3 *sM3, Hmm *hmm) {
    if (hmm->type != threeStateAsymmetric) {
        st_errAbort("Wrong hmm type");
    }
    //sM3->TRANSITION_MATCH_CONTINUE = log(hmm_getTransition(hmm, match, match));
    sM3->TRANSITION_MATCH_CONTINUE = log(hmm->getTransitionsExpFcn(hmm, match, match));
    sM3->TRANSITION_MATCH_FROM_GAP_X = log(hmm->getTransitionsExpFcn(hmm, shortGapX, match));
    sM3->TRANSITION_MATCH_FROM_GAP_Y = log(hmm->getTransitionsExpFcn(hmm, shortGapY, match));
    sM3->TRANSITION_GAP_OPEN_X = log(hmm->getTransitionsExpFcn(hmm, match, shortGapX));
    sM3->TRANSITION_GAP_OPEN_Y = log(hmm->getTransitionsExpFcn(hmm, match, shortGapY));
    sM3->TRANSITION_GAP_EXTEND_X = log(hmm->getTransitionsExpFcn(hmm, shortGapX, shortGapX));
    sM3->TRANSITION_GAP_EXTEND_Y = log(hmm->getTransitionsExpFcn(hmm, shortGapY, shortGapY));
    sM3->TRANSITION_GAP_SWITCH_TO_X = log(hmm->getTransitionsExpFcn(hmm, shortGapY, shortGapX));
    sM3->TRANSITION_GAP_SWITCH_TO_Y = log(hmm->getTransitionsExpFcn(hmm, shortGapX, shortGapY));
    emissions_em_loadMatchProbs(sM3->EMISSION_MATCH_MATRIX, hmm, match);
    int64_t xGapStates[1] = { shortGapX };
    int64_t yGapStates[1] = { shortGapY };
    emissions_loadGapProbs(sM3->EMISSION_GAP_X_PROBS, hmm, xGapStates, 1, NULL, 0);
    emissions_loadGapProbs(sM3->EMISSION_GAP_Y_MATRIX, hmm, NULL, 0, yGapStates, 1);
}

static void stateMachine3_loadSymmetric(StateMachine3 *sM3, Hmm *hmm) {
    if (hmm->type != threeState) {
        st_errAbort("Wrong hmm type");
    }
    sM3->TRANSITION_MATCH_CONTINUE = log(hmm_getTransition(hmm, match, match));
    sM3->TRANSITION_MATCH_FROM_GAP_X = log(
            (hmm_getTransition(hmm, shortGapX, match) + hmm_getTransition(hmm, shortGapY, match)) / 2.0);
    sM3->TRANSITION_MATCH_FROM_GAP_Y = sM3->TRANSITION_MATCH_FROM_GAP_X;
    sM3->TRANSITION_GAP_OPEN_X = log(
            (hmm_getTransition(hmm, match, shortGapX) + hmm_getTransition(hmm, match, shortGapY)) / 2.0);
    sM3->TRANSITION_GAP_OPEN_Y = sM3->TRANSITION_GAP_OPEN_X;
    sM3->TRANSITION_GAP_EXTEND_X = log(
            (hmm_getTransition(hmm, shortGapX, shortGapX) + hmm_getTransition(hmm, shortGapY, shortGapY)) / 2.0);
    sM3->TRANSITION_GAP_EXTEND_Y = sM3->TRANSITION_GAP_EXTEND_X;
    sM3->TRANSITION_GAP_SWITCH_TO_X = log(
            (hmm_getTransition(hmm, shortGapY, shortGapX) + hmm_getTransition(hmm, shortGapX, shortGapY)) / 2.0);
    sM3->TRANSITION_GAP_SWITCH_TO_Y = sM3->TRANSITION_GAP_SWITCH_TO_X;
    emissions_em_loadMatchProbsSymmetrically(sM3->EMISSION_MATCH_MATRIX, hmm, match);
    int64_t xGapStates[2] = { shortGapX };
    int64_t yGapStates[2] = { shortGapY };
    emissions_loadGapProbs(sM3->EMISSION_GAP_X_PROBS, hmm, xGapStates, 1, yGapStates, 1);
    emissions_em_loadGapProbs(sM3->EMISSION_GAP_Y_MATRIX, hmm, xGapStates, 1, yGapStates, 1);
}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Public functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
StateMachine *getStateMachine5(Hmm *hmmD, StateMachineFunctions *sMfs) {
    st_errAbort("5-state stateMachine not implemented\n");
    if (hmmD->type == fiveState) {
        StateMachine5 *sM5 = (StateMachine5 *) stateMachine5_construct(fiveState, hmmD->parameterSetSize,
                                                                       emissions_discrete_initEmissionsToZero,
                                                                       sMfs->gapXProbFcn,
                                                                       sMfs->gapYProbFcn,
                                                                       sMfs->matchProbFcn,
                                                                       cell_updateExpectations);
        stateMachine5_loadSymmetric(sM5, hmmD);
        return (StateMachine *) sM5;
    }
    if (hmmD->type == fiveStateAsymmetric) {
        StateMachine5 *sM5 = (StateMachine5 *) stateMachine5_construct(fiveState, hmmD->parameterSetSize,
                                                                       emissions_discrete_initEmissionsToZero,
                                                                       sMfs->gapXProbFcn,
                                                                       sMfs->gapYProbFcn,
                                                                       sMfs->matchProbFcn,
                                                                       cell_updateExpectations);
        stateMachine5_loadAsymmetric(sM5, hmmD);
        return (StateMachine *) sM5;
    }
    else {
        return 0;
    }
}

StateMachine *getStateMachine3_descaled(const char *modelFile, NanoporeReadAdjustmentParameters npp, bool scaleNoise) {
    if (!stFile_exists(modelFile)) {
        st_errAbort("getStateMachine3_descaled: Cannot find model file %s\n", modelFile);
    };
    StateMachine *sM = stateMachine3_loadFromFile(modelFile, threeState, emissions_kmer_getGapProb,
                                                  emissions_signal_strawManGetKmerEventMatchProbWithDescaling,
                                                  stateMachine3_loadTransitionsFromFile, NULL);

    sM->scale = npp.scale;
    sM->shift = npp.shift;
    sM->var = npp.var;

    if (scaleNoise) {
        emissions_signal_scaleNoise(sM, npp);
    }
    return sM;
}

StateMachine *getStateMachine3(const char *modelFile) {
    if (!stFile_exists(modelFile)) {
        st_errAbort("getStateMachine3: Cannot find model file %s\n", modelFile);
    };
    StateMachine *sM = stateMachine3_loadFromFile(modelFile, threeState, emissions_kmer_getGapProb,
                                                  emissions_signal_strawManGetKmerEventMatchProb,
                                                  stateMachine3_loadTransitionsFromFile, NULL);
    return sM;
}

StateMachine *getHdpStateMachine(NanoporeHDP *hdp, const char *modelFile, NanoporeReadAdjustmentParameters npp) {
    if (!stFile_exists(modelFile)) {
        st_errAbort("getHdpStateMachine: Cannot find model file %s\n", modelFile);
    };
    StateMachine *sM = stateMachine3_loadFromFile(modelFile, threeStateHdp, NULL, emissions_signal_getHdpKmerDensity,
                                                  stateMachine3_loadTransitionsFromFile, hdp);
    sM->scale = npp.scale;
    sM->shift = npp.shift;
    sM->var = npp.var;

    // set the level_mean table to the expected value from the HDP
    stateMachine3_setModelToHdpExpectedValues(sM, hdp);

    return sM;
}

static void stateMachine_destructHdpModel(StateMachine *sM) {
    StateMachine3_HDP *sMHdp = (StateMachine3_HDP *)sM;
    destroy_nanopore_hdp(sMHdp->hdpModel);
}

void stateMachine_destruct(StateMachine *sM) {
    if (sM->type == threeStateHdp) {
        stateMachine_destructHdpModel(sM);
    }
    free(sM->EMISSION_GAP_X_PROBS);
    free(sM->EMISSION_GAP_Y_MATRIX);
    free(sM->EMISSION_MATCH_MATRIX);
    free(sM->alphabet);
    free(sM);
}



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
#include <stdint.h>
#include <stateMachine.h>
#include "nanopore.h"
#include "nanopore_hdp.h"
#include "bioioC.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
#include "stateMachine.h"
#include "emissionMatrix.h"
#include "discreteHmm.h"

//////////////////////////////////////////////////////////////////////////////
// StateMachine Emission functions for discrete alignments (symbols and kmers)
//////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////// STATIC FUNCTIONS ////////////////////////////////////////////////////////
// moved to stateMachine.h 10/28
//typedef enum {
//    match = 0, shortGapX = 1, shortGapY = 2, longGapX = 3, longGapY = 4
//} State;

static void state_check(StateMachine *sM, State s) {
    assert(s >= 0 && s < sM->stateNumber);
}

static void index_check(int64_t c) {
    assert(c >= 0 && c < NUM_OF_KMERS);
}

static inline void emissions_discrete_initializeEmissionsMatrices(StateMachine *sM) {
    sM->EMISSION_GAP_X_PROBS = st_malloc(sM->parameterSetSize*sizeof(double));
    sM->EMISSION_GAP_Y_PROBS = st_malloc(sM->parameterSetSize*sizeof(double));
    sM->EMISSION_MATCH_PROBS = st_malloc(sM->parameterSetSize*sM->parameterSetSize*sizeof(double));
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

    memcpy(sM->EMISSION_MATCH_PROBS, M, sizeof(double)*SYMBOL_NUMBER_NO_N*SYMBOL_NUMBER_NO_N);

    // Set Gap probs to default values
    const double EMISSION_GAP = -1.6094379124341003; //log(0.2)
    const double G[4] = { EMISSION_GAP, EMISSION_GAP, EMISSION_GAP, EMISSION_GAP };
    memcpy(sM->EMISSION_GAP_X_PROBS, G, sizeof(double)*SYMBOL_NUMBER_NO_N);
    memcpy(sM->EMISSION_GAP_Y_PROBS, G, sizeof(double)*SYMBOL_NUMBER_NO_N);
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
    emissions_discrete_initMatchProbsToZero(sM->EMISSION_MATCH_PROBS, sM->parameterSetSize);
    // set gap matrix to zeros
    emissions_discrete_initGapProbsToZero(sM->EMISSION_GAP_X_PROBS, sM->parameterSetSize);
    emissions_discrete_initGapProbsToZero(sM->EMISSION_GAP_Y_PROBS, sM->parameterSetSize);
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

int64_t emissions_discrete_getKmerIndex(void *kmer) {
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

int64_t emissions_discrete_getKmerIndexFromKmer(void *kmer) {
    // make temp kmer meant to work with getKmer
    char *kmer_i = malloc((KMER_LENGTH) * sizeof(char));
    for (int64_t x = 0; x < KMER_LENGTH; x++) {
        kmer_i[x] = *((char *)kmer+x);
    }
    kmer_i[KMER_LENGTH] = '\0';

    int64_t i = emissions_discrete_getKmerIndex(kmer_i);
    //index_check(i);
    free(kmer_i);
    return i;
}

double emissions_symbol_getGapProb(const double *emissionGapProbs, void *base) {
    int64_t i = emissions_discrete_getBaseIndex(base);
    index_check(i);
    if(i == 4) {
        return -1.386294361; //log(0.25)
    }
    return emissionGapProbs[i];
}

double emissions_symbol_getMatchProb(const double *emissionMatchProbs, void *x, void *y) {
    int64_t iX = emissions_discrete_getBaseIndex(x);
    int64_t iY = emissions_discrete_getBaseIndex(y);
    index_check(iX);
    index_check(iY);
    if(iX == 4 || iY == 4) {
        return -2.772588722; //log(0.25**2)
    }
    return emissionMatchProbs[iX * SYMBOL_NUMBER_NO_N + iY];
}

double emissions_kmer_getGapProb(const double *emissionGapProbs, void *x_i) {
    // even though this shouldn't happen, check for null x_i and return logzero
    if (x_i == NULL) {
        return LOG_ZERO;
    }
    // make temp x_i meant to work with getKmer
    char *kmer_i = malloc((KMER_LENGTH) * sizeof(char));
    for (int64_t x = 0; x < KMER_LENGTH; x++) {
        kmer_i[x] = *((char *) x_i + x);
    }
    kmer_i[KMER_LENGTH] = '\0';
    int64_t i = emissions_discrete_getKmerIndex(kmer_i);
    //index_check(i);
    free(kmer_i);
    double p = i > NUM_OF_KMERS ? LOG_ZERO : emissionGapProbs[i];
    return p;
}

double emissions_kmer_getMatchProb(const double *emissionMatchProbs, void *x, void *y) {
    int64_t iX = emissions_discrete_getKmerIndex(x);
    int64_t iY = emissions_discrete_getKmerIndex(y);
    int64_t tableIndex = iX * NUM_OF_KMERS + iY;
    return emissionMatchProbs[tableIndex];
}
/////////////////////////////////////////
// functions for signal/kmer alignment //
/////////////////////////////////////////

/////////////////////////////////////////// STATIC FUNCTIONS ////////////////////////////////////////////////////////

static inline void emissions_vanilla_initializeEmissionsMatrices(StateMachine *sM, int64_t nbSkipParams) {
    // changed to 30 for skip prob bins
    // the kmer/gap and skip (f(|ui-1 - ui|)) probs have smaller tables, either 30 for the skip or parameterSetSize
    // for the naive kmer skip one
    // see note at emissions_signal_initEmissionsToZero about this change
    sM->EMISSION_GAP_X_PROBS = st_malloc(nbSkipParams * sizeof(double));

    // both the Iy and M - type states use the event/kmer match model so the matrices need to be the same size
    sM->EMISSION_GAP_Y_PROBS = st_malloc(1 + (sM->parameterSetSize * MODEL_PARAMS) * sizeof(double));
    sM->EMISSION_MATCH_PROBS = st_malloc(1 + (sM->parameterSetSize * MODEL_PARAMS) * sizeof(double));
}

static inline void emissions_signal_initMatchMatrixToZero(double *matchModel, int64_t parameterSetSize) {
    memset(matchModel, 0, (1 + (parameterSetSize * MODEL_PARAMS)) * sizeof(double));
}

static inline void emissions_signal_initKmerSkipTableToZero(double *skipModel, int64_t parameterSetSize) {
    memset(skipModel, 0, parameterSetSize * sizeof(double));
}

static inline double emissions_signal_getModelLevelMean(const double *eventModel, int64_t kmerIndex) {
    // 1 + i*MODEL PARAMS because the first element is the correlation parameter
    return kmerIndex > NUM_OF_KMERS ? 0.0 : eventModel[1 + (kmerIndex * MODEL_PARAMS)];
}

static inline double emissions_signal_getModelLevelSd(const double *eventModel, int64_t kmerIndex) {
    return kmerIndex > NUM_OF_KMERS ? 0.0 : eventModel[1 + (kmerIndex * MODEL_PARAMS + 1)];
}

static inline double emissions_signal_getModelFluctuationMean(const double *eventModel, int64_t kmerIndex) {
    return kmerIndex > NUM_OF_KMERS ? 0.0 : eventModel[1 + (kmerIndex * MODEL_PARAMS + 2)];
}

static inline double emissions_signal_getModelFluctuationSd(const double *eventModel, int64_t kmerIndex) {
    return kmerIndex > NUM_OF_KMERS ? 0.0 : eventModel[1 + (kmerIndex * MODEL_PARAMS + 3)];
}

static inline double emissions_signal_getModelFluctuationLambda(const double *eventModel, int64_t kmerIndex) {
    return kmerIndex > NUM_OF_KMERS ? 0.0 : eventModel[1 + (kmerIndex * MODEL_PARAMS + 4)];
}

static void emissions_signal_loadPoreModel(StateMachine *sM, const char *modelFile, StateMachineType type) {
    /*
     *  the model file has the format:
     *  line 1: [correlation coefficient] [level_mean] [level_sd] [noise_mean]
     *              [noise_sd] [noise_lambda ](.../kmer) \n
     *  line 2: [correlation coefficient] [level_mean] [level_sd, scaled]
     *              [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
     */

    FILE *fH = fopen(modelFile, "r");

    // Line 1: parse the match emissions line
    char *string = stFile_getLineFromFile(fH);
    stList *tokens = stString_split(string);
    // check to make sure that the model will fit in the stateMachine
    if (stList_length(tokens) != 1 + (sM->parameterSetSize * MODEL_PARAMS)) {
        st_errAbort("This stateMachine is not correct for signal model (match emissions) got %lld, should be %llf\n",
                    stList_length(tokens), 1 + (sM->parameterSetSize * MODEL_PARAMS));
    }
    // load the model into the state machine emissions
    for (int64_t i = 0; i < 1 + (sM->parameterSetSize * MODEL_PARAMS); i++) {
        int64_t j = sscanf(stList_get(tokens, i), "%lf", &(sM->EMISSION_MATCH_PROBS[i]));
        if (j != 1) {
            st_errAbort("emissions_signal_loadPoreModel: error loading pore model (match emissions)\n");
        }
    }
    // clean up match emissions line
    free(string);
    stList_destruct(tokens);

    // parse Y Gap emissions line
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    // check to make sure that the model will fit in the stateMachine
    if (stList_length(tokens) != 1 + (sM->parameterSetSize * MODEL_PARAMS)) {
        st_errAbort("This stateMachine is not correct for signal model (dupeEvent - Y emissions)\n");
    }
    // load the model into the state machine emissions
    for (int64_t i = 0; i < 1 + (sM->parameterSetSize * MODEL_PARAMS); i++) {
        int64_t j = sscanf(stList_get(tokens, i), "%lf", &(sM->EMISSION_GAP_Y_PROBS[i]));
        if (j != 1) {
            st_errAbort("emissions_signal_loadPoreModel: error loading pore model (dupeEvent - Y emissions)\n");
        }
    }
    // clean up match emissions line
    free(string);
    stList_destruct(tokens);

    // close file
    fclose(fH);
}

static inline double emissions_signal_logInvGaussPdf(double eventNoise, double modelNoiseMean,
                                                     double modelNoiseLambda) {
    double l_twoPi = 1.8378770664093453;// log(2*pi)
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

void emissions_signal_initEmissionsToZero(StateMachine *sM, int64_t nbSkipParams) {
    // initialize
    emissions_vanilla_initializeEmissionsMatrices(sM, nbSkipParams);

    // set kmer skip matrix to zeros
    emissions_signal_initKmerSkipTableToZero(sM->EMISSION_GAP_X_PROBS, nbSkipParams);

    // set extra event matrix to zeros
    emissions_signal_initMatchMatrixToZero(sM->EMISSION_GAP_Y_PROBS, sM->parameterSetSize);

    // set match matrix to zeros
    emissions_signal_initMatchMatrixToZero(sM->EMISSION_MATCH_PROBS, sM->parameterSetSize);
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
    int64_t k_i= emissions_discrete_getKmerIndex(kmer_i);
    int64_t k_im1 = emissions_discrete_getKmerIndex(kmer_im1);

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
    int64_t bin = emissions_signal_getKmerSkipBin(sM->EMISSION_MATCH_PROBS, kmers);
    return getAlpha ? sM->EMISSION_GAP_X_PROBS[bin+30] : sM->EMISSION_GAP_X_PROBS[bin];
}

double emissions_signal_getKmerSkipProb(StateMachine *sM, void *kmers) {
    StateMachine3Vanilla *sM3v = (StateMachine3Vanilla *) sM;
    // make kmer_i-1
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
    int64_t k_i= emissions_discrete_getKmerIndex(kmer_i);
    int64_t k_im1 = emissions_discrete_getKmerIndex(kmer_im1);

    // get the expected mean current for each one
    double u_ki = emissions_signal_getModelLevelMean(sM3v->model.EMISSION_MATCH_PROBS, k_i);
    double u_kim1 = emissions_signal_getModelLevelMean(sM3v->model.EMISSION_MATCH_PROBS, k_im1);

    // find the difference
    double d = fabs(u_ki - u_kim1);

    // get the 'bin' for skip prob, clamp to the last bin
    int64_t bin = (int64_t)(d / 0.5); // 0.5 pA bins right now
    bin = bin >= 30 ? 29 : bin;

    // for debugging
    //double pSkip = sM3v->model.EMISSION_GAP_X_PROBS[bin];
    //st_uglyf("SENTINAL - getKmerSkipProb; kmer_i: %s (index: %lld, model u=%f), kmer_im1: %s (index: %lld, model u=%f), d: %f, bin: %lld, pSkip: %f\n",
    //         kmer_i, k_i, u_ki, kmer_im1, k_im1, u_kim1, d, bin, pSkip);

    // clean up
    free(kmer_im1);
    free(kmer_i);
    // NOT log space
    return sM3v->model.EMISSION_GAP_X_PROBS[bin];
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
    int64_t kmerIndex = emissions_discrete_getKmerIndex(kmer_i);
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
    int64_t kmerIndex = emissions_discrete_getKmerIndex(kmer_i);

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

    int64_t kmerIndex = emissions_discrete_getKmerIndex(kmer_i);

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

double emissions_signal_strawManGetKmerEventMatchProb(const double *eventModel, void *x_i, void *e_j) {
    if (x_i == NULL) {
        return LOG_ZERO;
    }
    // this is meant to work with getKmer (NOT getKmer2)
    // wrangle e_j data
    double eventMean = *(double *) e_j;
    double eventNoise = *(double *) ((char *) e_j + sizeof(double)); // aaah pointers

    // make temp x_i
    char *kmer_i = malloc((KMER_LENGTH) * sizeof(char));
    for (int64_t x = 0; x < KMER_LENGTH; x++) {
        kmer_i[x] = *((char *) x_i + x);
    }
    kmer_i[KMER_LENGTH] = '\0';

    // get index
    int64_t kmerIndex = emissions_discrete_getKmerIndex(kmer_i);

    // get the µ and σ for the level and noise for the model
    double levelMean = emissions_signal_getModelLevelMean(eventModel, kmerIndex);
    double levelStdDev = emissions_signal_getModelLevelSd(eventModel, kmerIndex);
    double noiseMean = emissions_signal_getModelFluctuationMean(eventModel, kmerIndex);
    double noiseStdDev = emissions_signal_getModelFluctuationSd(eventModel, kmerIndex);

    double l_probEventMean = emissions_signal_logGaussPdf(eventMean, levelMean, levelStdDev);
    double l_probEventNoise = emissions_signal_logGaussPdf(eventNoise, noiseMean, noiseStdDev);

    // clean
    free(kmer_i);

    // debugging
    //double prob = l_probEventMean + l_probEventNoise;
    //st_uglyf("MATCHING--x_i:%s (index: %lld), e_j mean: %f, \n modelMean: %f, modelLsd: %f probEvent: %f probNoise: %f, combined: %f\n",
    //         kmer_i, kmerIndex, eventMean, levelMean, levelStdDev, l_probEventMean, l_probEventNoise, prob);

    return l_probEventMean + l_probEventNoise;
}

void emissions_signal_scaleModel(StateMachine *sM,
                                 double scale, double shift, double var,
                                 double scale_sd, double var_sd) {
    // model is arranged: level_mean, level_stdev, sd_mean, sd_stdev, sd_lambda per kmer
    // already been adjusted for correlation coeff.
    for (int64_t i = 1; i < (sM->parameterSetSize * MODEL_PARAMS) + 1; i += MODEL_PARAMS) {
        // Level adjustments
        // level_mean = mean * scale + shift
        sM->EMISSION_MATCH_PROBS[i] = sM->EMISSION_MATCH_PROBS[i] * scale + shift;
        // level_stdev = stdev * var
        sM->EMISSION_MATCH_PROBS[i+1] = sM->EMISSION_MATCH_PROBS[i+1] * var;

        // Fluctuation (noise) adjustments
        // noise_mean *= scale_sd
        sM->EMISSION_MATCH_PROBS[i+2] = sM->EMISSION_MATCH_PROBS[i+2] * scale_sd;
        // noise_lambda *= var_sd
        sM->EMISSION_MATCH_PROBS[i+4] = sM->EMISSION_MATCH_PROBS[i+4] * var_sd;
        // noise_sd = sqrt(adjusted_noise_mean**3 / adjusted_noise_lambda);
        sM->EMISSION_MATCH_PROBS[i+3] = sqrt(pow(sM->EMISSION_MATCH_PROBS[i+2], 3.0) / sM->EMISSION_MATCH_PROBS[i+4]);
    }
}

void emissions_signal_scaleModelNoiseOnly(StateMachine *sM,
                                          double scale, double shift, double var,
                                          double scale_sd, double var_sd) {
    // model is arranged: level_mean, level_stdev, sd_mean, sd_stdev, sd_lambda per kmer
    // already been adjusted for correlation coeff.
    for (int64_t i = 1; i < (sM->parameterSetSize * MODEL_PARAMS) + 1; i += MODEL_PARAMS) {
        // Level adjustments
        // level_mean = mean * scale + shift
        //sM->EMISSION_MATCH_PROBS[i] = sM->EMISSION_MATCH_PROBS[i] * scale + shift;
        // level_stdev = stdev * var
        sM->EMISSION_MATCH_PROBS[i+1] = sM->EMISSION_MATCH_PROBS[i+1] * var;

        // Fluctuation (noise) adjustments
        // noise_mean *= scale_sd
        sM->EMISSION_MATCH_PROBS[i+2] = sM->EMISSION_MATCH_PROBS[i+2] * scale_sd;
        // noise_lambda *= var_sd
        sM->EMISSION_MATCH_PROBS[i+4] = sM->EMISSION_MATCH_PROBS[i+4] * var_sd;
        // noise_sd = sqrt(adjusted_noise_mean**3 / adjusted_noise_lambda);
        sM->EMISSION_MATCH_PROBS[i+3] = sqrt(pow(sM->EMISSION_MATCH_PROBS[i+2], 3.0) / sM->EMISSION_MATCH_PROBS[i+4]);
    }
}

////////////////////////////
// EM emissions functions //
////////////////////////////

static void emissions_em_loadMatchProbs(double *emissionMatchProbs, Hmm *hmm, int64_t matchState) {
    //Load the matches
    for(int64_t x = 0; x < hmm->symbolSetSize; x++) {
        for(int64_t y = 0; y < hmm->symbolSetSize; y++) {
            emissionMatchProbs[x * hmm->symbolSetSize + y] = log(hmm->getEmissionExpFcn(hmm, matchState, x, y));
        }
    }
}

static void emissions_em_loadMatchProbsSymmetrically(double *emissionMatchProbs, Hmm *hmm, int64_t matchState) {
    //Load the matches
    for(int64_t x = 0; x < hmm->symbolSetSize; x++) {
        emissionMatchProbs[x * hmm->symbolSetSize + x] = log(hmm->getEmissionExpFcn(hmm, matchState, x, x));
        for(int64_t y=x+1; y<hmm->symbolSetSize; y++) {
            double d = log((hmm->getEmissionExpFcn(hmm, matchState, x, y) +
                    hmm->getEmissionExpFcn(hmm, matchState, y, x)) / 2.0);
            emissionMatchProbs[x * hmm->symbolSetSize + y] = d;
            emissionMatchProbs[y * hmm->symbolSetSize + x] = d;
        }
    }
}

static void emissions_em_collapseMatrixEmissions(Hmm *hmm, int64_t state, double *gapEmissions, bool collapseToX) {
    for(int64_t x=0; x<hmm->symbolSetSize; x++) {
        for(int64_t y=0; y<hmm->symbolSetSize; y++) {
            gapEmissions[collapseToX ? x : y] += hmm->getEmissionExpFcn(hmm, state, x, y);
        }
    }
}


static void emissions_em_loadGapProbs(double *emissionGapProbs, Hmm *hmm,
                                      int64_t *xGapStates, int64_t xGapStateNo,
                                      int64_t *yGapStates, int64_t yGapStateNo) {
    //Initialise to 0.0
    for(int64_t i=0; i < hmm->symbolSetSize; i++) {
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
    for(int64_t i=0; i < hmm->symbolSetSize; i++) {
        total += emissionGapProbs[i];
    }
    for(int64_t i=0; i< hmm->symbolSetSize; i++) {
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

static double stateMachine4_raggedStartStateProb(StateMachine *sM, int64_t state) {
    state_check(sM, state);
    return (state == longGapX || state == shortGapY) ? 0 : LOG_ZERO;
}

static double stateMachine4_endStateProb(StateMachine *sM, int64_t state) {
    StateMachine4 *sM4 = (StateMachine4 *)sM;
    state_check(sM, state);
    switch (state) {
        case match:
            return sM4->TRANSITION_MATCH_CONTINUE;
        case shortGapX:
            return sM4->TRANSITION_MATCH_FROM_SHORT_GAP_X;
        case shortGapY:
            return sM4->TRANSITION_MATCH_FROM_SHORT_GAP_Y;
        case longGapX:
            return sM4->TRANSITION_MATCH_FROM_LONG_GAP_X;
    }
    return 0.0;
}

static double stateMachine4_raggedEndStateProb(StateMachine *sM, int64_t state) {
    StateMachine4 *sM4 = (StateMachine4 *)sM;
    state_check(sM, state);
    switch (state) {
        case match:
            return sM4->TRANSITION_GAP_LONG_OPEN_X;
        case shortGapX:
            return sM4->TRANSITION_GAP_LONG_OPEN_X;
        case shortGapY:
            return sM4->TRANSITION_GAP_LONG_OPEN_X;
        case longGapX:
            return sM4->TRANSITION_GAP_LONG_EXTEND_X;
    }
    return 0.0;
}


static void stateMachine5_cellCalculate(StateMachine *sM,
                                        double *current, double *lower, double *middle, double *upper,
                                        void *cX, void *cY,
                                        void (*doTransition)(double *, double *, // fromCells, toCells
                                                             int64_t, int64_t,   // from, to
                                                             double, double,     // emissionProb, transitionProb
                                                             void *),            // extraArgs
                                        void *extraArgs) {
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
        double eP = sM5->getMatchProbFcn(sM5->model.EMISSION_MATCH_PROBS, cX, cY);
        doTransition(middle, current, match, match, eP, sM5->TRANSITION_MATCH_CONTINUE, extraArgs);
        doTransition(middle, current, shortGapX, match, eP, sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X, extraArgs);
        doTransition(middle, current, shortGapY, match, eP, sM5->TRANSITION_MATCH_FROM_SHORT_GAP_Y, extraArgs);
        doTransition(middle, current, longGapX, match, eP, sM5->TRANSITION_MATCH_FROM_LONG_GAP_X, extraArgs);
        doTransition(middle, current, longGapY, match, eP, sM5->TRANSITION_MATCH_FROM_LONG_GAP_Y, extraArgs);
    }
    if (upper != NULL) {
        double eP = sM5->getYGapProbFcn(sM5->model.EMISSION_GAP_Y_PROBS, cY);
        doTransition(upper, current, match, shortGapY, eP, sM5->TRANSITION_GAP_SHORT_OPEN_Y, extraArgs);
        doTransition(upper, current, shortGapY, shortGapY, eP, sM5->TRANSITION_GAP_SHORT_EXTEND_Y, extraArgs);
        //doTransition(upper, current, shortGapX, shortGapY, eP, sM5->TRANSITION_GAP_SHORT_SWITCH_TO_Y, extraArgs);
        doTransition(upper, current, match, longGapY, eP, sM5->TRANSITION_GAP_LONG_OPEN_Y, extraArgs);
        doTransition(upper, current, longGapY, longGapY, eP, sM5->TRANSITION_GAP_LONG_EXTEND_Y, extraArgs);
        //doTransition(upper, current, longGapX, longGapY, eP, sM5->TRANSITION_GAP_LONG_SWITCH_TO_Y, extraArgs);
    }
}

static void stateMachine4_cellCalculate(StateMachine *sM,
                                              double *current, double *lower, double *middle, double *upper,
                                            void *cX, void *cY,
                                            void (*doTransition)(double *, double *, // fromCells, toCells
                                                                 int64_t, int64_t,   // from, to
                                                                 double, double,     // emissionProb, transitionProb
                                                                 void *),            // extraArgs
                                            void *extraArgs) {
    StateMachine4 *sM4 = (StateMachine4 *) sM;
    if (lower != NULL) {
        double eP = sM4->getXGapProbFcn(sM4->model.EMISSION_GAP_X_PROBS, cX);
        doTransition(lower, current, match, shortGapX, eP, sM4->TRANSITION_GAP_SHORT_OPEN_X, extraArgs);
        doTransition(lower, current, shortGapX, shortGapX, eP, sM4->TRANSITION_GAP_SHORT_EXTEND_X, extraArgs);
        doTransition(lower, current, match, longGapX, eP, sM4->TRANSITION_GAP_LONG_OPEN_X, extraArgs);
        doTransition(lower, current, longGapX, longGapX, eP, sM4->TRANSITION_GAP_LONG_EXTEND_X, extraArgs);
        doTransition(lower, current, shortGapY, longGapX, eP, sM4->TRANSITION_GAP_LONG_SWITCH_TO_X, extraArgs);

    }
    if (middle != NULL) {
        double eP = sM4->getMatchProbFcn(sM4->model.EMISSION_MATCH_PROBS, cX, cY);
        doTransition(middle, current, match, match, eP, sM4->TRANSITION_MATCH_CONTINUE, extraArgs);
        doTransition(middle, current, shortGapX, match, eP, sM4->TRANSITION_MATCH_FROM_SHORT_GAP_X, extraArgs);
        doTransition(middle, current, shortGapY, match, eP, sM4->TRANSITION_MATCH_FROM_SHORT_GAP_Y, extraArgs);
        doTransition(middle, current, longGapX, match, eP, sM4->TRANSITION_MATCH_FROM_LONG_GAP_X, extraArgs);
    }
    if (upper != NULL) {
        double eP = sM4->getYGapProbFcn(sM4->model.EMISSION_GAP_Y_PROBS, cX, cY);
        doTransition(upper, current, match, shortGapY, eP, sM4->TRANSITION_GAP_SHORT_OPEN_Y, extraArgs);
        doTransition(upper, current, shortGapY, shortGapY, eP, sM4->TRANSITION_GAP_SHORT_EXTEND_Y, extraArgs);
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

StateMachine *stateMachine4_construct(StateMachineType type, int64_t parameterSetSize,
                                      void (*setEmissionsToDefaults)(StateMachine *, int64_t nbSkipParams),
                                      double (*gapXProbFcn)(const double *, void *),
                                      double (*gapYProbFcn)(const double *, void *, void *),
                                      double (*matchProbFcn)(const double *, void *, void *),
                                      void (*cellCalcUpdateFcn)(double *, double *, int64_t from, int64_t to,
                                                                double eP, double tP, void *)) {
    StateMachine4 *sM4 = st_malloc(sizeof(StateMachine4));
    if (type != fourState) {
        st_errAbort("Tried to make four-state stateMachine, with wrong type\n");
    }

    // setup parent class
    sM4->model.type = type,
    sM4->model.parameterSetSize = parameterSetSize;
    sM4->model.stateNumber = 4;
    sM4->model.matchState = match;
    // start
    sM4->model.startStateProb = stateMachine5_startStateProb;
    sM4->model.raggedStartStateProb = stateMachine4_raggedStartStateProb;
    // end
    sM4->model.endStateProb = stateMachine4_endStateProb;
    sM4->model.raggedEndStateProb = stateMachine4_raggedEndStateProb;
    sM4->model.cellCalculate = stateMachine4_cellCalculate;
    // cell calculate
    sM4->model.cellCalculateUpdateExpectations = cellCalcUpdateFcn;

    // setup functions
    sM4->getXGapProbFcn = gapXProbFcn;
    sM4->getYGapProbFcn = gapYProbFcn;
    sM4->getMatchProbFcn = matchProbFcn;


    // set transitions to defaults (these are from a template read)
    // from Match
    sM4->TRANSITION_MATCH_CONTINUE = -0.23552123624314988; // log(p_step) (0.79015888282447311)
    sM4->TRANSITION_GAP_SHORT_OPEN_X = -1.6269694202638481; // log(p_skip) 0.19652425498269727
    sM4->TRANSITION_GAP_SHORT_OPEN_Y = -4.7241893208381773; // log(1 - p_step - p_skip - p_longGapXOpen)
    sM4->TRANSITION_GAP_LONG_OPEN_X = -5.4173365013981227; // log((1 - p_step - p_skip)/3)

    // from shortGapX
    sM4->TRANSITION_GAP_SHORT_EXTEND_X = -1.6269694202638481; // log(p_skip)
    sM4->TRANSITION_MATCH_FROM_SHORT_GAP_X = -0.21880828092192281; // log(1 - skip_prob) (1 - 0.19652425498269727)

    // from longGapX
    sM4->TRANSITION_GAP_LONG_EXTEND_X = -0.003442492794189331; // 0.99656342579062f;
    sM4->TRANSITION_MATCH_FROM_LONG_GAP_X = -5.6732801731704612; // 1 - long gap extend

    // from shortGapY
    sM4->TRANSITION_MATCH_FROM_SHORT_GAP_Y = -0.013406326748077823; // log(1 - stay_prob) (0.98668313780708949)
    sM4->TRANSITION_GAP_SHORT_EXTEND_Y = -4.724189320832104; // log(0.008877908128607004)
    sM4->TRANSITION_GAP_LONG_SWITCH_TO_X = -5.4173365013920494; // log((1 - 0.98668313780708949)/3)

    /*
    // set transitions to defaults (these are from a complement read)
    // from Match
    sM4->TRANSITION_MATCH_CONTINUE = -0.38537441356664359; // log(p_step) (0.68019591394367218)
    sM4->TRANSITION_GAP_SHORT_OPEN_X = -1.1462174196841781; // log(p_skip) 0.31783674145479074
    sM4->TRANSITION_GAP_SHORT_OPEN_Y = -6.6365356716005186; // log(1 - p_step - p_skip - p_longGapXOpen)
    sM4->TRANSITION_GAP_LONG_OPEN_X = -7.329682852160464; // log((1 - p_step - p_skip)/3)

    // from shortGapX
    sM4->TRANSITION_GAP_SHORT_EXTEND_X = -1.1462174196841781; // log(p_skip)
    sM4->TRANSITION_MATCH_FROM_SHORT_GAP_X = -0.38248626775488298; // log(1 - skip_prob) (1 - 0.19652425498269727)

    // from longGapX
    sM4->TRANSITION_GAP_LONG_EXTEND_X = -0.003442492794189331; // 0.99656342579062f;
    sM4->TRANSITION_MATCH_FROM_LONG_GAP_X = -5.6732801731704612; // 1 - long gap extend

    // from shortGapY
    sM4->TRANSITION_MATCH_FROM_SHORT_GAP_Y = -0.0019692823655876384; // log(1 - stay_prob) (0.98668313780708949)
    sM4->TRANSITION_GAP_SHORT_EXTEND_Y = -6.6365356717310187; // log(0.008877908128607004)
    sM4->TRANSITION_GAP_LONG_SWITCH_TO_X = -7.3296828522909641; // log((1 - (p_matchFromShortGapY)/3)
    */

    // initialize emissions
    setEmissionsToDefaults((StateMachine *)sM4, parameterSetSize);

    return (StateMachine *)sM4;
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

    emissions_em_loadMatchProbs(sM5->model.EMISSION_MATCH_PROBS, hmm, match);
    int64_t xGapStates[2] = { shortGapX, longGapX };
    int64_t yGapStates[2] = { shortGapY, longGapY };
    emissions_em_loadGapProbs(sM5->model.EMISSION_GAP_X_PROBS, hmm, xGapStates, 2, NULL, 0);
    emissions_em_loadGapProbs(sM5->model.EMISSION_GAP_Y_PROBS, hmm, NULL, 0, yGapStates, 2);
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

    emissions_em_loadMatchProbsSymmetrically(sM5->model.EMISSION_MATCH_PROBS, hmm, match);
    int64_t xGapStates[2] = { shortGapX, longGapX };
    int64_t yGapStates[2] = { shortGapY, longGapY };
    emissions_em_loadGapProbs(sM5->model.EMISSION_GAP_X_PROBS, hmm, xGapStates, 2, yGapStates, 2);
    emissions_em_loadGapProbs(sM5->model.EMISSION_GAP_Y_PROBS, hmm, xGapStates, 2, yGapStates, 2);
}


//////////////////////////////////////////////////////////////////////
//Three state state-machine [StateMachine3 and StateMachine3Vanilla]
//////////////////////////////////////////////////////////////////////

/////////////////////////////////////////// STATIC FUNCTIONS ////////////////////////////////////////////////////////

// Transitions //
typedef enum {
    match0 = 0, match1 = 1, match2 = 2, match3 = 3, match4 = 4, match5 = 5, gapX = 6
} SignalState;

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

static double stateMachine3Vanilla_raggedEndStateProb(StateMachine *sM, int64_t state) {
    StateMachine3Vanilla *sM3v = (StateMachine3Vanilla *) sM;
    state_check(sM, state);
    switch (state) {
        case match:
            return (sM3v->DEFAULT_END_FROM_X_PROB + sM3v->DEFAULT_END_FROM_Y_PROB) / 2.0;
        case shortGapX:
            return sM3v->DEFAULT_END_FROM_X_PROB;
        case shortGapY:
            return sM3v->DEFAULT_END_FROM_Y_PROB;
    }
    return 0.0;
}

static double stateMachine3Vanilla_endStateProb(StateMachine *sM, int64_t state) {
    StateMachine3Vanilla *sM3v = (StateMachine3Vanilla *) sM;
    state_check(sM, state);
    switch (state) {
        case match:
            return sM3v->DEFAULT_END_MATCH_PROB;
        case shortGapX:
            return sM3v->DEFAULT_END_FROM_X_PROB;
        case shortGapY:
            return sM3v->DEFAULT_END_FROM_Y_PROB;
    }
    return 0.0;
}

static double stateMachineEchelon_startStateProb(StateMachine *sM, int64_t state) {
    //Match state is like going to a match.
    state_check(sM, state);
    return state == match1 ? 0 : LOG_ZERO;
}

static double stateMachineEchelon_raggedStartStateProb(StateMachine *sM, int64_t state) {
    state_check(sM, state);
    return state == gapX ? 0 : LOG_ZERO;
}

static double stateMachineEchelon_endStateProb(StateMachine *sM, int64_t state) {
    StateMachineEchelon *sMe = (StateMachineEchelon *) sM;
    state_check(sM, state);
    switch (state) {
        case match0:
        case match1:
        case match2:
        case match3:
        case match4:
        case match5:
            return sMe->DEFAULT_END_MATCH_PROB;
        case gapX:
            return sMe->DEFAULT_END_FROM_X_PROB;
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

void stateMachine3Vanilla_setStrandTransitionsToDefaults(StateMachine *sM, Strand strand) {
    StateMachine3Vanilla *sM3v = (StateMachine3Vanilla *)sM;
    switch (strand) {
        case template:
            sM3v->TRANSITION_M_TO_Y_NOT_X = 0.17f;
            sM3v->TRANSITION_E_TO_E = 0.55f;
            return;
        case complement:
            sM3v->TRANSITION_M_TO_Y_NOT_X = 0.14f;
            sM3v->TRANSITION_E_TO_E = 0.49f;
            return;
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
    StateMachine3 *sM3 = (StateMachine3 *) sM;
    HDCell *hdCurrent = current == NULL ? NULL : (HDCell *)current;
    HDCell *hdLower = lower == NULL ? NULL : (HDCell *)lower;
    HDCell *hdMiddle = middle == NULL ? NULL : (HDCell *)middle;
    HDCell *hdUpper = upper == NULL ? NULL : (HDCell *)upper;
    //st_uglyf("SENTINAL - at begining of cell calc\n");
    if (hdLower != NULL) {
        for (int64_t p = 0; p < hdCurrent->numberOfPaths; p ++) {
            Path *pathC = hdCell_getPath(hdCurrent, p);
            for (int64_t q = 0; q < hdLower->numberOfPaths; q++) {
                Path *pathL = hdCell_getPath(hdLower, q);
                if (path_checkLegal(pathL, pathC)) {
                    //st_uglyf("SENTINAL - legal LOWER : pathC kmer %s\n", pathC->kmer);
                    double *lowerCells = path_getCell(pathL);
                    double *currentCells = path_getCell(pathC);
                    double eP = sM3->getXGapProbFcn(sM3->model.EMISSION_GAP_X_PROBS, pathC->kmer);
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
                    double eP = sM3->getMatchProbFcn(sM3->model.EMISSION_MATCH_PROBS, pathC->kmer, cY);
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
                    double eP = sM3->getYGapProbFcn(sM3->model.EMISSION_GAP_Y_PROBS, pathC->kmer, cY);
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
                    double eP = sM3->getMatchProbFcn(sM3->hdpModel, pathC->kmer, cY);
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
                    double eP = sM3->getYGapProbFcn(sM3->hdpModel, pathC->kmer, cY);
                    doTransition(upperCells, currentCells, match, shortGapY, eP, sM3->TRANSITION_GAP_OPEN_Y, extraArgs);
                    doTransition(upperCells, currentCells, shortGapY, shortGapY, eP, sM3->TRANSITION_GAP_EXTEND_Y, extraArgs);
                    // shortGapX -> shortGapY not allowed, this would be going from a kmer skip to extra event?
                    //doTransition(upper, current, shortGapX, shortGapY, eP, sM3->TRANSITION_GAP_SWITCH_TO_Y, extraArgs);
                }
            }
        }
    }
}

static void stateMachine3Vanilla_cellCalculate(StateMachine *sM,
                                               double *current, double *lower, double *middle, double *upper,
                                               void *cX, void *cY,
                                               void (*doTransition)(double *, double *, int64_t, int64_t,
                                                                    double, double, void *),
                                               void *extraArgs) {

    StateMachine3Vanilla *sM3v = (StateMachine3Vanilla *) sM;
    // Establish transition probs (all adopted from Nanopolish by JTS)
    // from match
    double a_mx = sM3v->getKmerSkipProb((StateMachine *) sM3v, cX, 0); // get beta prob
    double a_my = (1 - a_mx) * sM3v->TRANSITION_M_TO_Y_NOT_X; // trans M to Y not X fudge factor
    double a_mm = 1.0f - a_my - a_mx;

    // from Y [Extra event state]
    double a_yy = sM3v->TRANSITION_E_TO_E;
    double a_ym = 1.0f - a_yy;

    // from X [Skipped event state]
    double a_xx = sM3v->getKmerSkipProb((StateMachine *)sM3v, cX, 1); // get alpha prob
    double a_xm = 1.0f - a_xx;

    if (lower != NULL) {
        doTransition(lower, current, match, shortGapX, 0, log(a_mx), extraArgs);
        doTransition(lower, current, shortGapX, shortGapX, 0, log(a_xx), extraArgs);
        // X to Y not allowed
        //doTransition(lower, current, shortGapY, shortGapX, eP, sM3->TRANSITION_GAP_SWITCH_TO_X, extraArgs);
    }
    if (middle != NULL) {
        double eP = sM3v->getMatchProbFcn(sM3v->model.EMISSION_MATCH_PROBS, cX, cY);
        doTransition(middle, current, match, match, eP, log(a_mm), extraArgs);
        doTransition(middle, current, shortGapX, match, eP, log(a_xm), extraArgs);
        doTransition(middle, current, shortGapY, match, eP, log(a_ym), extraArgs);
    }
    if (upper != NULL) {
        double eP = sM3v->getScaledMatchProbFcn(sM3v->model.EMISSION_GAP_Y_PROBS, cX, cY);
        doTransition(upper, current, match, shortGapY, eP, log(a_my), extraArgs);
        doTransition(upper, current, shortGapY, shortGapY, eP, log(a_yy), extraArgs);
        // Y to X not allowed
        //doTransition(upper, current, shortGapX, shortGapY, eP, sM3->TRANSITION_GAP_SWITCH_TO_Y, extraArgs);
    }
}

static void stateMachineEchelon_cellCalculate(StateMachine *sM,
                                              double *current, double *lower, double *middle, double *upper,
                                              void *cX, void *cY,
                                              void (*doTransition)(double *, double *, int64_t, int64_t,
                                                                   double, double, void *),
                                              void *extraArgs) {
    StateMachineEchelon *sMe = (StateMachineEchelon *) sM;
    // transitions
    // from M
    double a_mx = sMe->getKmerSkipProb((StateMachine *)sMe, cX, 0), la_mx = log(a_mx); // beta
    double a_mh = 1 - a_mx, la_mh = log(a_mh); // 1 - beta

    // from X (kmer skip)
    //double a_xx = a_mx, la_xx = log(a_xx); // alpha, to seperate alpha, need to change here
    double a_xx = sMe->getKmerSkipProb((StateMachine *)sMe, cX, 1), la_xx = log(a_xx);
    double a_xh = 1 - a_xx, la_xh = log(a_xh); // 1 - alpha

    if (lower != NULL) {
        // go from all of the match states to gapX
        for (int64_t n = 1; n < 6; n++) {
            doTransition(lower, current, n, gapX, 0, la_mx, extraArgs);
        }
        // gapX --> gapX
        doTransition(lower, current, gapX, gapX, 0, la_xx, extraArgs);
    }
    if (middle != NULL) {
        // first we handle going from all of the match states to match1 through match5
        for (int64_t n = 1; n < 6; n++) {
            for (int64_t from = 0; from < 6; from++) {
                doTransition(middle, current, from, n,
                             sMe->getMatchProbFcn(sMe->model.EMISSION_MATCH_PROBS, cX, cY, n),
                             (la_mh + sMe->getDurationProb(cY, n)), extraArgs);
            }
        }
        // now do from gapX to the match states
        for (int64_t n = 1; n < 6; n++) {
            doTransition(middle, current, gapX, n, sMe->getMatchProbFcn(sMe->model.EMISSION_MATCH_PROBS, cX, cY, n),
                         (la_xh + sMe->getDurationProb(cY, n)), extraArgs);
        }
    }
    if (upper != NULL) {
        // only allowed to go from match states to match0 (extra event state)
        for (int64_t n = 1; n < 6; n++) {
            doTransition(upper, current, n, match0,
                         sMe->getScaledMatchProbFcn(sMe->model.EMISSION_GAP_Y_PROBS, cX, cY),
                         //sMe->getDurationProb(cY, 0), extraArgs);
                         (la_mh + sMe->getDurationProb(cY, 0)), extraArgs);
        }
    }
}

///////////////////////////////////////////// CORE FUNCTIONS ////////////////////////////////////////////////////////
StateMachine *stateMachine3_construct(StateMachineType type, int64_t parameterSetSize,
                                      void (*setTransitionsToDefaults)(StateMachine *sM),
                                      void (*setEmissionsDefaults)(StateMachine *sM, int64_t nbSkipParams),
                                      double (*gapXProbFcn)(const double *, void *),
                                      double (*gapYProbFcn)(const double *, void *, void *),
                                      double (*matchProbFcn)(const double *, void *, void *),
                                      void (*cellCalcUpdateExpFcn)(double *fromCells, double *toCells,
                                                                   int64_t from, int64_t to,
                                                                   double eP, double tP, void *extraArgs)) {
    /*
     * Description of (potentially ambigious) arguments:
     * parameterSetSize = the number of kmers that we are using, of len(kmer) = 1, then the number is 4 (or 5 if we're
     * including N). It's 25 if len(kmer) = 2, it's 4096 in the 6-mer model.
     */
    StateMachine3 *sM3 = st_malloc(sizeof(StateMachine3));
    if (type != threeState && type != threeStateAsymmetric) {
        st_errAbort("Tried to create a three state state-machine with the wrong type");
    }

    // setup the parent class
    sM3->model.type = type;
    sM3->model.parameterSetSize = parameterSetSize;
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
    setEmissionsDefaults((StateMachine *) sM3, parameterSetSize);

    // set gap probs
    for (int64_t i = 0; i < parameterSetSize; i++) {
        sM3->model.EMISSION_GAP_X_PROBS[i] = -2.3025850929940455; // log(0.1)
    }

    return (StateMachine *) sM3;
}

StateMachine *stateMachine3Hdp_construct(StateMachineType type, int64_t parameterSetSize,
                                         void (*setTransitionsToDefaults)(StateMachine *sM),
                                         void (*setEmissionsDefaults)(StateMachine *sM, int64_t nbSkipParams),
                                         NanoporeHDP *hdpModel,
                                         double (*gapXProbFcn)(const double *, void *),
                                         double (*gapYProbFcn)(NanoporeHDP *, void *, void *),
                                         double (*matchProbFcn)(NanoporeHDP *, void *, void *),
                                         void (*cellCalcUpdateExpFcn)(double *fromCells, double *toCells,
                                                                      int64_t from, int64_t to,
                                                                      double eP, double tP, void *extraArgs)) {
    StateMachine3_HDP *sM3 = st_malloc(sizeof(StateMachine3_HDP));
    if (type != threeStateHdp) {
        st_errAbort("Tried to create a three state state-machine with the wrong type");
    }

    // setup the parent class
    sM3->model.type = type;
    sM3->model.parameterSetSize = parameterSetSize;
    sM3->model.stateNumber = 3;
    sM3->model.matchState = match;
    sM3->model.startStateProb = stateMachine3_startStateProb;
    sM3->model.endStateProb = stateMachine3_endStateProb;
    sM3->model.raggedStartStateProb = stateMachine3_raggedStartStateProb;
    sM3->model.raggedEndStateProb = stateMachine3_raggedEndStateProb;
    sM3->model.cellCalculate = stateMachine3HDP_cellCalculate;
    sM3->model.cellCalculateUpdateExpectations = cellCalcUpdateExpFcn;

    // setup functions
    sM3->getXGapProbFcn = gapXProbFcn;
    sM3->getYGapProbFcn = gapYProbFcn;
    sM3->getMatchProbFcn = matchProbFcn;

    // setup HDP
    sM3->hdpModel = hdpModel;

    // setup transitions
    setTransitionsToDefaults((StateMachine *) sM3);
    // set emissions to defaults or zeros
    setEmissionsDefaults((StateMachine *) sM3, parameterSetSize);

    // initialize kmer emission gap probs
    for (int64_t i = 0; i < parameterSetSize; i++) {
        sM3->model.EMISSION_GAP_X_PROBS[i] = -2.3025850929940455; // log(0.1)
    }
    return (StateMachine *) sM3;
}

StateMachine *stateMachine3Vanilla_construct(StateMachineType type, int64_t parameterSetSize,
                                             void (*setEmissionsDefaults)(StateMachine *sM, int64_t nbSkipParams),
                                             double (*xSkipProbFcn)(StateMachine *, void *, bool),
                                             double (*scaledMatchProbFcn)(const double *, void *, void *),
                                             double (*matchProbFcn)(const double *, void *, void *),
                                             void (*cellCalcUpdateExpFcn)(double *fromCells, double *toCells,
                                                                          int64_t from, int64_t to,
                                                                          double eP, double tP, void *extraArgs)) {
    StateMachine3Vanilla *sM3v = st_malloc(sizeof(StateMachine3Vanilla));
    // check
    if (type != vanilla) {
        st_errAbort("Tried to create a vanilla state machine with the wrong type?");
    }
    // setup end state prob defaults
    // defaults from a template nanopore file
    sM3v->TRANSITION_M_TO_Y_NOT_X = 0.17; //default for template
    sM3v->TRANSITION_E_TO_E = 0.55f; //default for template
    sM3v->DEFAULT_END_MATCH_PROB = -0.23552123624314988; // log(step_prob) (0.79015888282447311)
    sM3v->DEFAULT_END_FROM_X_PROB = -1.6269694202638481; // log(skip_prob) (0.19652425498269727)
    sM3v->DEFAULT_END_FROM_Y_PROB = -4.3187242127300092; // log(1 - (skip_prob + step_prob)) (0.013316862192829682)
    // setup the parent class
    sM3v->model.type = type;
    sM3v->model.parameterSetSize = parameterSetSize;
    sM3v->model.stateNumber = 3;
    sM3v->model.matchState = match;
    sM3v->model.startStateProb = stateMachine3_startStateProb;
    sM3v->model.raggedStartStateProb = stateMachine3_raggedStartStateProb;
    sM3v->model.endStateProb = stateMachine3Vanilla_endStateProb;
    sM3v->model.raggedEndStateProb = stateMachine3Vanilla_raggedEndStateProb;
    sM3v->model.cellCalculate = stateMachine3Vanilla_cellCalculate;
    sM3v->model.cellCalculateUpdateExpectations = cellCalcUpdateExpFcn;

    // stateMachine3Vanilla-specific functions
    sM3v->getKmerSkipProb = xSkipProbFcn;
    sM3v->getScaledMatchProbFcn = scaledMatchProbFcn;
    sM3v->getMatchProbFcn = matchProbFcn;

    // set emissions to defaults or zeros
    setEmissionsDefaults((StateMachine *) sM3v, 60);
    return (StateMachine *) sM3v;
}

StateMachine *stateMachineEchelon_construct(StateMachineType type, int64_t parameterSetSize,
                                            void (*setEmissionsToDefaults)(StateMachine *sM, int64_t nbSkipParams),
                                            double (*durationProbFcn)(void *event, int64_t n),
                                            double (*skipProbFcn)(StateMachine *sM, void *kmerList, bool),
                                            double (*matchProbFcn)(const double *, void *, void *, int64_t n),
                                            double (*scaledMatchProbFcn)(const double *, void *, void *),
                                            void (*cellCalcUpdateExpFcn)(double *fromCells, double *toCells,
                                                                         int64_t from, int64_t to,
                                                                         double eP, double tP, void *extraArgs)) {
    StateMachineEchelon *sMe = st_malloc(sizeof(StateMachineEchelon));

    if (type != echelon) {
        st_errAbort("Tried to create a echelon state machine with the wrong type?");
    }

    // setup end state probabilities // todo these aren't log and won't work
    sMe->DEFAULT_END_MATCH_PROB = 0.79015888282447311; // stride_prb
    sMe->DEFAULT_END_FROM_X_PROB = 0.19652425498269727; // skip_prob
    sMe->BACKGROUND_EVENT_PROB = -3.0;

    // parent class setup
    sMe->model.type = type;
    sMe->model.parameterSetSize = parameterSetSize;
    sMe->model.stateNumber = 7;
    sMe->model.matchState = match1; // take a look at this, might need to change
    sMe->model.startStateProb = stateMachineEchelon_startStateProb;
    sMe->model.raggedStartStateProb = stateMachineEchelon_raggedStartStateProb;
    sMe->model.endStateProb = stateMachineEchelon_endStateProb;
    sMe->model.raggedEndStateProb = stateMachineEchelon_endStateProb;
    sMe->model.cellCalculate = stateMachineEchelon_cellCalculate;
    sMe->model.cellCalculateUpdateExpectations = cellCalcUpdateExpFcn;

    // class functions
    sMe->getKmerSkipProb = skipProbFcn;
    sMe->getDurationProb = durationProbFcn;
    sMe->getMatchProbFcn = matchProbFcn;
    sMe->getScaledMatchProbFcn = scaledMatchProbFcn;

    setEmissionsToDefaults((StateMachine *) sMe, 60);
    return (StateMachine *) sMe;
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
    emissions_em_loadMatchProbs(sM3->EMISSION_MATCH_PROBS, hmm, match);
    int64_t xGapStates[1] = { shortGapX };
    int64_t yGapStates[1] = { shortGapY };
    emissions_loadGapProbs(sM3->EMISSION_GAP_X_PROBS, hmm, xGapStates, 1, NULL, 0);
    emissions_loadGapProbs(sM3->EMISSION_GAP_Y_PROBS, hmm, NULL, 0, yGapStates, 1);
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
    emissions_em_loadMatchProbsSymmetrically(sM3->EMISSION_MATCH_PROBS, hmm, match);
    int64_t xGapStates[2] = { shortGapX };
    int64_t yGapStates[2] = { shortGapY };
    emissions_loadGapProbs(sM3->EMISSION_GAP_X_PROBS, hmm, xGapStates, 1, yGapStates, 1);
    emissions_em_loadGapProbs(sM3->EMISSION_GAP_Y_PROBS, hmm, xGapStates, 1, yGapStates, 1);
}
*/

///////////////////////////////////
///////////////////////////////////
//Public functions
///////////////////////////////////
///////////////////////////////////

StateMachine *getStateMachine5(Hmm *hmmD, StateMachineFunctions *sMfs) {
    if (hmmD->type == fiveState) {
        StateMachine5 *sM5 = (StateMachine5 *) stateMachine5_construct(fiveState, hmmD->symbolSetSize,
                                                                       emissions_discrete_initEmissionsToZero,
                                                                       sMfs->gapXProbFcn,
                                                                       sMfs->gapYProbFcn,
                                                                       sMfs->matchProbFcn,
                                                                       cell_updateExpectations);
        stateMachine5_loadSymmetric(sM5, hmmD);
        return (StateMachine *) sM5;
    }
    if (hmmD->type == fiveStateAsymmetric) {
        StateMachine5 *sM5 = (StateMachine5 *) stateMachine5_construct(fiveState, hmmD->symbolSetSize,
                                                                       emissions_discrete_initEmissionsToZero,
                                                                       sMfs->gapXProbFcn,
                                                                       sMfs->gapYProbFcn,
                                                                       sMfs->matchProbFcn,
                                                                       cell_updateExpectations);
        stateMachine5_loadAsymmetric(sM5, hmmD);
        return (StateMachine *) sM5;
    }
    else {
        st_uglyf("SENTINAL - getStateMachine5 --> something didn't work\n");
        return 0;
    }
}

StateMachine *getStrawManStateMachine3(const char *modelFile) {
    if (!stFile_exists(modelFile)) {
        st_errAbort("getStrawManStateMachine3: Cannot find model file %s\n", modelFile);
    };
    StateMachine *sM3 = stateMachine3_construct(threeState, NUM_OF_KMERS,
                                                stateMachine3_setTransitionsToNanoporeDefaults,
                                                emissions_signal_initEmissionsToZero,
                                                emissions_kmer_getGapProb,
                                                emissions_signal_strawManGetKmerEventMatchProb,
                                                emissions_signal_strawManGetKmerEventMatchProb,
                                                cell_signal_updateTransAndKmerSkipExpectations);
    emissions_signal_loadPoreModel(sM3, modelFile, sM3->type);
    return sM3;
}


StateMachine *getHdpStateMachine3(NanoporeHDP *hdp, const char *modelFile) {
    StateMachine *sM3Hdp = stateMachine3Hdp_construct(threeStateHdp, NUM_OF_KMERS,
                                                      stateMachine3_setTransitionsToNanoporeDefaults,
                                                      emissions_signal_initEmissionsToZero,
                                                      hdp,
                                                      emissions_kmer_getGapProb,
                                                      get_nanopore_kmer_density,
                                                      get_nanopore_kmer_density,
                                                      cell_signal_updateTransAndKmerSkipExpectations2);
    emissions_signal_loadPoreModel(sM3Hdp, modelFile, sM3Hdp->type);
    return sM3Hdp;
}

StateMachine *getStateMachine4(const char *modelFile) {
    StateMachine *sM4 = stateMachine4_construct(fourState, NUM_OF_KMERS,
                                                emissions_signal_initEmissionsToZero,
                                                emissions_kmer_getGapProb,
                                                emissions_signal_strawManGetKmerEventMatchProb,
                                                emissions_signal_strawManGetKmerEventMatchProb,
                                                cell_signal_updateTransAndKmerSkipExpectations);
    emissions_signal_loadPoreModel(sM4, modelFile, sM4->type);
    return sM4;
}

StateMachine *getSignalStateMachine3Vanilla(const char *modelFile) {
    // construct a stateMachine3Vanilla then load the model
    StateMachine *sM3v = stateMachine3Vanilla_construct(vanilla, NUM_OF_KMERS,
                                                        emissions_signal_initEmissionsToZero,
                                                        emissions_signal_getBetaOrAlphaSkipProb,
                                                        emissions_signal_getEventMatchProbWithTwoDists,
                                                        emissions_signal_getEventMatchProbWithTwoDists,
                                                        cell_signal_updateBetaAndAlphaProb);
    emissions_signal_loadPoreModel(sM3v, modelFile, sM3v->type);
    return sM3v;
}

StateMachine *getStateMachineEchelon(const char *modelFile) {
    StateMachine *sMe = stateMachineEchelon_construct(echelon, NUM_OF_KMERS,
                                                      emissions_signal_initEmissionsToZero,
                                                      emissions_signal_getDurationProb,
                                                      //emissions_signal_getKmerSkipProb,
                                                      emissions_signal_getBetaOrAlphaSkipProb,
                                                      emissions_signal_multipleKmerMatchProb,
                                                      emissions_signal_getEventMatchProbWithTwoDists,
                                                      NULL); // cell update expectation, to be implemented
    emissions_signal_loadPoreModel(sMe, modelFile, sMe->type);
    return sMe;
}

void stateMachine_destruct(StateMachine *stateMachine) {
    free(stateMachine);
}


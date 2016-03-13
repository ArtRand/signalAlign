/*
 * stateMachine.h
 *
 *  Created on: Aug 8, 2014
 *      Author: benedictpaten
 */

#ifndef STATEMACHINE_H_
#define STATEMACHINE_H_

#include "sonLib.h"
#include "nanopore_hdp.h"

#define SYMBOL_NUMBER 5
#define SYMBOL_NUMBER_NO_N 6
#define SYMBOL_NUMBER_EPIGENETIC_C 6
#define MODEL_PARAMS 5 // level_mean, level_sd, fluctuation_mean, fluctuation_noise, fluctuation_lambda
#define SIX_LETTER_ALPHA "ACEGOT"


typedef enum {
    fiveState = 0,
    fiveStateAsymmetric = 1,
    threeState = 2,
    threeStateAsymmetric = 3,
    vanilla = 4,
    echelon = 5,
    fourState = 6,
    threeStateHdp = 7,
} StateMachineType;

typedef enum {
    match = 0, shortGapX = 1, shortGapY = 2, longGapX = 3, longGapY = 4
} State;

typedef enum _strand {
    template = 0,
    complement = 1,
} Strand;

typedef struct _stateMachine StateMachine;
typedef struct _hmm Hmm;

/*
 * Hmm for loading/unloading HMMs and storing expectations.
 * Maybe move these definitions to stateMachine.c to clean this up?
 */
struct _hmm {
    double likelihood;
    StateMachineType type;
    int64_t stateNumber;
    int64_t symbolSetSize;
    int64_t matrixSize;

    // 10/27 moved into HmmDiscrete subClass
    //double *transitions;
    //double *emissions;

    void (*addToTransitionExpectationFcn)(Hmm *hmm, int64_t from, int64_t to, double p);

    void (*setTransitionFcn)(Hmm *hmm, int64_t from, int64_t to, double p);

    double (*getTransitionsExpFcn)(Hmm *hmm, int64_t from, int64_t to);

};

struct _stateMachine {
    StateMachineType type;
    int64_t stateNumber;
    int64_t matchState;
    int64_t parameterSetSize;

    double *EMISSION_MATCH_MATRIX; //Match emission probs
    double *EMISSION_GAP_X_PROBS;  //Gap emission probs
    double *EMISSION_GAP_Y_MATRIX; //Gap emission probs

    double (*startStateProb)(StateMachine *sM, int64_t state);

    double (*endStateProb)(StateMachine *sM, int64_t state);

    double (*raggedEndStateProb)(StateMachine *sM, int64_t state);

    double (*raggedStartStateProb)(StateMachine *sM, int64_t state);

    //Cells (states at a given coordinate)
    void (*cellCalculate)(StateMachine *sM, void *current, void *lower, void *middle, void *upper,
                          void* cX, void* cY,
                          void(*doTransition)(double *, double *, int64_t, int64_t, double, double, void *),
                          void *extraArgs);

    void (*cellCalculateUpdateExpectations) (double *fromCells, double *toCells, int64_t from, int64_t to,
                                             double eP, double tP, void *extraArgs);
};

typedef struct _StateMachine5 StateMachine5;

struct _StateMachine5 {
    StateMachine model;
    double TRANSITION_MATCH_CONTINUE; //0.9703833696510062f
    double TRANSITION_MATCH_FROM_SHORT_GAP_X; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_MATCH_FROM_LONG_GAP_X; //1.0 - gapExtend = 0.00343657420938
    double TRANSITION_GAP_SHORT_OPEN_X; //0.0129868352330243
    double TRANSITION_GAP_SHORT_EXTEND_X; //0.7126062401851738f;
    double TRANSITION_GAP_SHORT_SWITCH_TO_X; //0.0073673675173412815f;
    double TRANSITION_GAP_LONG_OPEN_X; //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    double TRANSITION_GAP_LONG_EXTEND_X; //0.99656342579062f;
    double TRANSITION_GAP_LONG_SWITCH_TO_X; //0.0073673675173412815f;
    double TRANSITION_MATCH_FROM_SHORT_GAP_Y; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_MATCH_FROM_LONG_GAP_Y; //1.0 - gapExtend = 0.00343657420938
    double TRANSITION_GAP_SHORT_OPEN_Y; //0.0129868352330243
    double TRANSITION_GAP_SHORT_EXTEND_Y; //0.7126062401851738f;
    double TRANSITION_GAP_SHORT_SWITCH_TO_Y; //0.0073673675173412815f;
    double TRANSITION_GAP_LONG_OPEN_Y; //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    double TRANSITION_GAP_LONG_EXTEND_Y; //0.99656342579062f;
    double TRANSITION_GAP_LONG_SWITCH_TO_Y; //0.0073673675173412815f;

    double (*getXGapProbFcn)(const double *emissionXGapProbs, void *i);
    double (*getYGapProbFcn)(const double *emissionYGapProbs, void *i);
    double (*getMatchProbFcn)(const double *emissionMatchProbs, void *x, void *y);

};

typedef struct _StateMachine3 StateMachine3;

struct _StateMachine3 {
    // 3 state state machine, allowing for symmetry in x and y.
    StateMachine model;
    double TRANSITION_MATCH_CONTINUE; //0.9703833696510062f
    double TRANSITION_MATCH_FROM_GAP_X; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_MATCH_FROM_GAP_Y; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_GAP_OPEN_X; //0.0129868352330243
    double TRANSITION_GAP_OPEN_Y; //0.0129868352330243
    double TRANSITION_GAP_EXTEND_X; //0.7126062401851738f;
    double TRANSITION_GAP_EXTEND_Y; //0.7126062401851738f;
    double TRANSITION_GAP_SWITCH_TO_X; //0.0073673675173412815f;
    double TRANSITION_GAP_SWITCH_TO_Y; //0.0073673675173412815f;

    double (*getXGapProbFcn)(const double *emissionXGapProbs, void *i);

    double (*getYGapProbFcn)(const double *emissionYGapProbs, void *x, void *y);
    double (*getMatchProbFcn)(const double *emissionMatchProbs, void *x, void *y);
};

typedef struct _StateMachine3_HDP StateMachine3_HDP;

struct _StateMachine3_HDP {
    StateMachine model;
    double TRANSITION_MATCH_CONTINUE; //0.9703833696510062f
    double TRANSITION_MATCH_FROM_GAP_X; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_MATCH_FROM_GAP_Y; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_GAP_OPEN_X; //0.0129868352330243
    double TRANSITION_GAP_OPEN_Y; //0.0129868352330243
    double TRANSITION_GAP_EXTEND_X; //0.7126062401851738f;
    double TRANSITION_GAP_EXTEND_Y; //0.7126062401851738f;
    double TRANSITION_GAP_SWITCH_TO_X; //0.0073673675173412815f;
    double TRANSITION_GAP_SWITCH_TO_Y; //0.0073673675173412815f;

    double (*getXGapProbFcn)(const double *emissionXGapProbs, void *i);
    NanoporeHDP *hdpModel;
    double (*getYGapProbFcn)(NanoporeHDP *hdp, void *x, void *y);
    double (*getMatchProbFcn)(NanoporeHDP *hdp, void *x, void *y);
};

typedef struct _StateMachineEchelon {
    // 8-state general hmm
    StateMachine model;

    double BACKGROUND_EVENT_PROB;
    double DEFAULT_END_MATCH_PROB; //0.79015888282447311; // stride_prb
    double DEFAULT_END_FROM_X_PROB; //0.19652425498269727; // skip_prob

    double (*getKmerSkipProb)(StateMachine *sM, void *kmerList, bool); // beta
    double (*getDurationProb)(void *event, int64_t n); // P(dj|n)
    double (*getMatchProbFcn)(const double *eventModel, void *kmers, void *event, int64_t n); // P(ej|xi..xn)
    double (*getScaledMatchProbFcn)(const double *scaledEventModel, void *kmer, void *event);

} StateMachineEchelon;

typedef struct _stateMachineFunctions { // TODO rename these to be more general or might remove them all together
    double (*gapXProbFcn)(const double *, void *);
    double (*gapYProbFcn)(const double *, void *);
    double (*matchProbFcn)(const double *, void *, void *);
} StateMachineFunctions;

//////////////////
// stateMachine //
//////////////////

// StateMachine constructors
StateMachine *stateMachine5_construct(StateMachineType type, int64_t parameterSetSize,
                                      void (*setEmissionsDefaults)(StateMachine *sM),
                                      double (*gapXProbFcn)(const double *, void *),
                                      double (*gapYProbFcn)(const double *, void *),
                                      double (*matchProbFcn)(const double *, void *, void *),
                                      void (*cellCalcUpdateExpFcn)(double *fromCells, double *toCells,
                                                                   int64_t from, int64_t to,
                                                                   double eP, double tP, void *extraArgs));

StateMachine *stateMachine3Hdp_construct(StateMachineType type, int64_t parameterSetSize,
                                         void (*setTransitionsToDefaults)(StateMachine *sM),
                                         void (*setEmissionsDefaults)(StateMachine *sM, int64_t nbSkipParams),
                                         NanoporeHDP *hdpModel,
                                         double (*gapXProbFcn)(const double *, void *),
                                         double (*gapYProbFcn)(NanoporeHDP *, void *, void *),
                                         double (*matchProbFcn)(NanoporeHDP *, void *, void *),
                                         void (*cellCalcUpdateExpFcn)(double *fromCells, double *toCells,
                                                                      int64_t from, int64_t to,
                                                                      double eP, double tP, void *extraArgs));

StateMachine *stateMachine3_construct(StateMachineType type, int64_t parameterSetSize,
                                      void (*setTransitionsToDefaults)(StateMachine *sM),
                                      void (*setEmissionsDefaults)(StateMachine *sM, int64_t nbSkipParams),
                                      double (*gapXProbFcn)(const double *, void *),
                                      double (*gapYProbFcn)(const double *, void *, void *),
                                      double (*matchProbFcn)(const double *, void *, void *),
                                      void (*cellCalcUpdateExpFcn)(double *fromCells, double *toCells,
                                                                   int64_t from, int64_t to,
                                                                   double eP, double tP, void *extraArgs));

StateMachine *stateMachineEchelon_construct(StateMachineType type, int64_t parameterSetSize,
                                            void (*setEmissionsToDefaults)(StateMachine *sM, int64_t nbSkipParams),
                                            double (*durationProbFcn)(void *event, int64_t n),
                                            double (*skipProbFcn)(StateMachine *sM, void *kmerList, bool),
                                            double (*matchProbFcn)(const double *, void *, void *, int64_t n),
                                            double (*scaledMatchProbFcn)(const double *, void *, void *),
                                            void (*cellCalcUpdateExpFcn)(double *fromCells, double *toCells,
                                                                         int64_t from, int64_t to,
                                                                         double eP, double tP, void *extraArgs));

// indexing functions //
void stateMachine_index_check(int64_t c);
//Returns the index for a base, for use with matrices and emissions_discrete_getKmerIndex
int64_t emissions_discrete_getBaseIndex(void *base);

//Returns the index for a kmer from pointer to kmer string
int64_t emissions_discrete_getKmerIndex(void *kmer);

// Returns index of kmer from pointer to array
int64_t emissions_discrete_getKmerIndexFromKmer(void *kmer);

// transition defaults
void stateMachine3_setTransitionsToNucleotideDefaults(StateMachine *sM);

void stateMachine3_setTransitionsToNanoporeDefaults(StateMachine *sM);

// emissions defaults
void emissions_discrete_initEmissionsToZero(StateMachine *sM);

//void emissions_symbol_setEmissionsToDefaults(StateMachine *sM);

/*
* For a discrete HMM the gap and match matrices are defined by the number of symbols in the set (nK). The gap
* matrix is nK x 1 and the match matrix is nK x nK
* In the most simple case, with 4 nucleotides the gap matrix is 4x1 matrix and the match matrix is a 4x4 matrix.
*/
void emissions_signal_initEmissionsToZero(StateMachine *sM, int64_t nbSkipParams);

double emissions_symbol_getGapProb(const double *emissionGapProbs, void *base);

double emissions_symbol_getMatchProb(const double *emissionMatchProbs, void *x, void *y);

double emissions_kmer_getGapProb(const double *emissionGapProbs, void *x_i);

double emissions_kmer_getMatchProb(const double *emissionMatchProbs, void *x, void *y);

int64_t emissions_signal_getKmerSkipBin(double *matchModel, void *kmers);

double emissions_signal_getBetaOrAlphaSkipProb(StateMachine *sM, void *kmers, bool getAlpha);

double emissions_signal_getKmerSkipProb(StateMachine *sM, void *kmers);

double emissions_signal_logGaussMatchProb(const double *eventModel, void *kmer, void *event);

// returns log of the probability density function for a Gaussian distribution
double emissions_signal_getBivariateGaussPdfMatchProb(const double *eventModel, void *kmer, void *event);

double emissions_signal_getEventMatchProbWithTwoDists(const double *eventModel, void *kmer, void *event);

void emissions_signal_scaleModel(StateMachine *sM, double scale, double shift, double var,
                                 double scale_sd, double var_sd);

void emissions_signal_scaleEmissions(StateMachine *sM, double scale, double shift, double var);

double emissions_signal_getDurationProb(void *event, int64_t n);

StateMachine *getStrawManStateMachine3(const char *modelFile);

StateMachine *getHdpStateMachine3(NanoporeHDP *hdp, const char *modelFile);

StateMachine *getStateMachineEchelon(const char *modelFile);

// EM
StateMachine *getStateMachine5(Hmm *hmmD, StateMachineFunctions *sMfs);

void stateMachine_destruct(StateMachine *stateMachine);

#endif /* STATEMACHINE_H_ */

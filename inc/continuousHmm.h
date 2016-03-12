#ifndef CONTINUOUS_HMM_H
#define CONTINUOUS_HMM_H
#define NORMAL_DISTRIBUTION_PARAMS 2
#include "stateMachine.h"

typedef struct _strawManHmm {
    Hmm baseHmm;
    // threeState transitions matrix (3x3)
    double *transitions;
    // matrix formatted [x, y] * NUM_OF_KMERS. x = Sigma_k(p_i * obs_i) and y = Sigma_k(p_i * (obs - u_K)**2),
    // u_k is calculated as x / Sigma_k(p_k) and _k indicates the kmer
    double *eventExpectations;
    // matrix formatted [E(u), E(o)] * NUM_OF_KMERS, this is the expected parameters for the Normal Distributions
    // for each kmer
    double *eventModel;
    // matrix formatted [Sigma(p_k)] * NUM_OF_KMERS, a running sum of the posterior probabilities
    double *posteriors;
    // this read's scale parameter and shift parameter, these are used to descale the events when adding to
    // the expectation
    double scale;
    double shift;
    // gets the kmer index for a nucleotide string
    int64_t (*getElementIndexFcn)(void *);
    // descales the event mean
    double (*getDescaledEvent)(double scale, double shift, double eventMean);
    // updates the eventExpectations matrix
    void (*addToEmissionExpectationFcn)(Hmm *hmm, int64_t kmerIndex, double meanCurrent, double p);
    // sets a value in the eventExpectations matrix
    void (*setEmissionExpectationFcn)(Hmm *hmm, int64_t kmerIndex, double meanCurrent, double p);
    // these functions return pointers to their respective matrices
    double *(*getEmissionExpFcn)(Hmm *hmm, int64_t kmerIndex);
    double *(*getPosteriorExpFcn)(Hmm *hmm, int64_t kmerIndex);
    double *(*getEventModelEntry)(Hmm *hmm, int64_t kmerIndex);
    // a boolian mask for the kmers that have been observed
    bool *observed;
    // boolian indicating whether the HMM object has a an event model or not
    bool hasModel;
} ContinuousPairHmm;

typedef struct _hdpHmm {
    Hmm baseHmm;
    double *transitions;
    double threshold;
    void (*addToAssignments)(Hmm *, void *, void *);
    stList *eventAssignments;
    stList *kmerAssignments;
    int64_t numberOfAssignments;
    NanoporeHDP *nhdp;
} HdpHmm;

Hmm *continuousPairHmm_construct(double transitionPseudocount, double emissionPseudocount,
                                 int64_t stateNumber, int64_t symbolSetSize, StateMachineType type,
                                 double scale, double shift);

void continuousPairHmm_addToTransitionsExpectation(Hmm *hmm, int64_t from, int64_t to, double p);

void continuousPairHmm_setTransitionExpectation(Hmm *hmm, int64_t from, int64_t to, double p);

double continuousPairHmm_getTransitionExpectation(Hmm *hmm, int64_t from, int64_t to);

void continuousPairHmm_addToEmissionExpectation(Hmm *hmm, int64_t kmerIndex, double meanCurrent, double p);

void continuousPairHmm_setEmissionExpectation(Hmm *hmm, int64_t kmerIndex, double meanCurrent, double p);

double *continuousPairHmm_getEmissionExpectation(Hmm *hmm, int64_t kmerIndex);

double *continuousPairHmm_getEmissionPosterior(Hmm *hmm, int64_t kmerIndex);

void continuousPairHmm_loadEmissionsIntoStateMachine(StateMachine *sM, Hmm *hmm);

void continuousPairHmm_loadTransitionsIntoStateMachine(StateMachine *sM, Hmm *hmm);

void continuousPairHmm_normalize(Hmm *hmm);

void continuousPairHmm_randomize(Hmm *hmm);

void continuousPairHmm_destruct(Hmm *hmm);

void continuousPairHmm_dump(Hmm *hmm, FILE *fileHandle);

Hmm *continuousPairHmm_loadFromFile(const char *fileName, double transitionPseudocount, double emissionPseudocount);

void continuousPairHmm_loadModelFromFile(ContinuousPairHmm *hmm, const char *modelFile);

Hmm *hdpHmm_constructEmpty(double pseudocount, int64_t stateNumber, StateMachineType type, double threshold,
                           void (*addToTransitionExpFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
                           void (*setTransitionFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
                           double (*getTransitionsExpFcn)(Hmm *hmm, int64_t from, int64_t to));

void hdpHmm_loadTransitions(StateMachine *sM, Hmm *hmm);

void hdpHmm_writeToFile(Hmm *hmm, FILE *fileHandle);

Hmm *hdpHmm_loadFromFile(const char *fileName, NanoporeHDP *ndhp);

void hdpHmm_destruct(Hmm *hmm);

// CORE
void hmmContinuous_loadSignalHmm(const char *hmmFile, StateMachine *sM, StateMachineType type);

void hmmContinuous_destruct(Hmm *hmm, StateMachineType type);

Hmm *hmmContinuous_getEmptyHmm(StateMachineType type, double pseudocount, double threshold,
                               double scale, double shift);

void hmmContinuous_normalize(Hmm *hmm, StateMachineType type);

void hmmContinuous_writeToFile(const char *outFile, Hmm *hmm, StateMachineType type);

int64_t hmmContinuous_howManyAssignments(Hmm *hmm);

#endif

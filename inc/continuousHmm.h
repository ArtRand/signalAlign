#ifndef CONTINUOUS_HMM_H
#define CONTINUOUS_HMM_H
#define NORMAL_DISTRIBUTION_PARAMS 2
#include "stateMachine.h"

typedef struct _strawManHmm ContinuousPairHmm;

struct _strawManHmm {
    Hmm baseHmm;
    // threeState transitions matrix (3x3)
    // matrix formatted [x, y] * NUM_OF_KMERS. x = Sigma_k(p_i * obs_i) and y = Sigma_k(p_i * (obs - u_K)**2),
    // u_k is calculated as x / Sigma_k(p_k) and _k indicates the kmer
    double *eventExpectations;
    // matrix formatted [E(u), E(o)] for each kmer, this is the expected parameters for the Normal Distributions
    // for each kmer
    double *eventModel;
    // matrix formatted [Sigma(p_k)] * NUM_OF_KMERS, a running sum of the posterior probabilities
    double *posteriors;
    // this read's scale parameter and shift parameter, these are used to descale the events when adding to
    // the expectation
    double scale;
    double shift;
    double var;

    // gets the kmer index for a nucleotide string
    int64_t (*getElementIndexFcn)(char *kmer, char *alphabet, int64_t alphabet_size, int64_t kmer_length);
    // descales the event mean
    double (*getDescaledEvent)(double scaledEvent, double levelMean, double scale, double shift, double var);

    // updates the eventExpectations matrix
    void (*addToEmissionExpectationFcn)(Hmm *hmm, int64_t kmerIndex, double meanCurrent, double p);

    // sets a value in the eventExpectations matrix
    void (*setEmissionExpectationFcn)(Hmm *hmm, int64_t kmerIndex, double meanCurrent, double p);

    // these functions return pointers to their respective matrices
    double *(*getEmissionExpFcn)(Hmm *hmm, int64_t kmerIndex);
    double *(*getPosteriorExpFcn)(Hmm *hmm, int64_t kmerIndex);
    double *(*getEventModelEntry)(Hmm *hmm, int64_t kmerIndex);

    // a boolean mask for the kmers that have been observed
    bool *observed;
    // indicates whether the HMM object has a an event model or not
    bool hasModel;
};

typedef struct _hdpHmm HdpHmm;
struct _hdpHmm {
    Hmm baseHmm;

    // posterior probability threshold required to make assignment
    double threshold;

    // re-scaling parameters
    double scale;
    double shift;
    double var;

    // matrix formatted [E(u), E(o)] for each kmer, this is the expected parameters for the Normal Distributions
    // for each kmer
    double *eventModel;

    // gets the kmer index for a nucleotide string (for looking up the E(level_mean)
    int64_t (*getElementIndexFcn)(char* kmer, char* alphabet, int64_t alphabet_size, int64_t kmer_length);

    // descales the event mean
    double (*getDescaledEvent)(double scaledEvent, double levelMean, double scale, double shift, double var);

    // add to assignments
    void (*addToAssignments)(HdpHmm *, void *, void *);
    stList *eventAssignments; // list of *double
    stList *kmerAssignments; // list of *char (reference sequence usually)
    int64_t numberOfAssignments; // total
    NanoporeHDP *nhdp;
};

Hmm *continuousPairHmm_construct(double transitionPseudocount, double emissionPseudocount,
                                 int64_t stateNumber, int64_t alphabetSize, char *alphabet, int64_t kmerLength,
                                 StateMachineType type, double scale, double shift, double var);

void hmm_addToTransitionsExpectation(Hmm *hmm, int64_t from, int64_t to, double p);

void hmm_setTransitionExpectation(Hmm *hmm, int64_t from, int64_t to, double p);

double hmm_getTransitionExpectation(Hmm *hmm, int64_t from, int64_t to);

void continuousPairHmm_addToEmissionExpectation(Hmm *hmm, int64_t kmerIndex, double meanCurrent, double p);

void continuousPairHmm_setEmissionExpectation(Hmm *hmm, int64_t kmerIndex, double meanCurrent, double p);

double *continuousPairHmm_getEmissionExpectation(Hmm *hmm, int64_t kmerIndex);

double *continuousPairHmm_getEmissionPosterior(Hmm *hmm, int64_t kmerIndex);

void continuousPairHmm_loadEmissionsIntoStateMachine(StateMachine *sM, Hmm *hmm);

void continuousPairHmm_loadTransitionsIntoStateMachine(StateMachine *sM, Hmm *hmm);

void continuousPairHmm_normalize(Hmm *hmm);

void continuousPairHmm_randomize(Hmm *hmm);

void continuousPairHmm_destruct(Hmm *hmm);

void continuousPairHmm_writeToFile(Hmm *hmm, FILE *fileHandle);

Hmm *continuousPairHmm_loadFromFile(const char *fileName, StateMachineType type,
                                    double transitionPseudocount, double emissionPseudocount);

Hmm *continuousPairHmm_makeExpectationsHmm(StateMachine *sM);

void continuousPairHmm_loadModelFromFile(ContinuousPairHmm *hmm, const char *modelFile);

Hmm *hdpHmm_constructEmpty(double transitionPseudocount,
                           int64_t stateNumber, int64_t alphabetSize, char *alphabet, int64_t kmerLength,
                           StateMachineType type,
                           double threshold, double scale, double shift, double var);

Hmm *hdpHmm_makeExpectationsHmm(StateMachine *sM, double threshold);

void hdpHmm_writeToFile(Hmm *hmm, FILE *fileHandle);

Hmm *hdpHmm_loadFromFile(const char *fileName, StateMachineType type, NanoporeHDP *nHdp);

void hdpHmm_destruct(Hmm *hmm);

// CORE
void hmmContinuous_loadSignalHmm(const char *hmmFile, StateMachine *sM, StateMachineType type, Hmm *expectations);

void hmmContinuous_destruct(Hmm *hmm, StateMachineType type);

Hmm *hmmContinuous_getExpectationsHmm(StateMachine *sM, double threshold);

void hmmContinuous_normalize(Hmm *hmm, StateMachineType type);

void hmmContinuous_loadModelPrior(Hmm *prior, Hmm *expectations);

void hmmContinuous_writeToFile(const char *outFile, Hmm *hmm, StateMachineType type);

int64_t hmmContinuous_howManyAssignments(Hmm *hmm);

#endif

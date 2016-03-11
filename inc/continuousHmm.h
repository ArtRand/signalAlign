#ifndef CONTINUOUS_HMM_H
#define CONTINUOUS_HMM_H
#define NORMAL_DISTRIBUTION_PARAMS 2
#include "stateMachine.h"


typedef struct _hmmContinuous {
    Hmm baseHmm;
    //double *matchModel;
    //double *extraEventMatchModel;
} HmmContinuous;

typedef struct _strawManHmm {
    Hmm baseHmm; // TODO TODO COMMENT HOW THIS WORKS
    double *transitions;
    double *eventExpectations;
    double *eventModel;
    double *posteriors;
    int64_t (*getElementIndexFcn)(void *);
    void (*addToEmissionExpectationFcn)(Hmm *hmm, int64_t kmerIndex, double meanCurrent, double p);
    void (*setEmissionExpectationFcn)(Hmm *hmm, int64_t kmerIndex, double meanCurrent, double p);
    double *(*getEmissionExpFcn)(Hmm *hmm, int64_t kmerIndex);
    double *(*getPosteriorExpFcn)(Hmm *hmm, int64_t kmerIndex);
    double *(*getEventModelEntry)(Hmm *hmm, int64_t kmerIndex);
    bool hasModel;
    bool hasExpectations;
} ContinuousPairHmm;

typedef struct _vanillaHmm {
    HmmContinuous baseContinuousHmm;
    double *matchModel;
    double *scaledMatchModel;
    double *kmerSkipBins;
    int64_t (*getKmerSkipBin)(double *matchModel, void *cX);
} VanillaHmm;

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
                                 int64_t stateNumber, int64_t symbolSetSize, StateMachineType type);

void continuousPairHmm_addToTransitionsExpectation(Hmm *hmm, int64_t from, int64_t to, double p);

void continuousPairHmm_setTransitionExpectation(Hmm *hmm, int64_t from, int64_t to, double p);

double continuousPairHmm_getTransitionExpectation(Hmm *hmm, int64_t from, int64_t to);

void continuousPairHmm_addToEmissionExpectation(Hmm *hmm, int64_t kmerIndex, double meanCurrent, double p);

void continuousPairHmm_setEmissionExpectation(Hmm *hmm, int64_t kmerIndex, double meanCurrent, double p);

double *continuousPairHmm_getEmissionExpectation(Hmm *hmm, int64_t kmerIndex);

double *continuousPairHmm_getEmissionPosterior(Hmm *hmm, int64_t kmerIndex);

void continuousPairHmm_loadIntoStateMachine(StateMachine *sM, Hmm *hmm);

void continuousPairHmm_normalize(Hmm *hmm);

void continuousPairHmm_randomize(Hmm *hmm);

void continuousPairHmm_destruct(Hmm *hmm);

void continuousPairHmm_dump(Hmm *hmm, FILE *fileHandle);

Hmm *continuousPairHmm_loadFromFile(const char *fileName);

void continuousPairHmm_loadModelFromFile(ContinuousPairHmm *hmm, const char *modelFile);

Hmm *hdpHmm_constructEmpty(double pseudocount, int64_t stateNumber, StateMachineType type, double threshold,
                           void (*addToTransitionExpFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
                           void (*setTransitionFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
                           double (*getTransitionsExpFcn)(Hmm *hmm, int64_t from, int64_t to));

void hdpHmm_loadTransitions(StateMachine *sM, Hmm *hmm);

void hdpHmm_writeToFile(Hmm *hmm, FILE *fileHandle);

Hmm *hdpHmm_loadFromFile(const char *fileName, NanoporeHDP *ndhp);
Hmm *hdpHmm_loadFromFile2(const char *fileName, NanoporeHDP *nHdp);

void hdpHmm_destruct(Hmm *hmm);

// CORE
void hmmContinuous_loadSignalHmm(const char *hmmFile, StateMachine *sM, StateMachineType type);

void hmmContinuous_destruct(Hmm *hmm, StateMachineType type);

Hmm *hmmContinuous_getEmptyHmm(StateMachineType type, double pseudocount, double threshold);

void hmmContinuous_normalize(Hmm *hmm, StateMachineType type);

void hmmContinuous_writeToFile(const char *outFile, Hmm *hmm, StateMachineType type);

int64_t hmmContinuous_howManyAssignments(Hmm *hmm);

#endif

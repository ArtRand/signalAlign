#ifndef DISCRETE_HMM_H_
#define DISCRETE_HMM_H_

#include <stdio.h>
#include <stdint.h>
#include "stateMachine.h"

typedef struct _hmmDiscrete {
    Hmm baseHmm;
    double *transitions;
    double *emissions;
} HmmDiscrete;

// Construct
Hmm *hmmDiscrete_constructEmpty(double pseudocount, int64_t stateNumber, int64_t symbolSetSize, StateMachineType type,
                                void (*addToTransitionExpFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
                                void (*setTransitionFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
                                double (*getTransitionsExpFcn)(Hmm *hmm, int64_t from, int64_t to),
                                void (*addEmissionsExpFcn)(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p),
                                void (*setEmissionExpFcn)(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p),
                                double (*getEmissionExpFcn)(Hmm *hmm, int64_t state, int64_t x, int64_t y),
                                int64_t (*getElementIndexFcn)(void *));

// Transitions
void hmmDiscrete_addToTransitionExpectation(Hmm *hmm, int64_t from, int64_t to, double p);
void hmmDiscrete_setTransitionExpectation(Hmm *hmm, int64_t from, int64_t to, double p);
double hmmDiscrete_getTransitionExpectation(Hmm *hmm, int64_t from, int64_t to);

// Emissions
void hmmDiscrete_addToEmissionExpectation(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p);
void hmmDiscrete_setEmissionExpectation(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p);
double hmmDiscrete_getEmissionExpectation(Hmm *hmm, int64_t state, int64_t x, int64_t y);

// Randomize/Normalize
void hmmDiscrete_randomizeTransitions(Hmm *hmm);
void hmmDiscrete_randomizeEmissions(Hmm *hmm);
void hmmDiscrete_randomize(Hmm *hmmD);
void hmmDiscrete_normalize(Hmm *hmmD);
void hmmDiscrete_normalize2(Hmm *hmm, bool normalizeEmissions);
// Writers
void hmmDiscrete_write(Hmm *hmmD, FILE *fileHandle);

// Loaders
Hmm *hmmDiscrete_loadFromFile(const char *fileName);

// Housekeeping
void hmmDiscrete_destruct(Hmm *hmmD);

// stateMachine interface
StateMachineFunctions *stateMachineFunctions_construct(double (*gapXProbFcn)(const double *, void *),
                                                       double (*gapYProbFcn)(const double *, void *),
                                                       double (*matchProbFcn)(const double *, void *, void *));

// Not yet implemented
//StateMachineFunctions *stateMachineFunctions_constructFromType(int64_t stateMachineType);

#endif

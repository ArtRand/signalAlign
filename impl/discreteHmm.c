#include <stdlib.h>
#include <stdio.h>
#include <stateMachine.h>
#include "bioioC.h"
#include "pairwiseAligner.h"
#include "stateMachine.h"
#include "discreteHmm.h"

// Construct
Hmm *hmmDiscrete_constructEmpty(double pseudocount, int64_t stateNumber, int64_t symbolSetSize, StateMachineType type,
                                void (*addToTransitionExpFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
                                void (*setTransitionFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
                                double (*getTransitionsExpFcn)(Hmm *hmm, int64_t from, int64_t to),
                                void (*addEmissionsExpFcn)(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p),
                                void (*setEmissionExpFcn)(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p),
                                double (*getEmissionExpFcn)(Hmm *hmm, int64_t state, int64_t x, int64_t y),
                                int64_t (*getElementIndexFcn)(void *)) {
    // malloc
    HmmDiscrete *hmmD = st_malloc(sizeof(HmmDiscrete));

    // Set up constants
    hmmD->baseHmm.stateNumber = stateNumber;
    hmmD->baseHmm.symbolSetSize = symbolSetSize;
    hmmD->baseHmm.matrixSize = symbolSetSize*symbolSetSize; // working with symmetric matrices
    hmmD->baseHmm.type = type;

    // Set up transitions matrix
    hmmD->transitions = st_malloc(hmmD->baseHmm.stateNumber * hmmD->baseHmm.stateNumber * sizeof(double));
    for (int64_t i = 0; i < hmmD->baseHmm.stateNumber * hmmD->baseHmm.stateNumber; i++) {
        hmmD->transitions[i] = pseudocount;
    }

    // Set up emissions matrix
    hmmD->emissions = st_malloc(hmmD->baseHmm.stateNumber * hmmD->baseHmm.matrixSize * sizeof(double));
    for (int64_t i = 0; i < hmmD->baseHmm.stateNumber * hmmD->baseHmm.matrixSize; i++) {
        hmmD->emissions[i] = pseudocount;
    }

    // Initialize likelihood
    hmmD->baseHmm.likelihood = 0.0;

    // Set up functions
    // transitions
    hmmD->baseHmm.addToTransitionExpectationFcn = addToTransitionExpFcn; // add
    hmmD->baseHmm.setTransitionFcn = setTransitionFcn;                   // set
    hmmD->baseHmm.getTransitionsExpFcn = getTransitionsExpFcn;           // get
    // emissions
    hmmD->baseHmm.addToEmissionExpectationFcn = addEmissionsExpFcn;      // add
    hmmD->baseHmm.setEmissionExpectationFcn = setEmissionExpFcn;         // set
    hmmD->baseHmm.getEmissionExpFcn = getEmissionExpFcn;                 // get
    // indexing
    hmmD->baseHmm.getElementIndexFcn = getElementIndexFcn;               // indexing

    return (Hmm *) hmmD;
}
// Transitions
void hmmDiscrete_addToTransitionExpectation(Hmm *hmm, int64_t from, int64_t to, double p) {
    HmmDiscrete *hmmD = (HmmDiscrete *) hmm;
    hmmD->transitions[from * hmmD->baseHmm.stateNumber + to] += p;
}

void hmmDiscrete_setTransitionExpectation(Hmm *hmm, int64_t from, int64_t to, double p) {
    HmmDiscrete *hmmD = (HmmDiscrete *) hmm;
    hmmD->transitions[from * hmmD->baseHmm.stateNumber + to] = p;
}

double hmmDiscrete_getTransitionExpectation(Hmm *hmm, int64_t from, int64_t to) {
    HmmDiscrete *hmmD = (HmmDiscrete *) hmm;
    return hmmD->transitions[from * hmmD->baseHmm.stateNumber + to];
}

// Emissions
void hmmDiscrete_addToEmissionExpectation(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p) {
    HmmDiscrete *hmmD = (HmmDiscrete *) hmm;
    int64_t tableIndex = x * hmmD->baseHmm.symbolSetSize + y;
    hmmD->emissions[(state * hmmD->baseHmm.matrixSize) + tableIndex] += p;
}

void hmmDiscrete_setEmissionExpectation(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p) {
    HmmDiscrete *hmmD = (HmmDiscrete *) hmm;
    int64_t tableIndex = x * hmmD->baseHmm.symbolSetSize + y;
    hmmD->emissions[(state * hmmD->baseHmm.matrixSize) + tableIndex] = p;
}

double hmmDiscrete_getEmissionExpectation(Hmm *hmm, int64_t state, int64_t x, int64_t y) {
    HmmDiscrete *hmmD = (HmmDiscrete *) hmm;
    int64_t tableIndex = x * hmmD->baseHmm.symbolSetSize + y;
    return hmmD->emissions[(state * hmmD->baseHmm.matrixSize) + tableIndex];
}

// Randomize/Normalize
void hmmDiscrete_randomizeTransitions(Hmm *hmm) {
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            //hmmDiscrete_setTransitionExpectation((Hmm *)hmmD, from, to, st_random());
            hmm->setTransitionFcn(hmm, from, to, st_random());
        }
    }
}

void hmmDiscrete_randomizeEmissions(Hmm *hmm) {
    if (hmm->symbolSetSize <= 0) {
        st_errAbort("hmmDiscrete_randomizeEmissions: got NULL for symbolSetSize\n");
    }
    for (int64_t state = 0; state < hmm->stateNumber; state++) {
        for (int64_t x = 0; x < hmm->symbolSetSize; x++) {
            for (int64_t y = 0; y < hmm->symbolSetSize; y++) {
                //hmmDiscrete_setEmissionExpectation((Hmm *)hmmD, state, x, y, st_random());
                hmm->setEmissionExpectationFcn(hmm, state, x, y, st_random());
            }
        }
    }
}
void hmmDiscrete_randomize(Hmm *hmm) {
    HmmDiscrete *hmmD = (HmmDiscrete *) hmm;
    // Transitions
    hmmDiscrete_randomizeTransitions(hmm);

    // Emissions
    hmmDiscrete_randomizeEmissions(hmm);

    hmmDiscrete_normalize2((Hmm *)hmmD, TRUE);
}

void hmmDiscrete_normalize2(Hmm *hmm, bool normalizeEmissions) {
    // Transitions
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        double total = 0.0;
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            total += hmm->getTransitionsExpFcn(hmm, from, to);
        }
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            double newProb = hmm->getTransitionsExpFcn(hmm, from, to) / total;
            hmm->setTransitionFcn(hmm, from, to, newProb);
        }
    }
    if (normalizeEmissions) {
        for (int64_t state = 0; state < hmm->stateNumber; state++) {
            double total = 0.0;
            for (int64_t x = 0; x < hmm->symbolSetSize; x++) {
                for (int64_t y = 0; y < hmm->symbolSetSize; y++) {
                    total += hmm->getEmissionExpFcn(hmm, state, x, y);
                }
            }
            for (int64_t x = 0; x < hmm->symbolSetSize; x ++) {
                for (int64_t y = 0; y < hmm->symbolSetSize; y++) {
                    double newProb = hmm->getEmissionExpFcn(hmm, state, x, y) / total;
                    hmm->setEmissionExpectationFcn(hmm, state, x, y, newProb);
                }
            }
        }
    }
}

// writers
void hmmDiscrete_write(Hmm *hmm, FILE *fileHandle) {
    /*
     * Format:
     * type \t stateNumber \t symbolSetSize \n
     * transition \t transition ... likelihood \n
     * emission \t emission \t ... \n
     */
    HmmDiscrete *hmmD = (HmmDiscrete *) hmm;
    // basics
    fprintf(fileHandle, "%i\t", hmmD->baseHmm.type); // 0
    fprintf(fileHandle, "%lld\t", hmmD->baseHmm.stateNumber); // 1
    fprintf(fileHandle, "%lld\t", hmmD->baseHmm.symbolSetSize); // 2
    fprintf(fileHandle, "\n");
    // transitions
    for (int64_t i = 0; i < hmmD->baseHmm.stateNumber * hmmD->baseHmm.stateNumber; i++) {
        fprintf(fileHandle, "%f\t", hmmD->transitions[i]);
    }
    // likelihood
    fprintf(fileHandle, "%f\n", hmmD->baseHmm.likelihood);
    // emissions
    for (int64_t i = 0; i < hmmD->baseHmm.stateNumber * hmmD->baseHmm.matrixSize; i++) {
        fprintf(fileHandle, "%f\t", hmmD->emissions[i]);
    }
    fprintf(fileHandle, "\n");
}

// Loaders
Hmm *hmmDiscrete_loadFromFile(const char *fileName) {
    // open the file
    FILE *fH = fopen(fileName, "r");
    // read the first line of the file and split at whitespace
    char *string = stFile_getLineFromFile(fH);
    stList *tokens = stString_split(string);

    // improper input check
    if (stList_length(tokens) < 2) {
        st_errAbort("Got an empty line in the input state machine file %s\n", fileName); // this is the wrong error message?
    }
    // setup
    int type;
    int64_t stateNumber, parameterSetSize;
    // get StateMachineType
    int64_t j = sscanf(stList_get(tokens, 0), "%i", &type);
    if (j != 1) {
        st_errAbort("Failed to parse type (int) from string: %s\n", string);
    }
    // get stateNumber
    int64_t s = sscanf(stList_get(tokens, 1), "%lld", &stateNumber);
    if (s != 1) {
        st_errAbort("Failed to parse state number (int) from string: %s\n", string);
    }
    // get parameterSetSize
    int64_t n = sscanf(stList_get(tokens, 2), "%lld", &parameterSetSize);
    if (n != 1) {
        st_errAbort("Failed to parse symbol set size (int) from string: %s\n", string);
    }

    // make empty Hmm
    Hmm *hmm = hmmDiscrete_constructEmpty(0.0, stateNumber, parameterSetSize, type,
                                           hmmDiscrete_addToTransitionExpectation,
                                           hmmDiscrete_setTransitionExpectation,
                                           hmmDiscrete_getTransitionExpectation,
                                           hmmDiscrete_addToEmissionExpectation,
                                           hmmDiscrete_setEmissionExpectation,
                                           hmmDiscrete_getEmissionExpectation,
                                           emissions_discrete_getBaseIndex);
    // cleanup setup line
    free(string);
    stList_destruct(tokens);

    // downcast
    HmmDiscrete *hmmD = (HmmDiscrete *) hmm;

    // Transitions
    // parse transitions line
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);

    // check for the correct number of transitions
    if (stList_length(tokens) != hmmD->baseHmm.stateNumber * hmmD->baseHmm.stateNumber + 1) { // + 1 bc. likelihood is also on that line
        st_errAbort(
                "Got the wrong number of transitions in the input state machine file %s, got %" PRIi64 " instead of %" PRIi64 "\n",
                fileName, stList_length(tokens), hmmD->baseHmm.stateNumber * hmmD->baseHmm.stateNumber + 1);
    }
    // load transitions
    for (int64_t i = 0; i < hmmD->baseHmm.stateNumber * hmmD->baseHmm.stateNumber; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(hmmD->transitions[i]));
        if (j != 1) {
            st_errAbort("Failed to parse transition prob (float) from string: %s\n", string);
        }
    }
    // load likelihood
    j = sscanf(stList_get(tokens, stList_length(tokens) - 1), "%lf", &(hmmD->baseHmm.likelihood));
    if (j != 1) {
        st_errAbort("Failed to parse likelihood (float) from string: %s\n", string);
    }
    // Cleanup transitions line
    free(string);
    stList_destruct(tokens);

    // Now parse the emissions line
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);

    // check for the correct number of emissions
    if (stList_length(tokens) != hmmD->baseHmm.stateNumber * hmmD->baseHmm.matrixSize) {
        st_errAbort(
                "Got the wrong number of emissions in the input state machine file %s, got %" PRIi64 " instead of %" PRIi64 "\n",
                fileName, stList_length(tokens), hmmD->baseHmm.matrixSize);
    }

    // load emissions
    for (int64_t i = 0; i < hmmD->baseHmm.stateNumber * hmmD->baseHmm.matrixSize; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(hmmD->emissions[i]));
        if (j != 1) {
            st_errAbort("Failed to parse emission prob (float) from string: %s\n", string);
        }
    }
    // Final cleanup
    free(string);
    stList_destruct(tokens);
    fclose(fH);

    return (Hmm *)hmmD;
}

// Housekeeping
void hmmDiscrete_destruct(Hmm *hmm) {
    HmmDiscrete *hmmD = (HmmDiscrete *) hmm;
    free(hmmD->transitions);
    free(hmmD->emissions);
    free(hmmD);
}

// stateMachine interface
StateMachineFunctions *stateMachineFunctions_construct(double (*gapXProbFcn)(const double *, void *),
                                                       double (*gapYProbFcn)(const double *, void *),
                                                       double (*matchProbFcn)(const double *, void *, void *)) {
    StateMachineFunctions *sMfs = malloc(sizeof(StateMachineFunctions));
    sMfs->gapXProbFcn = gapXProbFcn;
    sMfs->gapYProbFcn = gapYProbFcn;
    sMfs->matchProbFcn = matchProbFcn;
    return sMfs;
}


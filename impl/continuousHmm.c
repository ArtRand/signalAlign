#include <stdio.h>
#include <inttypes.h>
//#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <stateMachine.h>
#include "hdp_math_utils.h"
#include "discreteHmm.h"
#include "pairwiseAligner.h"


static bool hmmContinuous_checkTransitions(double *transitions, int64_t nbTransitions) {
    for (int64_t i = 0; i < nbTransitions; i++) {
        if (isnan(transitions[i])) {
            fprintf(stderr, "GOT NaN TRANS\n");
            return FALSE;
        } else {
            continue;
        }
    }
    return TRUE;
}

static inline Hmm *hmmContinuous_loadSignalHmmFromFile(const char *fileName, StateMachineType type,
                                                       double transitionPseudocount, double emissionPseudocount) {
    if (type == threeState) {
        Hmm *hmm = continuousPairHmm_loadFromFile(fileName, type, transitionPseudocount, emissionPseudocount);
        return hmm;
    }
    if (type == threeStateHdp) {
        Hmm *hmm = hdpHmm_loadFromFile(fileName, type, NULL);
        return hmm;
    }
    return 0;
}

static inline void hmmContinuous_loadIntoStateMachine(StateMachine *sM, Hmm *hmm, StateMachineType type) {
    if (type == threeState) {
        continuousPairHmm_loadTransitionsIntoStateMachine(sM, hmm);
        continuousPairHmm_loadEmissionsIntoStateMachine(sM, hmm);
    }
    if (type == threeStateHdp) {
        continuousPairHmm_loadTransitionsIntoStateMachine(sM, hmm);
    }
}

static void hmmContinuous_loadEventModel(double *stateMachineEventModel, double *hmmEventModel,
                                         int64_t parameterSetSize) {
    for (int64_t i = 0; i < parameterSetSize * MODEL_PARAMS; i++) {
        //int64_t eventMeanIdx = i * MODEL_PARAMS;
        //int64_t eventSdIdx = (i * MODEL_PARAMS) + 1;
        //hmmEventModel[eventMeanIdx] = stateMachineEventModel[eventMeanIdx];
        //hmmEventModel[eventSdIdx] = stateMachineEventModel[eventSdIdx];
        hmmEventModel[i] = stateMachineEventModel[i];
        hmmEventModel[i] = stateMachineEventModel[i];
    }
}

///////////////////////////////////////////// Continuous Pair HMM /////////////////////////////////////////////////////
static double continuousPairHmm_getEventModelMean(ContinuousPairHmm *hmm, int64_t kmerIndex) {
    stateMachine_index_check(kmerIndex);
    int64_t tableIndex = kmerIndex * MODEL_PARAMS;
    return hmm->eventModel[tableIndex];
}

static double continuousPairHmm_getEventModelSD(ContinuousPairHmm *hmm, int64_t kmerIndex) {
    stateMachine_index_check(kmerIndex);
    int64_t tableIndex = (kmerIndex * MODEL_PARAMS) + 1;
    return hmm->eventModel[tableIndex];
}

static double *continuousPairHmm_getEventModelEntry(Hmm *hmm, int64_t kmerIndex) {
    stateMachine_index_check(kmerIndex);
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *)hmm;
    return &(cpHmm->eventModel[kmerIndex * MODEL_PARAMS]);
}

static double continuousPairHmm_descaleEventMean_JordanStyle(ContinuousPairHmm *self,
                                                             double eventMean, double levelMean) {
    // (x  + var * level_mean - scale * level_mean - shift) / var
    return (eventMean + self->var * levelMean - self->scale * levelMean - self->shift) / self->var;
}

Hmm *continuousPairHmm_construct(double transitionPseudocount, double emissionPseudocount,
                                 int64_t stateNumber, int64_t alphabetSize, char *alphabet, int64_t kmerLength,
                                 StateMachineType type, double scale, double shift, double var) {
    if (type != threeState && type != threeStateHdp) {
        st_errAbort("ContinuousPair HMM construct: Wrong HMM type for this function got: %i", type);
    }
    ContinuousPairHmm *cpHmm = st_malloc(sizeof(ContinuousPairHmm));
    cpHmm->baseHmm.type = type;
    cpHmm->baseHmm.stateNumber = stateNumber;
    cpHmm->baseHmm.parameterSetSize = intPow(alphabetSize, kmerLength);
    cpHmm->baseHmm.matrixSize = MODEL_PARAMS;
    cpHmm->baseHmm.likelihood = 0.0;

    // transitions
    cpHmm->baseHmm.addToTransitionExpectationFcn = hmm_addToTransitionsExpectation;
    cpHmm->baseHmm.setTransitionFcn = hmm_setTransitionExpectation;
    cpHmm->baseHmm.getTransitionsExpFcn = hmm_getTransitionExpectation;
    cpHmm->baseHmm.transitions = st_malloc(stateNumber * stateNumber * sizeof(double));
    for (int64_t i = 0; i < (stateNumber * stateNumber); i++) {
        cpHmm->baseHmm.transitions[i] = transitionPseudocount;
    }

    cpHmm->getElementIndexFcn = kmer_id;

    // kmer stuff
    cpHmm->baseHmm.alphabetSize = alphabetSize;
    cpHmm->baseHmm.alphabet = sequence_prepareAlphabet(alphabet, alphabetSize);
    cpHmm->baseHmm.kmerLength = kmerLength;

    // scaling factors
    cpHmm->scale = scale;
    cpHmm->shift = shift;
    cpHmm->var = var;
    cpHmm->getDescaledEvent = emissions_signal_descaleEventMean_JordanStyle;

    // emissions
    cpHmm->addToEmissionExpectationFcn = continuousPairHmm_addToEmissionExpectation;
    cpHmm->setEmissionExpectationFcn = continuousPairHmm_setEmissionExpectation;
    cpHmm->getEmissionExpFcn = continuousPairHmm_getEmissionExpectation;
    cpHmm->getPosteriorExpFcn = continuousPairHmm_getEmissionPosterior;
    cpHmm->getEventModelEntry = continuousPairHmm_getEventModelEntry;
    // matrix to hold event expectations (just the means right now)
    cpHmm->eventExpectations = st_calloc(cpHmm->baseHmm.parameterSetSize * NORMAL_DISTRIBUTION_PARAMS, sizeof(double));

    // matrix to hold the old model
    cpHmm->eventModel = st_calloc(cpHmm->baseHmm.parameterSetSize * MODEL_PARAMS, sizeof(double));

    // matrix to hold the 'weights' of each event expectation
    cpHmm->posteriors = st_calloc(cpHmm->baseHmm.parameterSetSize, sizeof(double));
    for (int64_t i = 0; i < cpHmm->baseHmm.parameterSetSize; i++) {
        cpHmm->posteriors[i] = emissionPseudocount;
    }

    cpHmm->observed = st_malloc(cpHmm->baseHmm.parameterSetSize * sizeof(bool));
    for (int64_t i = 0; i < cpHmm->baseHmm.parameterSetSize; i++) {
        cpHmm->observed[i] = FALSE;
    }

    cpHmm->hasModel = FALSE;

    return (Hmm *)cpHmm;
}
// TODO remove all this needless casting
// transitions
void hmm_addToTransitionsExpectation(Hmm *hmm, int64_t from, int64_t to, double p) {
    hmm->transitions[from * hmm->stateNumber + to] += p;
}

void hmm_setTransitionExpectation(Hmm *hmm, int64_t from, int64_t to, double p) {
    hmm->transitions[from * hmm->stateNumber + to] = p;
}

double hmm_getTransitionExpectation(Hmm *hmm, int64_t from, int64_t to) {
    return hmm->transitions[from * hmm->stateNumber + to];
}
// TODO double check that this is right
void continuousPairHmm_addToEmissionExpectation(Hmm *hmm, int64_t kmerIndex, double meanCurrent, double p) {
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *) hmm;
    stateMachine_index_check(kmerIndex);
    int64_t tableIndex = kmerIndex * NORMAL_DISTRIBUTION_PARAMS;
    cpHmm->eventExpectations[tableIndex] += (p * meanCurrent);  // μ_k
    cpHmm->posteriors[kmerIndex] += p;
    double uK = cpHmm->eventExpectations[tableIndex] / cpHmm->posteriors[kmerIndex];
    cpHmm->eventExpectations[tableIndex + 1] += p * (meanCurrent - uK) * (meanCurrent - uK); // σ_k
    cpHmm->observed[kmerIndex] = TRUE;
}

void continuousPairHmm_setEmissionExpectation(Hmm *hmm, int64_t kmerIndex, double meanCurrent, double p) {
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *) hmm;
    int64_t tableIndex = kmerIndex * NORMAL_DISTRIBUTION_PARAMS;
    cpHmm->eventExpectations[tableIndex] = (p * meanCurrent);  // μ_k
    cpHmm->posteriors[kmerIndex] = p;
    double uK = cpHmm->eventExpectations[tableIndex] / cpHmm->posteriors[kmerIndex];
    cpHmm->eventExpectations[tableIndex + 1] = p * (meanCurrent - uK) * (meanCurrent - uK); // σ_k
    cpHmm->observed[kmerIndex] = TRUE;
}

double *continuousPairHmm_getEmissionExpectation(Hmm *hmm, int64_t kmerIndex) {
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *) hmm;
    int64_t tableIndex = kmerIndex * NORMAL_DISTRIBUTION_PARAMS;
    return &(cpHmm->eventExpectations[tableIndex]);
}

double *continuousPairHmm_getEmissionPosterior(Hmm *hmm, int64_t kmerIndex) {
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *) hmm;
    return &(cpHmm->posteriors[kmerIndex]);
}

void continuousPairHmm_destruct(Hmm *hmm) {
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *) hmm;
    free(cpHmm->baseHmm.transitions);
    free(cpHmm->eventExpectations);
    free(cpHmm->eventModel);
    free(cpHmm->posteriors);
    free(cpHmm);
}

Hmm *continuousPairHmm_makeExpectationsHmm(StateMachine *sM,
                                           double transitionsPseudocount, double emissionsPseudocount) {
    // make a hmm object that 'fits' with the stateMachine add pseudocount
    Hmm *hmm = continuousPairHmm_construct(transitionsPseudocount, emissionsPseudocount,
                                           sM->stateNumber, sM->alphabetSize, sM->alphabet, sM->kmerLength,
                                           sM->type, sM->scale, sM->shift, sM->var);
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *)hmm;
    // load the event model into the HMM
    hmmContinuous_loadEventModel(sM->EMISSION_MATCH_MATRIX, cpHmm->eventModel, cpHmm->baseHmm.parameterSetSize);
    cpHmm->hasModel = TRUE;
    return (Hmm *)cpHmm;
}

void continuousPairHmm_loadModelFromFile(ContinuousPairHmm *hmm, const char *modelFile) {
    /*
     *  the model file has the format:
     *  line 0: alphabetSize, alphabet, kmerLength
     *  line 1: [correlation coefficient] [level_mean] [level_sd] [noise_mean]
     *              [noise_sd] [noise_lambda ](.../kmer) \n
     *  line 2: [correlation coefficient] [level_mean] [level_sd, scaled]
     *              [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
     */
    st_errAbort("Model format has changed (7/18/16). This function should be depreciated\n");
    FILE *fH = fopen(modelFile, "r");

    // Line 0: parse the match emissions line
    char *string = stFile_getLineFromFile(fH);
    stList *tokens = stString_split(string);
    if (stList_length(tokens) != 3) {
        st_errAbort("continuousPairHmm_loadModelFromFile: Model file %s does not have the correct number of stateMachine "
                            "parameters should be 3 got %"PRId64"\n", modelFile, stList_length(tokens));
    }
    int64_t alphabetSize, kmerLength, j;
    j = sscanf(stList_get(tokens, 0), "%"SCNd64, &alphabetSize);
    if (j != 1) {
        st_errAbort("stateMachine3_loadFromFile: error parsing alphabet size\n");
    }
    j = sscanf(stList_get(tokens, 2), "%"SCNd64, &kmerLength);
    if (j != 1) {
        st_errAbort("stateMachine3_loadFromFile: error parsing kmer length\n");
    }
    if (intPow(alphabetSize, kmerLength) != hmm->baseHmm.parameterSetSize) {
        st_errAbort("continuousPairHmm_loadModelFromFile: This model won't fit with this HMM.\n");
    }
    free(string);
    stList_destruct(tokens);
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);

    // Line 1: emissions match matrix
    // check to make sure that the model will fit in the hmm
    if (stList_length(tokens) != 1 + (hmm->baseHmm.parameterSetSize * MODEL_PARAMS)) {
        st_errAbort("cpHMM_setEmissionsDefaults: incorrect for signal model (match emissions) got %lld, should be %llf\n",
                    stList_length(tokens), 1 + (hmm->baseHmm.parameterSetSize * MODEL_PARAMS));
    }
    // load means into the event expectations
    int64_t h, l;
    for (int64_t i = 0; i < hmm->baseHmm.parameterSetSize; i++) {
        h = sscanf(stList_get(tokens, (1 + (i * MODEL_PARAMS))), "%lf",
                   &(hmm->eventModel[i * NORMAL_DISTRIBUTION_PARAMS]));

        l = sscanf(stList_get(tokens, (1 + (i * MODEL_PARAMS + 1))), "%lf",
                   &(hmm->eventModel[i * NORMAL_DISTRIBUTION_PARAMS + 1]));

        if (h != 1) {
            st_errAbort("cpHmm_setEmissionsToDefaults: error loading pore model (event model)\n");
        }

        if (l != 1) {
            st_errAbort("cpHmm_setEmissionsToDefaults: error loading pore model (event SDs)\n");
        }
    }

    // clean up match emissions line
    free(string);
    stList_destruct(tokens);

    // close file
    fclose(fH);
    hmm->hasModel = TRUE;
}

void continuousPairHmm_normalize(Hmm *hmm) {
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *)hmm;
    if (cpHmm->baseHmm.type != threeState) {
        st_errAbort("continuousPairHmm_normalize: got invalid HMM type: %lld", hmm->type);
    }
    // normalize transitions
    hmmDiscrete_normalizeTransitions((Hmm *)cpHmm);

    // u_k = Sigma((p * obs)) / Sigma(p)
    // o_k = Sigma(p * (obs - model)^2) / Sigma(p) = Sigma(pu) / Sigma(p)
    for (int64_t i = 0; i < cpHmm->baseHmm.parameterSetSize; i++) {
        if (cpHmm->observed[i]) {
            int64_t tableIndex = i * MODEL_PARAMS;

            double sigmaKmerPosteriorProb = *(cpHmm->getPosteriorExpFcn((Hmm *)cpHmm, i));
            double u_k = *(cpHmm->getEmissionExpFcn((Hmm *)cpHmm, i)) / sigmaKmerPosteriorProb;
            double o_k = sqrt((*(cpHmm->getEmissionExpFcn((Hmm *)cpHmm, i) + 1)) / sigmaKmerPosteriorProb);

            double prevMean = cpHmm->eventModel[tableIndex];
            double prevSD = cpHmm->eventModel[tableIndex + 1];
            cpHmm->eventModel[tableIndex] = u_k == 0.0 ? prevMean : u_k;
            cpHmm->eventModel[tableIndex + 1] = o_k == 0.0 ? prevSD : o_k;
        } else {
            continue;
        }
    }
}

void continuousPairHmm_randomize(Hmm *hmm) {
    // set all the transitions to random numbers
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            hmm->setTransitionFcn(hmm, from, to, st_random());
        }
    }
    hmmDiscrete_normalizeTransitions(hmm);
}

void continuousPairHmm_loadTransitionsIntoStateMachine(StateMachine *sM, Hmm *hmm) {
    StateMachine3 *sM3 = (StateMachine3 *)sM;
    // load transitions
    // from match
    sM3->TRANSITION_MATCH_CONTINUE = log(hmm->getTransitionsExpFcn(hmm, match, match)); // stride
    sM3->TRANSITION_GAP_OPEN_X = log(hmm->getTransitionsExpFcn(hmm, match, shortGapX)); // skip
    sM3->TRANSITION_GAP_OPEN_Y = log(hmm->getTransitionsExpFcn(hmm, match, shortGapY));

    // from shortGapX (kmer skip)
    sM3->TRANSITION_MATCH_FROM_GAP_X = log(hmm->getTransitionsExpFcn(hmm, shortGapX, match));
    sM3->TRANSITION_GAP_EXTEND_X = log(hmm->getTransitionsExpFcn(hmm, shortGapX, shortGapX));
    sM3->TRANSITION_GAP_SWITCH_TO_Y = LOG_ZERO;

    // from shortGapY (extra event)
    sM3->TRANSITION_MATCH_FROM_GAP_Y = log(hmm->getTransitionsExpFcn(hmm, shortGapY, match));
    sM3->TRANSITION_GAP_EXTEND_Y = log(hmm->getTransitionsExpFcn(hmm, shortGapY, shortGapY));
    sM3->TRANSITION_GAP_SWITCH_TO_X = log(hmm->getTransitionsExpFcn(hmm, shortGapY, shortGapX));
    //sM3->TRANSITION_GAP_SWITCH_TO_X = LOG_ZERO;
}

void continuousPairHmm_loadEmissionsIntoStateMachine(StateMachine *sM, Hmm *hmm) {
    StateMachine3 *sM3 = (StateMachine3 *)sM;
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *)hmm;
    // load emissions
    for (int64_t i = 0; i < hmm->parameterSetSize; i++) {
        sM3->model.EMISSION_MATCH_MATRIX[(i * MODEL_PARAMS)] = continuousPairHmm_getEventModelMean(cpHmm, i);
        sM3->model.EMISSION_MATCH_MATRIX[(i * MODEL_PARAMS + 1)] = continuousPairHmm_getEventModelSD(cpHmm, i);
        sM3->model.EMISSION_GAP_Y_MATRIX[(i + MODEL_PARAMS)] = continuousPairHmm_getEventModelMean(cpHmm, i);
        sM3->model.EMISSION_GAP_Y_MATRIX[(i * MODEL_PARAMS + 1)] = (continuousPairHmm_getEventModelSD(cpHmm, i)
                                                                        * 1.75);
    }
}

void continuousPairHmm_writeToFile(Hmm *hmm, FILE *fileHandle) {
    /*
     *  the model file has the format:
     *  line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
     *  line 1: match->match \t match->gapX \t match->gapY \t
     *          gapX->match \t gapX->gapX \t gapX->gapY \t
     *          gapY->match \t gapY->gapX \t gapY->gapY \n
     *  line 1: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
     */
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *)hmm;

    int64_t nb_transitions = (cpHmm->baseHmm.stateNumber * cpHmm->baseHmm.stateNumber);

    bool check = hmmContinuous_checkTransitions(cpHmm->baseHmm.transitions, nb_transitions);
    if (check) {
        // write the basic stuff to disk (positions are line:item#, not line:col)
        fprintf(fileHandle, "%"PRId64"\t", cpHmm->baseHmm.stateNumber); // stateNumber 0:0
        fprintf(fileHandle, "%"PRId64"\t", cpHmm->baseHmm.alphabetSize); // alphabetSize 0:1
        fprintf(fileHandle, "%s\t", cpHmm->baseHmm.alphabet); // alphabet
        fprintf(fileHandle, "%"PRId64"\t", cpHmm->baseHmm.kmerLength); // parameterSetSize 0:2
        fprintf(fileHandle, "\n");

        // write out transitions
        for (int64_t i = 0; i < nb_transitions; i++) {
            fprintf(fileHandle, "%f\t", cpHmm->baseHmm.transitions[i]); // transitions 1:(0-9)
        }

        // write the likelihood
        fprintf(fileHandle, "%f\n", cpHmm->baseHmm.likelihood); // likelihood 1:10, newLine

        // write out event model
        for (int64_t i = 0; i < (cpHmm->baseHmm.parameterSetSize * MODEL_PARAMS); i++) {
            fprintf(fileHandle, "%lf\t", cpHmm->eventModel[i]);
        }
        fprintf(fileHandle, "\n"); // newLine

        // write out event expectations
        for (int64_t i = 0; i < (cpHmm->baseHmm.parameterSetSize * NORMAL_DISTRIBUTION_PARAMS); i++) {
            fprintf(fileHandle, "%lf\t", cpHmm->eventExpectations[i]);
        }
        fprintf(fileHandle, "\n");

        // write out event posteriors
        for (int64_t i = 0; i < cpHmm->baseHmm.parameterSetSize; i++) {
            fprintf(fileHandle, "%lf\t", cpHmm->posteriors[i]);
        }
        fprintf(fileHandle, "\n");

        // write out observation mask
        for (int64_t i = 0; i < cpHmm->baseHmm.parameterSetSize; i++) {
            fprintf(fileHandle, "%d\t", cpHmm->observed[i]);
        }
        fprintf(fileHandle, "\n");
    }
}

Hmm *continuousPairHmm_loadFromFile(const char *fileName, StateMachineType type,
                                    double transitionPseudocount, double emissionPseudocount) {
    /*
     *  the model file has the format:
     *  line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
     *  line 1: match->match \t match->gapX \t match->gapY \t
     *          gapX->match \t gapX->gapX \t gapX->gapY \t
     *          gapY->match \t gapY->gapX \t gapY->gapY \n
     *  line 1: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
     */
    FILE *fH = fopen(fileName, "r");

    // line 0
    char *string = stFile_getLineFromFile(fH);
    stList *tokens = stString_split(string);

    if (stList_length(tokens) != 4) {
        st_errAbort("cpHmm_loadFromFile: Got incorrect header line, has %lld tokens should have 4\n",
                    stList_length(tokens));
    }

    int64_t stateNumber, alphabetSize, kmerLength;
    char *alphabet;

    int64_t j = sscanf(stList_get(tokens, 0), "%"SCNd64, &stateNumber); // type
    if (j != 1) {
        st_errAbort("Failed to parse stateNumber (int) from string: %s\n", string);
    }
    j = sscanf(stList_get(tokens, 1), "%"SCNd64, &alphabetSize); // stateNumber
    if (j != 1) {
        st_errAbort("Failed to parse state number (int) from string: %s\n", string);
    }
    alphabet = (char *)stList_get(tokens, 2);

    j = sscanf(stList_get(tokens, 3), "%"SCNd64, &kmerLength); // parameterSetSize
    if (j != 1) {
        st_errAbort("Failed to parse symbol set size (int) from string: %s\n", string);
    }

    // make empty cpHMM
    Hmm *hmm = continuousPairHmm_construct(transitionPseudocount, emissionPseudocount,
                                           stateNumber, alphabetSize, alphabet, kmerLength, type,
                                           0.0, 0.0, 0.0);

    // Downcast
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *)hmm;

    // cleanup
    free(string);
    stList_destruct(tokens);

    // Transitions
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);

    int64_t nb_transitions = (cpHmm->baseHmm.stateNumber * cpHmm->baseHmm.stateNumber);

    // check for the correct number of transitions
    if (stList_length(tokens) != nb_transitions + 1) { // + 1 bc. likelihood is also on that line
        st_errAbort(
                "Incorrect number of transitions in the input HMM file %s, got %" PRIi64 " instead of %" PRIi64 "\n",
                fileName, stList_length(tokens), nb_transitions + 1);
    }
    // load them
    for (int64_t i = 0; i < nb_transitions; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(cpHmm->baseHmm.transitions[i]));
        if (j != 1) {
            st_errAbort("Failed to parse transition prob (float) from string: %s\n", string);
        }
    }
    // load likelihood
    j = sscanf(stList_get(tokens, stList_length(tokens) - 1), "%lf", &(cpHmm->baseHmm.likelihood));
    if (j != 1) {
        st_errAbort("Failed to parse likelihood (float) from string: %s\n", string);
    }
    // Cleanup transitions line
    free(string);
    stList_destruct(tokens);

    // Event Model
    string = stFile_getLineFromFile(fH);
    if (string == NULL) {
        st_errAbort("cpHmm_loadFromFile: Got to end of file when looking for event model\n");
    }
    tokens = stString_split(string);
    if (stList_length(tokens) != cpHmm->baseHmm.parameterSetSize * MODEL_PARAMS) {
        st_errAbort("cpHmm_loadFromFile: Didn't find the correct number of tokens for event model\n");
    }

    for (int64_t i = 0; i < cpHmm->baseHmm.parameterSetSize * MODEL_PARAMS; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(cpHmm->eventModel[i]));
        if (j != 1) {
            st_errAbort("cpHmm_loadFromFile: error parsing event model for event %lld\n", i);
        }
    }
    free(string);
    stList_destruct(tokens);
    return (Hmm *)cpHmm;
}

/////////////////////////////////////////////////// HDP HMM  //////////////////////////////////////////////////////////
static void hdpHmm_addToAssignment(HdpHmm *self, void *kmer, void *event) {
    stList_append(self->kmerAssignments, kmer);
    stList_append(self->eventAssignments, event);
    self->numberOfAssignments += 1;
}

static bool hdpHmm_checkAssignments(HdpHmm *hdpHmm) {
    int64_t nb_kmerAssignments = stList_length(hdpHmm->kmerAssignments);
    int64_t nb_eventAssignments = stList_length(hdpHmm->eventAssignments);
    return nb_eventAssignments == nb_kmerAssignments;
}

Hmm *hdpHmm_constructEmpty(double transitionPseudocount,
                           int64_t stateNumber, int64_t alphabetSize, char *alphabet, int64_t kmerLength,
                           StateMachineType type,
                           double threshold, double scale, double shift, double var) {
    HdpHmm *hmm = st_malloc(sizeof(HdpHmm));
    // setup base Hmm
    hmm->baseHmm.type = type;
    hmm->baseHmm.stateNumber = stateNumber;
    hmm->baseHmm.parameterSetSize = intPow(alphabetSize, kmerLength);
    hmm->baseHmm.matrixSize = MODEL_PARAMS;
    hmm->baseHmm.likelihood = 0.0;

    // transitions
    hmm->baseHmm.addToTransitionExpectationFcn = hmm_addToTransitionsExpectation;   // add
    hmm->baseHmm.setTransitionFcn = hmm_setTransitionExpectation;                   // set
    hmm->baseHmm.getTransitionsExpFcn = hmm_getTransitionExpectation;               // get

    // kmer stuff
    hmm->baseHmm.alphabetSize = alphabetSize;
    hmm->baseHmm.alphabet = sequence_prepareAlphabet(alphabet, alphabetSize);
    hmm->baseHmm.kmerLength = kmerLength;

    // scaling parameters
    hmm->scale = scale;
    hmm->shift = shift;
    hmm->var = var;

    // transitions
    hmm->baseHmm.transitions = st_malloc(stateNumber * stateNumber * sizeof(double));
    for (int64_t i = 0; i < (stateNumber * stateNumber); i++) {
        hmm->baseHmm.transitions[i] = transitionPseudocount;
    }

    // event model (for rescaling events)
    hmm->eventModel = st_calloc(hmm->baseHmm.parameterSetSize * MODEL_PARAMS, sizeof(double));

    // HDP specific stuff
    hmm->threshold = threshold;  // threshold for assignments (must be >= threshold to be an assignment)
    hmm->getElementIndexFcn = kmer_id;
    hmm->getDescaledEvent = emissions_signal_descaleEventMean_JordanStyle;
    hmm->addToAssignments = hdpHmm_addToAssignment;  // function to add to assignment tally
    hmm->kmerAssignments = stList_construct3(0, &free);  // list of kmers that are assigned
    hmm->eventAssignments = stList_construct3(0, &free);  // to the list of events
    hmm->numberOfAssignments = 0;  // total number of assignments
    hmm->nhdp = NULL;  // initialized to NULL

    return (Hmm *)hmm;
}

void hdpHmm_writeToFile(Hmm *hmm, FILE *fileHandle) {
    /*
     *  the model file has the format:
     *  line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
     *  line 1: match->match \t match->gapX \t match->gapY \t
     *          gapX->match \t gapX->gapX \t gapX->gapY \t
     *          gapY->match \t gapY->gapX \t gapY->gapY \n
     *  line 1: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
     */
    HdpHmm *hdpHmm = (HdpHmm *)hmm;

    int64_t nb_transitions = (hdpHmm->baseHmm.stateNumber * hdpHmm->baseHmm.stateNumber);

    bool transitionCheck = hmmContinuous_checkTransitions(hdpHmm->baseHmm.transitions, nb_transitions);
    bool assignmentCheck = hdpHmm_checkAssignments(hdpHmm);
    if (transitionCheck && assignmentCheck) {
        // write the basic stuff to disk (positions are line:item#, not line:col)
        fprintf(fileHandle, "%"PRId64"\t", hdpHmm->baseHmm.stateNumber); // stateNumber 0:0
        fprintf(fileHandle, "%"PRId64"\t", hdpHmm->baseHmm.alphabetSize); // alphabetSize 0:1
        fprintf(fileHandle, "%s\t", hdpHmm->baseHmm.alphabet); // alphabet
        fprintf(fileHandle, "%"PRId64"\t", hdpHmm->baseHmm.kmerLength); // parameterSetSize 0:2
        fprintf(fileHandle, "\n");

        // write out transitions
        for (int64_t i = 0; i < nb_transitions; i++) {
            // transitions 1:(0-9)
            fprintf(fileHandle, "%f\t", hdpHmm->baseHmm.transitions[i]);
        }

        // likelihood 1:10, newLine
        fprintf(fileHandle, "%f\n", hdpHmm->baseHmm.likelihood);

        // write out event model
        for (int64_t i = 0; i < (hdpHmm->baseHmm.parameterSetSize * MODEL_PARAMS); i++) {
            fprintf(fileHandle, "%lf\t", hdpHmm->eventModel[i]);
        }
        fprintf(fileHandle, "\n"); // newLine

        // write out the assignment events
        for (int64_t i = 0; i < hdpHmm->numberOfAssignments; i++) {
            double *meanCurrent = stList_get(hdpHmm->eventAssignments, i);
            fprintf(fileHandle, "%lf\t", *meanCurrent);
        }
        fprintf(fileHandle, "\n"); // newLine

        // write out the assignment kmers
        for (int64_t i = 0; i < hdpHmm->numberOfAssignments; i++) {
            // get the starting position in the kmer sequence
            char *kmer = stList_get(hdpHmm->kmerAssignments, i);
            for (int64_t n = 0; n < hdpHmm->baseHmm.kmerLength; n++) {
                fprintf(fileHandle, "%c", kmer[n]);
            }
            //fprintf(fileHandle, " "); // space
            fprintf(fileHandle, "\t"); // tab?
        }
        fprintf(fileHandle, "\n"); // newLine
    }
}

Hmm *hdpHmm_loadFromFile(const char *fileName, StateMachineType type, NanoporeHDP *nHdp) {
    /*
     *  the model file has the format:
     *  line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
     *  line 1: match->match \t match->gapX \t match->gapY \t
     *          gapX->match \t gapX->gapX \t gapX->gapY \t
     *          gapY->match \t gapY->gapX \t gapY->gapY \n
     *  line 1: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n

     */
    FILE *fH = fopen(fileName, "r");

    // line 0
    char *string = stFile_getLineFromFile(fH);
    stList *tokens = stString_split(string);

    // should have type, stateNumber, threshold, and number of assignments
    if (stList_length(tokens) != 4) {
        st_errAbort("ERROR loading hdpHmm, got %lld tokens should get 4\nLine is %s", stList_length(tokens), string);

    }

    int64_t stateNumber, alphabetSize, kmerLength;
    char *alphabet;

    int64_t j = sscanf(stList_get(tokens, 0), "%"SCNd64, &stateNumber); // type
    if (j != 1) {
        st_errAbort("Failed to parse stateNumber (int) from string: %s\n", string);
    }
    j = sscanf(stList_get(tokens, 1), "%"SCNd64, &alphabetSize); // stateNumber
    if (j != 1) {
        st_errAbort("Failed to parse state number (int) from string: %s\n", string);
    }
    alphabet = (char *)stList_get(tokens, 2);  // alphabet

    j = sscanf(stList_get(tokens, 3), "%"SCNd64, &kmerLength); // kmerLength
    if (j != 1) {
        st_errAbort("Failed to parse number of assignments (int) from string: %s\n", string);
    }

    // make empty hdpHmm object
    Hmm *hmm = hdpHmm_constructEmpty(0.0, stateNumber, alphabetSize, alphabet, kmerLength, type,
                                     0.0, 0.0, 0.0, 0.0);

    HdpHmm *hdpHmm = (HdpHmm *)hmm;  // downcast
    hdpHmm->nhdp = nHdp;

    // cleanup
    free(string);
    stList_destruct(tokens);

    // transitions
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);

    int64_t nb_transitions = (hdpHmm->baseHmm.stateNumber * hdpHmm->baseHmm.stateNumber);

    // check for the correct number of transitions
    if (stList_length(tokens) != nb_transitions + 1) { // + 1 bc. likelihood is also on that line
        st_errAbort(
                "Incorrect number of transitions in the input HMM file %s, got %" PRIi64 " instead of %" PRIi64 "\n",
                fileName, stList_length(tokens), nb_transitions + 1);
    }

    // load them
    for (int64_t i = 0; i < nb_transitions; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(hdpHmm->baseHmm.transitions[i]));
        if (j != 1) {
            st_errAbort("Failed to parse transition prob (float) from string: %s\n", string);
        }
    }

    // load likelihood
    j = sscanf(stList_get(tokens, stList_length(tokens) - 1), "%lf",
               &(hdpHmm->baseHmm.likelihood));
    if (j != 1) {
        st_errAbort("Failed to parse likelihood (float) from string: %s\n", string);
    }

    // Cleanup transitions line
    free(string);
    stList_destruct(tokens);

    // Event Model
    string = stFile_getLineFromFile(fH);
    if (string == NULL) {
        st_errAbort("hdpHmm_loadFromFile: Got to end of file when looking for event model\n");
    }
    tokens = stString_split(string);
    if (stList_length(tokens) != hdpHmm->baseHmm.parameterSetSize * MODEL_PARAMS) {
        st_errAbort("hdpHmm_loadFromFile: Didn't find the correct number of tokens for event model\n");
    }

    for (int64_t i = 0; i < hdpHmm->baseHmm.parameterSetSize * MODEL_PARAMS; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(hdpHmm->eventModel[i]));
        if (j != 1) {
            st_errAbort("hdpHmm_loadFromFile: error parsing event model for event %lld\n", i);
        }
    }
    free(string);
    stList_destruct(tokens);

    // load the assignments into the Nanopore HDP [this is basically the same as Jordan's code for updating from
    // an alignment]
    if (nHdp != NULL) {
        // make stLists to hold things temporarily
        stList *signalList = stList_construct3(0, &free);

        // parse the events (current means)
        string = stFile_getLineFromFile(fH);
        tokens = stString_split(string);

        int64_t nb_assignments = stList_length(tokens);
        hdpHmm->numberOfAssignments = nb_assignments;

        // parse the events into a list
        char *signal_str;
        for (int64_t i = 0; i < hdpHmm->numberOfAssignments; i++) {
            signal_str = (char *)stList_get(tokens, i);
            double *signal_ptr = (double *)st_malloc(sizeof(double));
            sscanf(signal_str, "%lf", signal_ptr);
            stList_append(signalList, signal_ptr);
        }

        // cleanup event line
        free(string);
        stList_destruct(tokens);

        // parse the kmer assignment line
        string = stFile_getLineFromFile(fH);
        tokens = stString_split(string);
        if (stList_length(tokens) != hdpHmm->numberOfAssignments) {
            st_errAbort("Incorrect number of events got %lld, should be %lld\n",
                        stList_length(tokens), hdpHmm->numberOfAssignments);
        }
        char *assignedKmer;
        int64_t *dp_ids = st_malloc(sizeof(int64_t) * hdpHmm->numberOfAssignments);
        for (int64_t i = 0; i < hdpHmm->numberOfAssignments; i++) {
            assignedKmer = (char *)stList_get(tokens, i);
            dp_ids[i] = kmer_id(assignedKmer, nHdp->alphabet, nHdp->alphabet_size, nHdp->kmer_length);
        }
        // cleanup
        free(string);
        stList_destruct(tokens);

        // convert to arrays
        int64_t dataLength;
        double *signal = stList_toDoublePtr(signalList, &dataLength);
        stList_destruct(signalList);
        reset_hdp_data(nHdp->hdp);
        pass_data_to_hdp(nHdp->hdp, signal, dp_ids, dataLength);
    }
    // close file
    fclose(fH);
    return (Hmm *)hdpHmm;
}

Hmm *hdpHmm_makeExpectationsHmm(StateMachine *sM, double threshold, double transitionsPseudocount) {
    Hmm *hmm = hdpHmm_constructEmpty(transitionsPseudocount,
                                     sM->stateNumber, sM->alphabetSize, sM->alphabet, sM->kmerLength,
                                     sM->type, threshold, sM->scale, sM->shift, sM->var);
    HdpHmm *hdpHmm = (HdpHmm *)hmm;
    hmmContinuous_loadEventModel(sM->EMISSION_MATCH_MATRIX, hdpHmm->eventModel, hdpHmm->baseHmm.parameterSetSize);
    return (Hmm *)hdpHmm;
}

void hdpHmm_destruct(Hmm *hmm) {
    HdpHmm *hdpHmm = (HdpHmm *)hmm;
    free(hdpHmm->kmerAssignments);
    free(hdpHmm->eventAssignments);
    free(hdpHmm->baseHmm.transitions);
    // todo check if we need destructers here
    //stList_destruct(hdpHmm->kmerAssignments);
    free(hdpHmm);
}

///////////////////////////////////////////////// CORE FUNCTIONS //////////////////////////////////////////////////////
void hmmContinuous_loadSignalHmm(const char *hmmFile, StateMachine *sM, StateMachineType type, Hmm *expectations) {
    // check for supported types
    if ((type != threeStateHdp) && (type != threeState)) {
        st_errAbort("hmmContinuous_loadSignalHmm - ERROR: got unsupported HMM type %i\n", type);
    }

    // load hmm from file
    Hmm *hmm = hmmContinuous_loadSignalHmmFromFile(hmmFile, type, 0.0, 0.001);
    // probably wont need this in the future. the StateMachine will already have the 'latest' parameters and the
    // expectations object should be initialized to those parameters. this is just here to load the emissions model
    // into the expectations object, which should have already happened.
    if ((expectations != NULL) && (type == threeState)) {
        hmmContinuous_loadModelPrior(hmm, expectations);  // change this to StateMachine function
    }
    // shouldn't need to do this either. the stateMachine should be loaded with the correct parameters
    hmmContinuous_loadIntoStateMachine(sM, hmm, type); // todo remove types here, they're in the hmm object
    hmmContinuous_destruct(hmm, type);
}

void hmmContinuous_loadModelPrior(Hmm *prior, Hmm *expectations) {
    ContinuousPairHmm *priorHmm = (ContinuousPairHmm *)prior;
    ContinuousPairHmm *newHmm = (ContinuousPairHmm *)expectations;
    if (!priorHmm->hasModel) {
        st_errAbort("hmmContinuous_loadModelPrior: Prior HMM needs to have model");
    }
    if (priorHmm->baseHmm.parameterSetSize != newHmm->baseHmm.parameterSetSize) {
        st_errAbort("hmmContinuous_loadModelPrior: These models aren't the same size\n");
    }
    for (int64_t i = 0; i < priorHmm->baseHmm.parameterSetSize * NORMAL_DISTRIBUTION_PARAMS; i++) {
        newHmm->eventModel[i] = priorHmm->eventModel[i];
    }
    newHmm->hasModel = TRUE;
}

Hmm *hmmContinuous_getExpectationsHmm(StateMachine *sM, double threshold,
                                      double transitionsPseudocount, double emissionsPseudocount) {
    if ((sM->type != threeStateHdp) && (sM->type != threeState)) {
        st_errAbort("hmmContinuous_getExpectationsHmm - ERROR: got unsupported HMM type %i\n", sM->type);
    }
    if (sM->type == threeState) {
        Hmm *hmm = continuousPairHmm_makeExpectationsHmm(sM, transitionsPseudocount, emissionsPseudocount);
        return hmm;
    }
    if (sM->type == threeStateHdp) {
        Hmm *hmm = hdpHmm_makeExpectationsHmm(sM, threshold, transitionsPseudocount);
        return hmm;
    }
    return 0;
}
int64_t hmmContinuous_howManyAssignments(Hmm *hmm) {
    if (hmm->type != threeStateHdp) {
        st_errAbort("hmmContinuous: this type of Hmm doesn't have assignments got type: %lld", hmm->type);
    }
    HdpHmm *hdpHmm = (HdpHmm *)hmm;
    return hdpHmm->numberOfAssignments;
}

void hmmContinuous_normalize(Hmm *hmm, StateMachineType type) {
    assert((type == vanilla) || (type == threeState));
    if (type == threeState) {
        continuousPairHmm_normalize(hmm);
    }
}

void hmmContinuous_writeToFile(const char *outFile, Hmm *hmm, StateMachineType type) {
    if ((type != threeStateHdp) && (type != threeState)) {
        st_errAbort("hmmContinuous_writeToFile - ERROR: got unsupported HMM type %i\n", type);
    }
    FILE *fH = fopen(outFile, "w");
    if (type == threeState) {
        continuousPairHmm_writeToFile(hmm, fH);
    }
    if (type == threeStateHdp) {
        hdpHmm_writeToFile(hmm, fH);
    }
    fclose(fH);
}
void hmmContinuous_destruct(Hmm *hmm, StateMachineType type) {
    if (type == threeState) {
        continuousPairHmm_destruct(hmm);
    }
    if (type == threeStateHdp) {
        hdpHmm_destruct(hmm);
    }
}
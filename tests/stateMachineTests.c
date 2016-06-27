// Tests for stateMachine and alignments

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include "stateMachine.h"
#include "CuTest.h"
#include "pairwiseAligner.h"
#include "discreteHmm.h"
#include "continuousHmm.h"
#include "randomSequences.h"

Sequence *getTestReferenceSequence() {
    char *ZymoReferenceFilePath = stString_print("../../signalAlign/tests/test_npReads/ZymoRef.txt");
    FILE *fH = fopen(ZymoReferenceFilePath, "r");
    char *ZymoReferenceSeq = stFile_getLineFromFile(fH);
    int64_t lX = sequence_correctSeqLength(strlen(ZymoReferenceSeq), event);
    Sequence *refSeq = sequence_constructKmerSequence(lX, ZymoReferenceSeq,
                                                      sequence_getKmer, sequence_sliceNucleotideSequence,
            //sequence_getKmerWithBoundsCheck, sequence_sliceNucleotideSequence,
                                                      "CEO", 3, kmer);
    free(ZymoReferenceFilePath);
    return refSeq;
}

double absPercentDiff(double obs, double exp) {
    double percentDiff = ((obs - exp) / exp) * 100;
    return percentDiff > 0 ? percentDiff : -percentDiff;
}

NanoporeRead *loadTestNanoporeRead() {
    char *npReadFile = stString_print("../../signalAlign/tests/test_npReads/ZymoC_ch_1_file1.npRead");
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);
    free(npReadFile);
    return npRead;
}


StateMachine *loadScaledStateMachine3(NanoporeRead *npRead) {
    char *templateModelFile = stString_print("../../signalAlign/models/testModel_template.model");
    StateMachine *sMt = getStateMachine3(templateModelFile);

    // scale model
    emissions_signal_scaleModel(sMt, npRead->templateParams.scale, npRead->templateParams.shift,
                                npRead->templateParams.var, npRead->templateParams.scale_sd,
                                npRead->templateParams.var_sd);
    free(templateModelFile);
    return sMt;
}

StateMachine *loadDescaledStateMachine3(NanoporeRead *npRead) {
    // load stateMachine from model file
    char *templateModelFile = stString_print("../../signalAlign/models/testModel_template.model");
    StateMachine *sM = getStateMachine3_descaled(templateModelFile, npRead->templateParams);
    free(templateModelFile);
    return sM;
}

StateMachine *loadStateMachineHdp(NanoporeRead *npRead) {
    char *modelFile = stString_print("../../signalAlign/models/testModel_template.model");
    NanoporeHDP *nHdp = deserialize_nhdp("../../signalAlign/models/templateSingleLevelFixed.nhdp");
    StateMachine *sM = getHdpStateMachine(nHdp, modelFile, npRead->templateParams);
    free(modelFile);
    return sM;
}

stList *getRemappedAnchors(Sequence *ref, NanoporeRead *npRead, PairwiseAlignmentParameters *p) {
    stList *anchorPairs = getBlastPairsForPairwiseAlignmentParameters((char *)ref->elements, npRead->twoDread, p);

    // remap and filter
    stList *remappedAnchors = nanopore_remapAnchorPairs(anchorPairs, npRead->templateEventMap);
    stList *filteredRemappedAnchors = filterToRemoveOverlap(remappedAnchors);
    return filteredRemappedAnchors;
}

Sequence *replaceBasesInSequence(Sequence *orig, const char *toReplace, const char *replacement) {
    Sequence *copy = sequence_deepCopyNucleotideSequence(orig);
    char *newElements = stString_replace(copy->elements, toReplace, replacement);
    free(copy->elements);
    copy->elements = newElements;
    return copy;
}

static void checkAlignedPairs(CuTest *testCase, stList *blastPairs, int64_t lX, int64_t lY) {
    st_logInfo("I got %" PRIi64 " pairs to check\n", stList_length(blastPairs));
    stSortedSet *pairs = stSortedSet_construct3((int (*)(const void *, const void *)) stIntTuple_cmpFn,
                                                (void (*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(blastPairs); i++) {
        stIntTuple *j = stList_get(blastPairs, i);
        CuAssertTrue(testCase, stIntTuple_length(j) == 4);

        int64_t x = stIntTuple_get(j, 1);
        int64_t y = stIntTuple_get(j, 2);
        int64_t score = stIntTuple_get(j, 0);
        CuAssertTrue(testCase, score > 0);
        CuAssertTrue(testCase, score <= PAIR_ALIGNMENT_PROB_1);

        CuAssertTrue(testCase, x >= 0);
        CuAssertTrue(testCase, y >= 0);
        CuAssertTrue(testCase, x < lX);
        CuAssertTrue(testCase, y < lY);

        //Check is unique
        stIntTuple *pair = stIntTuple_construct2(x, y);
        CuAssertTrue(testCase, stSortedSet_search(pairs, pair) == NULL);
        stSortedSet_insert(pairs, pair);
    }
    stSortedSet_destruct(pairs);
}

// same as above except allows for one event to be aligned to multiple positions
static void checkAlignedPairsWithOverlap(CuTest *testCase, stList *blastPairs, int64_t lX, int64_t lY) {
    st_logInfo("I got %" PRIi64 " pairs to check\n", stList_length(blastPairs));
    stSortedSet *pairs = stSortedSet_construct3((int (*)(const void *, const void *)) stIntTuple_cmpFn,
                                                (void (*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(blastPairs); i++) {
        stIntTuple *j = stList_get(blastPairs, i);
        CuAssertTrue(testCase, stIntTuple_length(j) == 4);

        int64_t x = stIntTuple_get(j, 1);
        int64_t y = stIntTuple_get(j, 2);
        int64_t score = stIntTuple_get(j, 0);
        CuAssertTrue(testCase, score > 0);
        CuAssertTrue(testCase, score <= PAIR_ALIGNMENT_PROB_1);

        CuAssertTrue(testCase, x >= 0);
        CuAssertTrue(testCase, y >= 0);
        CuAssertTrue(testCase, x < lX);
        CuAssertTrue(testCase, y < lY);

        stIntTuple *pair = stIntTuple_construct2(x, y);
        stSortedSet_insert(pairs, pair);
    }
    stSortedSet_destruct(pairs);
}

void updateStateMachineHDP(const char *expectationsFile, StateMachine *sM) {
    StateMachine3_HDP *sM3Hdp = (StateMachine3_HDP *)sM;
    Hmm *transitionsExpectations = hdpHmm_loadFromFile(expectationsFile, sM3Hdp->hdpModel);
    continuousPairHmm_loadTransitionsIntoStateMachine((StateMachine *) sM3Hdp, transitionsExpectations);
    hdpHmm_destruct(transitionsExpectations);
}

void implantModelFromStateMachineIntoHmm(StateMachine *sM, ContinuousPairHmm *hmm) {
    for (int64_t i = 0; i < hmm->baseHmm.symbolSetSize; i++) {
        hmm->eventModel[i * NORMAL_DISTRIBUTION_PARAMS] = sM->EMISSION_MATCH_MATRIX[1 + (i * MODEL_PARAMS)];
        hmm->eventModel[1 + (1 * NORMAL_DISTRIBUTION_PARAMS)] = sM->EMISSION_MATCH_MATRIX[2 + (i * MODEL_PARAMS)];
    }
    hmm->hasModel = TRUE;
}

static void test_loadPoreModel(CuTest *testCase) {
    char *tempFile = stString_print("./tempModel.model");
    CuAssertTrue(testCase, !stFile_exists(tempFile));
    FILE *fH = fopen(tempFile, "w");

    int64_t matrixSize = (1 + NUM_OF_KMERS * MODEL_PARAMS);

    for (int64_t i = 0; i < matrixSize; i++) {
        fprintf(fH, "%"PRId64"\t", i);
    }
    fprintf(fH, "\n");
    //for (int64_t i = matrixSize; i > -1; i--) {
    for (int64_t i = 0; i < matrixSize; i++) {
        fprintf(fH, "%"PRId64"\t", i);
    }
    fprintf(fH, "\n");
    fclose(fH);

    StateMachine *sM = getStateMachine3(tempFile);

    for (int64_t i = 0; i < matrixSize; i++) {
        double x = sM->EMISSION_MATCH_MATRIX[i];
        double y = sM->EMISSION_GAP_Y_MATRIX[i];
        CuAssertDblEquals(testCase, x, (double )i, 0.0);
        CuAssertDblEquals(testCase, y, (double )i, 0.0);
    }

    stFile_rmrf(tempFile);
    free(tempFile);
    stateMachine_destruct(sM);
}

static void test_nanoporeScaleParams(CuTest *testCase) {
    NanoporeRead *npRead = loadTestNanoporeRead();
    Sequence *eventSequence = sequence_construct2(npRead->nbTemplateEvents, npRead->templateEvents, sequence_getEvent,
                                                  sequence_sliceEventSequence, event);
    Sequence *referenceSequence = getTestReferenceSequence();
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
    p->constraintDiagonalTrim = 18;
    stList *filteredRemappedAnchors = getRemappedAnchors(referenceSequence, npRead, p);
    stList *map = nanopore_getAnchorKmersToEventsMap(filteredRemappedAnchors, eventSequence->elements,
                                                     referenceSequence->elements);
    st_uglyf("Map has %lld assignments\n", stList_length(map));
    // test the event/kmer tuples
    for (int64_t i = 0; i < stList_length(map); i++) {
        EventKmerTuple *t = stList_get(map, i);
        stIntTuple *pair = stList_get(filteredRemappedAnchors, i);
        int64_t checkKmerIndex = emissions_discrete_getKmerIndexFromKmer(
                referenceSequence->get(referenceSequence->elements, stIntTuple_get(pair, 0)));
        double *event = eventSequence->get(eventSequence->elements, stIntTuple_get(pair, 1));
        double eventMean = *event;
        double eventSd = *(1 + event);
        double eventDuration = *(2 + event);
        CuAssertIntEquals(testCase, checkKmerIndex, t->kmerIndex);
        CuAssertDblEquals(testCase, eventMean, t->eventMean, 0.0);
        CuAssertDblEquals(testCase, eventSd, t->eventSd, 0.0);
        CuAssertDblEquals(testCase, eventDuration, t->eventDuration, 0.0);
    }
    NanoporeReadAdjustmentParameters *params = nanopore_readAdjustmentParametersConstruct();
    StateMachine *sM = loadDescaledStateMachine3(npRead);
    nanopore_compute_scale_params(sM->EMISSION_MATCH_MATRIX, map, params, FALSE, TRUE);

    //CuAssertTrue(testCase, absPercentDiff(params->scale, npRead->templateParams.scale) < 50.0);
    //CuAssertTrue(testCase, absPercentDiff(params->shift, npRead->templateParams.shift) < 50.0);
    //CuAssertTrue(testCase, absPercentDiff(params->var, npRead->templateParams.var) < 50.0);

    st_uglyf("npRead scale: %f, estimated: %f diff %f\n", npRead->templateParams.scale, params->scale,
             absPercentDiff(params->scale, npRead->templateParams.scale));
    st_uglyf("npRead shift: %f, estimated: %f diff %f\n", npRead->templateParams.shift, params->shift,
             absPercentDiff(params->shift, npRead->templateParams.shift));
    st_uglyf("npRead var: %f, estimated: %f diff %f\n", npRead->templateParams.var, params->var,
             absPercentDiff(params->var, npRead->templateParams.var));
}

static void test_sm3_diagonalDPCalculations(CuTest *testCase) {
    // make some DNA sequences and fake nanopore read data
    char *sX = "ACGATAXGGACAT";
    double sY[21] = {
            58.743435, 0.887833, 0.0571, //ACGATA 0
            53.604965, 0.816836, 0.0571, //CGATAC 1
            58.432015, 0.735143, 0.0571, //GATACG 2
            63.684352, 0.795437, 0.0571, //ATACGG 3
            //63.520262, 0.757803, 0.0571, //TACGGA skip
            58.921430, 0.812959, 0.0571, //ACGGAC 4
            59.895882, 0.740952, 0.0571, //CGGACA 5
            61.684303, 0.722332, 0.0571, //GGACAT 6
    };

    // make variables for the (corrected) length of the sequences
    int64_t lX = sequence_correctSeqLength(strlen(sX), kmer);
    int64_t lY = 7;

    // make Sequence objects
    //Sequence *SsX = makeKmerSequence(sX);
    Sequence *SsX = sequence_constructKmerSequence(lX, sX, sequence_getKmer, sequence_sliceNucleotideSequence,
                                                   THREE_CYTOSINES, NB_CYTOSINE_OPTIONS, kmer);
    Sequence *SsY = sequence_construct(lY, sY, sequence_getEvent, event);

    // make stateMachine, forward and reverse DP matrices and banding stuff
    char *modelFile = stString_print("../../signalAlign/models/testModel_template.model");
    StateMachine *sM = getStateMachine3(modelFile);

    DpMatrix *dpMatrixForward = dpMatrix_construct(lX + lY, sM->stateNumber);
    DpMatrix *dpMatrixBackward = dpMatrix_construct(lX + lY, sM->stateNumber);
    stList *anchorPairs = stList_construct();
    Band *band = band_construct(anchorPairs, SsX->length, SsY->length, 2);
    BandIterator *bandIt = bandIterator_construct(band);

    // Initialize Matrices
    for (int64_t i = 0; i <= lX + lY; i++) {
        Diagonal d = bandIterator_getNext(bandIt);
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(dpMatrixBackward, d, SsX));
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(dpMatrixForward, d, SsX));
    }
    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(dpMatrixForward, 0), sM, sM->startStateProb);
    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(dpMatrixBackward, lX + lY), sM, sM->endStateProb);

    //Forward algorithm
    for (int64_t i = 1; i <= lX + lY; i++) {
        diagonalCalculationForward(sM, i, dpMatrixForward, SsX, SsY);
    }
    //Backward algorithm
    for (int64_t i = lX + lY; i > 0; i--) {
        diagonalCalculationBackward(sM, i, dpMatrixBackward, SsX, SsY);
    }

    //Calculate total probabilities
    double totalProbForward = LOG_ZERO;
    DpDiagonal *dpDiagonalF = dpMatrix_getDiagonal(dpMatrixForward, lX + lY);
    HDCell *cellF = dpDiagonal_getCell(dpDiagonalF, lX - lY);
    for (int64_t p = 0; p < cellF->numberOfPaths; p++) {
        Path *path = hdCell_getPath(cellF, p);
        totalProbForward = logAdd(totalProbForward, cell_dotProduct2(path->cells, sM, sM->endStateProb));
    }

    double totalProbBackward = LOG_ZERO;
    DpDiagonal *dpDiagonalB = dpMatrix_getDiagonal(dpMatrixBackward, 0);
    HDCell *cellB = dpDiagonal_getCell(dpDiagonalB, 0);
    for (int64_t p = 0; p < cellB->numberOfPaths; p++) {
        Path *path = hdCell_getPath(cellB, p);
        totalProbBackward = logAdd(totalProbBackward, cell_dotProduct2(path->cells, sM, sM->startStateProb));
    }
    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.001);

    // Test the posterior probabilities along the diagonals of the matrix.
    for (int64_t i = 0; i <= lX + lY; i++) {
        double totalDiagonalProb = diagonalCalculationTotalProbability(sM, i,
                                                                       dpMatrixForward,
                                                                       dpMatrixBackward,
                                                                       SsX, SsY);
        //Check the forward and back probabilities are about equal
        CuAssertDblEquals(testCase, totalProbForward, totalDiagonalProb, 0.01);
    }


    // Now do the posterior probabilities, get aligned pairs with posterior match probs above threshold
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    void *extraArgs[1] = { alignedPairs };
    for (int64_t i = 1; i <= lX + lY; i++) {
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        p->threshold = 0.2;
        diagonalCalculationPosteriorMatchProbs(sM, i, dpMatrixForward, dpMatrixBackward, SsX, SsY,
                                               totalProbForward, p, extraArgs);
        pairwiseAlignmentBandingParameters_destruct(p);
    }

    // Make a list of the correct anchor points
    stSortedSet *alignedPairsSet = stSortedSet_construct3((int (*)(const void *, const void *)) stIntTuple_cmpFn,
                                                          (void (*)(void *)) stIntTuple_destruct);

    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(0, 0));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(1, 1));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(2, 2));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(3, 3));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(4, 3));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(5, 4));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(6, 5));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(7, 6));

    // make sure alignedPairs is correct
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *pair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(pair, 1), y = stIntTuple_get(pair, 2);
        //char *pairKmer = (char *)stIntTuple_get(pair, 3);
        //st_uglyf("Pair %f %" PRIi64 " %" PRIi64 " kmer %s\n", (float) stIntTuple_get(pair, 0) / PAIR_ALIGNMENT_PROB_1,
        //         x, y, pairKmer);
        CuAssertTrue(testCase, stSortedSet_search(alignedPairsSet, stIntTuple_construct2(x, y)) != NULL);
    }
    CuAssertIntEquals(testCase, 14, (int) stList_length(alignedPairs));

    // clean up
    stateMachine_destruct(sM);
    sequence_destruct(SsX);
    sequence_destruct(SsY);
}

static void test_stateMachine3_getAlignedPairsWithBanding(CuTest *testCase) {
    // load the reference sequence and the nanopore read
    Sequence *refSeq = getTestReferenceSequence();
    NanoporeRead *npRead = loadTestNanoporeRead();

    // load stateMachine(s)
    StateMachine *sM = loadScaledStateMachine3(npRead);
    StateMachine *sMdescaled = loadDescaledStateMachine3(npRead);
    // parameters for pairwise alignment using defaults
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

    // get anchors using lastz
    stList *filteredRemappedAnchors = getRemappedAnchors(refSeq, npRead, p);

    // make the event sequence object
    Sequence *eventSequence = sequence_construct2(npRead->nbTemplateEvents, npRead->templateEvents, sequence_getEvent,
                                                  sequence_sliceEventSequence, event);

    // do alignment with scaled stateMachine
    stList *alignedPairs = getAlignedPairsUsingAnchors(sM, refSeq, eventSequence, filteredRemappedAnchors, p,
                                                       diagonalCalculationPosteriorMatchProbs,
                                                       0, 0);
    checkAlignedPairs(testCase, alignedPairs, refSeq->length, npRead->nbTemplateEvents);
    // for ch1_file1 template there should be this many aligned pairs with banding
    //st_uglyf("got %lld alignedPairs with anchors\n", stList_length(alignedPairs));
    CuAssertTrue(testCase, stList_length(alignedPairs) == 1069);

    stList *alignedPairs_descaled = getAlignedPairsUsingAnchors(sMdescaled,
                                                                refSeq, eventSequence, filteredRemappedAnchors, p,
                                                                diagonalCalculationPosteriorMatchProbs,
                                                                0, 0);

    // for ch1_file1 template there should be this many aligned pairs with banding
    CuAssertTrue(testCase, stList_length(alignedPairs_descaled) == 1069);

    // check against alignment without banding
    stList *alignedPairs_noband = getAlignedPairsWithoutBanding(sM, refSeq->elements, npRead->templateEvents,
                                                                refSeq->length, npRead->nbTemplateEvents, p,
                                                                sequence_getKmer, sequence_getEvent,
                                                                diagonalCalculationPosteriorMatchProbs,
                                                                0, 0);

    checkAlignedPairs(testCase, alignedPairs_noband, refSeq->length, npRead->nbTemplateEvents);
    //st_uglyf("got %lld alignedPairs with anchors\n", stList_length(alignedPairs_noband));
    CuAssertTrue(testCase, stList_length(alignedPairs_noband) == 1069);

    // clean
    pairwiseAlignmentBandingParameters_destruct(p);
    nanopore_nanoporeReadDestruct(npRead);
    sequence_destruct(refSeq);
    sequence_destruct(eventSequence);
    stList_destruct(alignedPairs);
    stList_destruct(alignedPairs_descaled);
    stList_destruct(alignedPairs_noband);
    stateMachine_destruct(sM);
    stateMachine_destruct(sMdescaled);
}

static void test_sm3Hdp_getAlignedPairsWithBanding(CuTest *testCase) {
    // load the reference sequence and the nanopore read
    Sequence *refSeq = getTestReferenceSequence();
    NanoporeRead *npRead = loadTestNanoporeRead();

    // this is a hack for this test so that I don't have to load the 200MB hdp file
    nanopore_descaleNanoporeRead(npRead);

    StateMachine *sM = loadStateMachineHdp(npRead);

    // parameters for pairwise alignment using defaults
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
    p->threshold = 0.1;

    // get anchors using lastz
    stList *filteredRemappedAnchors = getRemappedAnchors(refSeq, npRead, p);

    Sequence *eventSequence = sequence_construct2(npRead->nbTemplateEvents, npRead->templateEvents, sequence_getEvent,
                                                  sequence_sliceEventSequence, event);
    // do alignment of template events
    stList *alignedPairs = getAlignedPairsUsingAnchors(sM, refSeq, eventSequence, filteredRemappedAnchors, p,
                                                       diagonalCalculationPosteriorMatchProbs,
                                                       0, 0);
    checkAlignedPairs(testCase, alignedPairs, refSeq->length, npRead->nbTemplateEvents);

    // for ch1_file1 template there should be this many aligned pairs with banding
    //st_uglyf("got %lld alignedPairs with anchors\n", stList_length(alignedPairs));
    CuAssertTrue(testCase, stList_length(alignedPairs) == 1213);

    // clean
    pairwiseAlignmentBandingParameters_destruct(p);
    nanopore_nanoporeReadDestruct(npRead);
    sequence_destruct(refSeq);
    sequence_destruct(eventSequence);
    stList_destruct(alignedPairs);
    stateMachine_destruct(sM);
}

static void test_DegenerateNucleotides(CuTest *testCase) {
    Sequence *refSeq = getTestReferenceSequence();
    NanoporeRead *npRead = loadTestNanoporeRead();

    StateMachine *sM = loadDescaledStateMachine3(npRead);

    // parameters for pairwise alignment using defaults
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

    stList *filteredRemappedAnchors = getRemappedAnchors(refSeq, npRead, p);

    // make Sequences for reference and template events
    Sequence *eventSequence = sequence_construct2(npRead->nbTemplateEvents, npRead->templateEvents, sequence_getEvent,
                                                  sequence_sliceEventSequence, event);

    // do alignment of template events
    stList *alignedPairs = getAlignedPairsUsingAnchors(sM, refSeq, eventSequence, filteredRemappedAnchors, p,
                                                       diagonalCalculationPosteriorMatchProbs,
                                                       0, 0);

    checkAlignedPairs(testCase, alignedPairs, refSeq->length, npRead->nbTemplateEvents);
    CuAssertTrue(testCase, stList_length(alignedPairs) == 1069);
    //st_uglyf("got %lld alignedPairs with anchors\n", stList_length(alignedPairs));
    Sequence *degenerateSequence = replaceBasesInSequence(refSeq, "C", "E");
    stList *alignedPairs_methyl = getAlignedPairsUsingAnchors(sM, degenerateSequence,
                                                              eventSequence, filteredRemappedAnchors, p,
                                                              diagonalCalculationPosteriorMatchProbs,
                                                              0, 0);
    checkAlignedPairs(testCase, alignedPairs_methyl, degenerateSequence->length, npRead->nbTemplateEvents);
    CuAssertTrue(testCase, stList_length(alignedPairs_methyl) == 1069);
    sequence_destruct(degenerateSequence);
    degenerateSequence = replaceBasesInSequence(refSeq, "C", "O");
    stList *alignedPairs_hydroxy = getAlignedPairsUsingAnchors(sM, degenerateSequence,
                                                               eventSequence, filteredRemappedAnchors, p,
                                                               diagonalCalculationPosteriorMatchProbs,
                                                               0, 0);
    checkAlignedPairs(testCase, alignedPairs_hydroxy, degenerateSequence->length, npRead->nbTemplateEvents);
    CuAssertTrue(testCase, stList_length(alignedPairs_hydroxy) == 1069);


    sequence_destruct(degenerateSequence);
    degenerateSequence = replaceBasesInSequence(refSeq, "C", "X");
    stList *alignedPairs_degenerate = getAlignedPairsUsingAnchors(sM, degenerateSequence,
                                                                  eventSequence, filteredRemappedAnchors, p,
                                                                  diagonalCalculationPosteriorMatchProbs,
                                                                  0, 0);
    checkAlignedPairsWithOverlap(testCase, alignedPairs_degenerate,
                                 degenerateSequence->length, npRead->nbTemplateEvents);
    CuAssertTrue(testCase, stList_length(alignedPairs_degenerate) == 7256);
    // for ch1_file1 template there should be this many aligned pairs with banding
    //st_uglyf("got %lld alignedPairs on the degenerate sequence\n", stList_length(alignedPairs_degenerate));

    // clean
    pairwiseAlignmentBandingParameters_destruct(p);
    nanopore_nanoporeReadDestruct(npRead);
    sequence_destruct(refSeq);
    sequence_destruct(eventSequence);
    sequence_destruct(degenerateSequence);
    stList_destruct(alignedPairs);
    stList_destruct(alignedPairs_methyl);
    stList_destruct(alignedPairs_hydroxy);
    stList_destruct(alignedPairs_degenerate);
    stateMachine_destruct(sM);
}

static void test_cpHmmEmissionsAgainstStateMachine(CuTest *testCase, StateMachine *sM, ContinuousPairHmm *cpHmm) {
    for (int64_t i = 0; i < NUM_OF_KMERS; i++) {
        double E_mean = *(cpHmm->getEventModelEntry((Hmm *)cpHmm, i));
        double E_noise = *(cpHmm->getEventModelEntry((Hmm *)cpHmm, i) + 1);
        CuAssertDblEquals(testCase, sM->EMISSION_MATCH_MATRIX[1 + (i * MODEL_PARAMS)],
                          cpHmm->eventModel[i * NORMAL_DISTRIBUTION_PARAMS], 0.0);
        CuAssertDblEquals(testCase, sM->EMISSION_MATCH_MATRIX[1 + (i * MODEL_PARAMS + 1)],
                          cpHmm->eventModel[i * NORMAL_DISTRIBUTION_PARAMS + 1], 0.0);
        CuAssertDblEquals(testCase, sM->EMISSION_MATCH_MATRIX[1 + (i * MODEL_PARAMS)], E_mean, 0.0);
        CuAssertDblEquals(testCase, sM->EMISSION_MATCH_MATRIX[1 + (i * MODEL_PARAMS + 1)], E_noise, 0.0);
    }
}

static void test_makeAndCheckModels(CuTest *testCase) {
    // make a blank HMM
    Hmm *hmm = continuousPairHmm_construct(0.0, 0.0, 3, NUM_OF_KMERS, threeState, 0.0, 0.0, 0.0);
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *)hmm;

    // this is the lookup table with default values
    const char *templateLookupTableFile = "../../signalAlign/models/testModel_template.model";
    if (!stFile_exists(templateLookupTableFile)) {
        st_errAbort("Didn't find model file %s\n", templateLookupTableFile);
    }

    // load the table into the HMM emissions
    continuousPairHmm_loadModelFromFile(cpHmm, templateLookupTableFile);
    if (!cpHmm->hasModel) {
        st_errAbort("Problem loading match model\n");
    }

    // make a stateMachine based on the same table
    StateMachine *sM = getStateMachine3(templateLookupTableFile);

    // check that the emissions are correct
    test_cpHmmEmissionsAgainstStateMachine(testCase, sM, cpHmm);
    CuAssertTrue(testCase, sM->stateNumber == cpHmm->baseHmm.stateNumber);

    // this is bogus
}

static void test_continuousPairHmm(CuTest *testCase) {
    // construct the empty hmm
    Hmm *hmm = continuousPairHmm_construct(0.0, 0.0, 3, NUM_OF_KMERS, threeState, 0.0, 0.0, 0.0);
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *)hmm;

    // test that it is empty
    CuAssertTrue(testCase, cpHmm->hasModel == FALSE);

    // Add some transition expectations
    int64_t nStates = cpHmm->baseHmm.stateNumber;
    CuAssertTrue(testCase, nStates == 3);
    for (int64_t from = 0; from < nStates; from++) {
        for (int64_t to = 0; to < nStates; to++) {
            double dummy = from * nStates + to;
            cpHmm->baseHmm.addToTransitionExpectationFcn((Hmm *)cpHmm, from, to, dummy);
        }
    }

    // Add event model
    const char *model = "../../signalAlign/models/testModel_template.model";
    if (!stFile_exists(model)) {
        st_errAbort("Didn't find model file %s\n", model);
    }
    continuousPairHmm_loadModelFromFile(cpHmm, model);
    CuAssertTrue(testCase, cpHmm->hasModel);

    // test that it loaded correctly
    StateMachine *sM = getStateMachine3(model);
    test_cpHmmEmissionsAgainstStateMachine(testCase, sM, cpHmm);

    // dump to file
    char *tempFile = stString_print("./temp%" PRIi64 ".hmm", st_randomInt(0, INT64_MAX));
    CuAssertTrue(testCase, !stFile_exists(tempFile));  // Quick check that we don't write over anything.
    FILE *fH = fopen(tempFile, "w");
    continuousPairHmm_dump((Hmm *)cpHmm, fH);
    continuousPairHmm_destruct((Hmm *)cpHmm);

    // load
    hmm = continuousPairHmm_loadFromFile(tempFile, 0.0, 0.0);
    stFile_rmrf(tempFile);
    cpHmm = (ContinuousPairHmm *)hmm;
    CuAssertTrue(testCase, cpHmm->hasModel);

    // (re)check the transitions
    for (int64_t from = 0; from < nStates; from++) {
        for (int64_t to = 0; to < nStates; to++) {
            double retrievedProb = cpHmm->baseHmm.getTransitionsExpFcn((Hmm *)cpHmm, from, to);
            double correctProb = from * nStates + to;
            CuAssertTrue(testCase, retrievedProb == correctProb);
        }
    }

    // recheck that the model loaded correctly
    test_cpHmmEmissionsAgainstStateMachine(testCase, sM, cpHmm);

    // change the expectations
    for (int64_t i = 0; i < cpHmm->baseHmm.symbolSetSize; i++) {
        double mean = *(cpHmm->getEventModelEntry((Hmm *)cpHmm, i));
        cpHmm->addToEmissionExpectationFcn((Hmm *)cpHmm, i, mean, 1);
    }

    // check it
    for (int64_t i = 0; i < cpHmm->baseHmm.symbolSetSize; i++) {
        double x = *(cpHmm->getEmissionExpFcn((Hmm *)cpHmm, i));
        double y = *(cpHmm->getEmissionExpFcn((Hmm *)cpHmm, i) + 1);
        double E_mean = *(cpHmm->getEventModelEntry((Hmm *)cpHmm, i));
        double E_sd = 0;
        CuAssertDblEquals(testCase, E_mean, x, 0.0);
        CuAssertDblEquals(testCase, E_sd, y, 0.0);
    }

    // check posteriors
    for (int64_t i = 0; i < cpHmm->baseHmm.symbolSetSize; i++) {
        CuAssertDblEquals(testCase, 1, *(cpHmm->getPosteriorExpFcn((Hmm *)cpHmm, i)), 0.0); // todo fails
    }

    // add new expectations at 2 * E(mean)
    for (int64_t i = 0; i < cpHmm->baseHmm.symbolSetSize; i++) {
        double mean = 2 * *(cpHmm->getEventModelEntry((Hmm *)cpHmm, i));
        cpHmm->addToEmissionExpectationFcn((Hmm *)cpHmm, i, mean, 1);
    }

    // check posteriors
    for (int64_t i = 0; i < cpHmm->baseHmm.symbolSetSize; i++) {
        CuAssertDblEquals(testCase, 2, *(cpHmm->getPosteriorExpFcn((Hmm *)cpHmm, i)), 0.0); // todo fails
    }

    // check that they got added
    for (int64_t i = 0; i < cpHmm->baseHmm.symbolSetSize; i++) {
        double x = *(cpHmm->getEmissionExpFcn((Hmm *)cpHmm, i));
        double y = *(cpHmm->getEmissionExpFcn((Hmm *)cpHmm, i) + 1);
        double E_mean = *(cpHmm->getEventModelEntry((Hmm *)cpHmm, i)) * 3;
        double twoX = 2 * *(cpHmm->getEventModelEntry((Hmm *)cpHmm, i));
        double onePointFiveX = 1.5 * *(cpHmm->getEventModelEntry((Hmm *)cpHmm, i));
        double E_sd = (twoX - onePointFiveX) * (twoX - onePointFiveX);
        CuAssertDblEquals(testCase, E_mean, x, 0.0);
        CuAssertDblEquals(testCase, E_sd, y, 0.0); // TODO fails
    }

    // check posteriors
    for (int64_t i = 0; i < cpHmm->baseHmm.symbolSetSize; i++) {
        CuAssertDblEquals(testCase, 2, *(cpHmm->getPosteriorExpFcn((Hmm *)cpHmm, i)), 0.0); // todo fails
    }

    // normalize it
    continuousPairHmm_normalize((Hmm *)cpHmm);

    // check that the posteriors didn't change
    for (int64_t i = 0; i < cpHmm->baseHmm.symbolSetSize; i++) {
        CuAssertDblEquals(testCase, 2, *(cpHmm->getPosteriorExpFcn((Hmm *)cpHmm, i)), 0.0); // todo fails
    }

    // Recheck transitions
    for (int64_t from = 0; from < nStates; from++) {
        for (int64_t to = 0; to < nStates; to++) {
            double z = from * nStates * nStates + (nStates * (nStates - 1)) / 2;
            double retrievedNormedProb = cpHmm->baseHmm.getTransitionsExpFcn((Hmm *)cpHmm, from, to);
            CuAssertDblEquals(testCase, (from * nStates + to) / z, retrievedNormedProb, 0.0);
        }
    }

    // check the emissions
    for (int64_t i = 0; i < cpHmm->baseHmm.symbolSetSize; i++) {
        double origMean = sM->EMISSION_MATCH_MATRIX[1 + (i * MODEL_PARAMS)];
        double newMean = *(cpHmm->getEventModelEntry((Hmm *)cpHmm, i));

        double twoX = 2 * sM->EMISSION_MATCH_MATRIX[1 + (i * MODEL_PARAMS)];
        double onePointFiveX = 1.5 * sM->EMISSION_MATCH_MATRIX[1 + (i * MODEL_PARAMS)];
        double sq = (twoX - onePointFiveX) * (twoX - onePointFiveX);
        double E_sd = sqrt(sq / 2);

        double sd = *(cpHmm->getEventModelEntry((Hmm *)cpHmm, i) + 1);
        CuAssertDblEquals(testCase, 1.5 * origMean, newMean, 0.0); // todo fails
        CuAssertDblEquals(testCase, E_sd, sd, 0.0); // todo fails

    }
    continuousPairHmm_destruct((Hmm *)cpHmm);
}

static void test_hdpHmmWithoutAssignments(CuTest *testCase) {
    // make the hmm object
    Hmm *hmm = hdpHmm_constructEmpty(0.0, 3, threeStateHdp, 0.02);
    HdpHmm *hdpHmm = (HdpHmm *) hmm;

    // Add some transition expectations
    int64_t nStates = hdpHmm->baseHmm.stateNumber;

    for (int64_t from = 0; from < nStates; from++) {
        for (int64_t to = 0; to < nStates; to++) {
            double dummy = from * nStates + to;
            hdpHmm->baseHmm.addToTransitionExpectationFcn((Hmm *) hdpHmm, from, to, dummy);
        }
    }

    // make a simple nucleotide and event sequence
    char *sequence = "ACGTCATACATGACTATA";
    double fakeMeans[3] = { 65.0, 64.0, 63.0 };

    for (int64_t a = 0; a < 3; a++) {
        hdpHmm->addToAssignments(hdpHmm, sequence + (a * KMER_LENGTH), fakeMeans + a);
    }

    CuAssertTrue(testCase, hdpHmm->numberOfAssignments == 3);

    // dump the HMM to a file
    char *tempFile = stString_print("./temp%" PRIi64 ".hmm", st_randomInt(0, INT64_MAX));
    CuAssertTrue(testCase, !stFile_exists(tempFile)); //Quick check that we don't write over anything.
    FILE *fH = fopen(tempFile, "w");
    hdpHmm_writeToFile((Hmm *) hdpHmm, fH);
    fclose(fH);
    hdpHmm_destruct((Hmm *) hdpHmm);

    //Load from a file
    hmm = hdpHmm_loadFromFile(tempFile, NULL);
    hdpHmm = (HdpHmm *)hmm;
    stFile_rmrf(tempFile);

    CuAssertTrue(testCase, (3 == hdpHmm->numberOfAssignments));  // recheck number of assignments

    // Check the transition expectations
    for (int64_t from = 0; from < nStates; from++) {
        for (int64_t to = 0; to < nStates; to++) {
            double retrievedProb = hdpHmm->baseHmm.getTransitionsExpFcn((Hmm *) hdpHmm, from, to);
            double correctProb = from * nStates + to;
            CuAssertTrue(testCase, retrievedProb == correctProb);
        }
    }

    // normalize
    hmmDiscrete_normalizeTransitions((Hmm *) hdpHmm);

    // Recheck transitions
    for (int64_t from = 0; from < nStates; from++) {
        for (int64_t to = 0; to < nStates; to++) {
            double z = from * nStates * nStates + (nStates * (nStates - 1)) / 2;
            double retrievedNormedProb = hdpHmm->baseHmm.getTransitionsExpFcn((Hmm *) hdpHmm, from, to);
            CuAssertDblEquals(testCase, (from * nStates + to) / z, retrievedNormedProb, 0.0);
        }
    }
    hdpHmm_destruct((Hmm *) hdpHmm);
}

static void test_continuousPairHmm_em(CuTest *testCase) {
    Sequence *refSeq = getTestReferenceSequence();
    NanoporeRead *npRead = loadTestNanoporeRead();
    StateMachine *sM = loadDescaledStateMachine3(npRead);
    double pLikelihood = -INFINITY;

    // parameters for pairwise alignment using defaults
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
    stList *filteredRemappedAnchors = getRemappedAnchors(refSeq, npRead, p);

    Sequence *eventSequence = sequence_construct2(npRead->nbTemplateEvents,
                                                  npRead->templateEvents, sequence_getEvent,
                                                  sequence_sliceEventSequence, event);

    for (int64_t iter = 0; iter < 10; iter++) {
        Hmm *hmm = continuousPairHmm_construct(0.001, 0.001, 3, NUM_OF_KMERS, threeState,
                                               npRead->templateParams.scale,
                                               npRead->templateParams.shift,
                                               npRead->templateParams.var);

        //print_stateMachine_transitions(sM);
        ContinuousPairHmm *cpHmm = (ContinuousPairHmm *)hmm;

        implantModelFromStateMachineIntoHmm(sM, cpHmm);

        // get expectations
        getExpectationsUsingAnchors(sM, (Hmm *)cpHmm, refSeq, eventSequence, filteredRemappedAnchors,
                                    p, diagonalCalculation_Expectations, 0, 0);

        // normalize
        continuousPairHmm_normalize((Hmm *)cpHmm);

        //st_uglyf("->->-> Got expected likelihood %f for iteration %" PRIi64 "\n", cpHmm->baseHmm.likelihood, iter);
        // M step
        continuousPairHmm_loadTransitionsIntoStateMachine(sM, (Hmm *)cpHmm);
        continuousPairHmm_loadEmissionsIntoStateMachine(sM, (Hmm *)cpHmm);

        // Tests
        if (iter > 1) {
            assert(pLikelihood <= cpHmm->baseHmm.likelihood * 0.95);
            CuAssertTrue(testCase, pLikelihood <= cpHmm->baseHmm.likelihood * 0.85);
        }

        // update
        pLikelihood = cpHmm->baseHmm.likelihood;

        // per iteration clean up
        continuousPairHmm_destruct((Hmm *)cpHmm);
    }
    sequence_destruct(refSeq);
    sequence_destruct(eventSequence);
    stList_destruct(filteredRemappedAnchors);
    nanopore_nanoporeReadDestruct(npRead);
    pairwiseAlignmentBandingParameters_destruct(p);
    stateMachine_destruct(sM);
}

static void test_hdpHmm_emTransitions(CuTest *testCase) {
    Sequence *refSeq = getTestReferenceSequence();
    NanoporeRead *npRead = loadTestNanoporeRead();
    nanopore_descaleNanoporeRead(npRead);  // hack, remember
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

    double pLikelihood = -INFINITY;
    StateMachine *sM = loadStateMachineHdp(npRead);

    stList *filteredRemappedAnchors = getRemappedAnchors(refSeq, npRead, p);
    Sequence *templateSeq = sequence_construct2(npRead->nbTemplateEvents, npRead->templateEvents, sequence_getEvent,
                                                sequence_sliceEventSequence, event);

    for (int64_t iter = 0; iter < 10; iter++) {
        Hmm *hmmExpectations = hdpHmm_constructEmpty(0.0001, 3, threeStateHdp, p->threshold);

        getExpectationsUsingAnchors(sM, hmmExpectations, refSeq, templateSeq, filteredRemappedAnchors,
                                    p, diagonalCalculation_Expectations, 0, 0);

        // norm
        hmmDiscrete_normalizeTransitions(hmmExpectations);

        // for debugging
        //st_uglyf("->->-> Got expected likelihood %f for iteration %" PRIi64 "\n", hmmExpectations->likelihood, iter);
        //print_stateMachine_transitions(sM);

        // M step
        continuousPairHmm_loadTransitionsIntoStateMachine(sM, hmmExpectations);

        // Tests
        assert(pLikelihood <= hmmExpectations->likelihood * 0.95);
        if (iter > 1) {
            CuAssertTrue(testCase, pLikelihood <= hmmExpectations->likelihood * 0.95);
        }
        // update
        pLikelihood = hmmExpectations->likelihood;

        // per iteration clean up
        //continuousPairHmm_destruct(hmmExpectations);
        hdpHmm_destruct(hmmExpectations);
    }

    sequence_destruct(refSeq);
    sequence_destruct(templateSeq);
    stList_destruct(filteredRemappedAnchors);
    nanopore_nanoporeReadDestruct(npRead);
    pairwiseAlignmentBandingParameters_destruct(p);
    stateMachine_destruct(sM);
}

CuSuite *stateMachineAlignmentTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_nanoporeScaleParams);

    SUITE_ADD_TEST(suite, test_loadPoreModel);
    SUITE_ADD_TEST(suite, test_sm3_diagonalDPCalculations);
    SUITE_ADD_TEST(suite, test_stateMachine3_getAlignedPairsWithBanding);
    SUITE_ADD_TEST(suite, test_sm3Hdp_getAlignedPairsWithBanding);
    SUITE_ADD_TEST(suite, test_DegenerateNucleotides);
    SUITE_ADD_TEST(suite, test_makeAndCheckModels);
    SUITE_ADD_TEST(suite, test_hdpHmmWithoutAssignments);
    SUITE_ADD_TEST(suite, test_continuousPairHmm);
    SUITE_ADD_TEST(suite, test_continuousPairHmm_em);
    SUITE_ADD_TEST(suite, test_hdpHmm_emTransitions);

    return suite;
}
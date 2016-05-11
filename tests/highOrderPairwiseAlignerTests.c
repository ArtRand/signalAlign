// Tests for high order HMM

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include "randomSequences.h"
#include "stateMachine.h"
#include "CuTest.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
#include "multipleAligner.h"
#include "emissionMatrix.h"
#include "discreteHmm.h"
#include "nanopore.h"
#include "continuousHmm.h"

Sequence *makeTestKmerSequence() {
    char *s = "ATGXAXA"; // has 2 6mers
    int64_t lX = sequence_correctSeqLength(strlen(s), kmer);
    Sequence *seq = sequence_construct(lX, s, sequence_getKmer3, kmer);
    seq->cytosines = "CEO";
    seq->nbCytosineOptions = 3;
    return seq;
}

Sequence *makeKmerSequence(char *nucleotides) {
    int64_t lX = sequence_correctSeqLength(strlen(nucleotides), kmer);
    Sequence *seq = sequence_constructReferenceSequence(lX, nucleotides,
                                                        sequence_getKmer3, sequence_sliceNucleotideSequence2,
                                                        "CEO", 3, kmer);
    return seq;
}

static void breakpoint(char *msg) {
    st_errAbort("Breakpoint %s\n", msg);
}

static void test_findPotentialMethylation(CuTest *testCase) {
    stList *noMethyls = path_findPotentialMethylation("ATGCAC");
    CuAssertTrue(testCase, stList_length(noMethyls) == 0);
    stList *methyls = path_findPotentialMethylation("ATGXAX");
    CuAssertTrue(testCase, stList_length(methyls) == 2);
    int64_t methyl_1 = *(int64_t *)stList_get(methyls, 0);
    int64_t methyl_2 = *(int64_t *)stList_get(methyls, 1);
    CuAssertTrue(testCase, methyl_1 == 3);
    CuAssertTrue(testCase, methyl_2 == 5);
}

static void test_pathLegalTransitions(CuTest *testCase) {
    Path *path0 = path_construct(NULL, 3);
    Path *path1 = path_construct("AAETTT", 3);
    Path *path2 = path_construct("AETTTC", 3);
    Path *path3 = path_construct("AETTTC", 2);
    Path *path4 = path_construct("ACTTTC", 3);
    CuAssertTrue(testCase, path_checkLegal(path0, path1) == TRUE);
    CuAssertTrue(testCase, path_checkLegal(path1, path2) == TRUE);
    CuAssertTrue(testCase, path_checkLegal(path1, path3) == FALSE);
    CuAssertTrue(testCase, path_checkLegal(path1, path4) == FALSE);
    path_destruct(path0);
    path_destruct(path1);
    path_destruct(path2);
    path_destruct(path3);
    path_destruct(path4);
}

static void test_getKmer4(CuTest *testCase) {
    Sequence *sX = makeKmerSequence("ATGCATAGC");
                                   //ATGCAT
                                   // TGCATA
                                   //  GCATAG
                                   //   CATAGC
    void *checkMinusOne = sequence_getKmer4(sX, -1);
    CuAssertTrue(testCase, checkMinusOne == NULL);
    void *checkGreaterThanLength = sequence_getKmer4(sX, (sX->length + 1));
    CuAssertTrue(testCase, checkGreaterThanLength == NULL);
}

static void test_methylPermutations(CuTest *testCase) {
    CuAssertTrue(testCase, path_getMehtylPermutations(0, 3, "CEO") == NULL);
    for (int64_t i = 1; i < 6; i++) {
        stList *patterns = path_getMehtylPermutations(i, 3, "CEO");
        int64_t correctLength = intPow(3, i);
        CuAssertTrue(testCase, stList_length(patterns) == correctLength);
        stList_destruct(patterns);
    }
}

static void test_substitutedKmers(CuTest *testCase) {
    char *ambigKmer = "ATGXAX";
    stList *positions = path_findPotentialMethylation(ambigKmer);
    char *pattern = "CE";
    char *kmer = hdCell_getSubstitutedKmer(positions, 2, pattern, ambigKmer);
    CuAssertTrue(testCase, stList_length(positions) == 2);
    CuAssertStrEquals(testCase, kmer, "ATGCAE");
    stList_destruct(positions);
    free(kmer);
}

static void test_hdCellConstruct(CuTest *testCase) {
    char *ambigKmer = "ATGXAXAAAAAA";
    int64_t nbCytosines = 3;
    char *cytosines = "CEO";
    HDCell *cell = hdCell_construct(ambigKmer, 3, nbCytosines, cytosines);
    Path *path2 = hdCell_getPath(cell, 8);
    Path *path = hdCell_getPath(cell, 0);
    CuAssertTrue(testCase, hdCell_getPath(cell, 9) == NULL);
    CuAssertIntEquals(testCase, 9, (int )cell->numberOfPaths);
    CuAssertStrEquals(testCase, path->kmer, "ATGCAC");
    CuAssertStrEquals(testCase, path2->kmer, "ATGOAO");
    hdCell_destruct(cell);
}

static void test_hdCellConstructWorstCase(CuTest *testCase) {
    char *ambigKmer = "XXXXXX";
    int64_t nbCytosines = 3;
    char *cytosines = "CEO";
    HDCell *cell = hdCell_construct(ambigKmer, 3, nbCytosines, cytosines);
    Path *path = hdCell_getPath(cell, 0);
    Path *path2 = hdCell_getPath(cell, 728);
    CuAssertIntEquals(testCase, 729, (int )cell->numberOfPaths);
    CuAssertStrEquals(testCase, path->kmer, "CCCCCC");
    CuAssertStrEquals(testCase, path2->kmer, "OOOOOO");
    hdCell_destruct(cell);
}

static void test_dpDiagonal(CuTest *testCase) {
    // load model and make stateMachine
    char *testModel = stString_print("../../signalAlign/models/testModel_template.model");
    char *canonicalAlphabet = "ACGT\0";

    NanoporeHDP *nHdp = flat_hdp_model(canonicalAlphabet, 4, KMER_LENGTH, 1.0, 1.0, 40.0, 100.0, 100, testModel);
    StateMachine *sM = getHdpStateMachine3(nHdp, testModel);
    Diagonal diagonal = diagonal_construct(3, -1, 1); // makes a diagonal with 2 cells

    Sequence *seq = makeTestKmerSequence();  // ATGXAXA

    DpDiagonal *dpDiagonal = dpDiagonal_construct(diagonal, sM->stateNumber, seq);

    //Get cell
    HDCell *c1 = dpDiagonal_getCell(dpDiagonal, -1);
    CuAssertTrue(testCase, c1 != NULL);

    HDCell *c2 = dpDiagonal_getCell(dpDiagonal, 1); // gets cell 1
    CuAssertTrue(testCase, c2 != NULL);

    CuAssertTrue(testCase, dpDiagonal_getCell(dpDiagonal, 3) == NULL);
    CuAssertTrue(testCase, dpDiagonal_getCell(dpDiagonal, -3) == NULL);

    dpDiagonal_initialiseValues(dpDiagonal, sM, sM->endStateProb); //Test initialise values
    double totalProb = LOG_ZERO;

    for (int64_t p = 0; p < c1->numberOfPaths; p++) {
        Path *path1 = hdCell_getPath(c1, p);
        Path *path2 = hdCell_getPath(c2, p);
        for (int64_t s = 0; s < path1->stateNumber; s++) {
            CuAssertDblEquals(testCase, path1->cells[s], sM->endStateProb(sM, s), 0.0);
            CuAssertDblEquals(testCase, path2->cells[s], sM->endStateProb(sM, s), 0.0);
            //totalProb = logAdd(totalProb, (1/c1->numberOfPaths) * path1->cells[s]);
            //totalProb = logAdd(totalProb, (1/c1->numberOfPaths) * path2->cells[s]);
            totalProb = logAdd(totalProb, 2 * path1->cells[s]);
            totalProb = logAdd(totalProb, 2 * path2->cells[s]);
        }
    }

    DpDiagonal *dpDiagonal2 = dpDiagonal_clone(dpDiagonal);
    CuAssertTrue(testCase, dpDiagonal_equals(dpDiagonal, dpDiagonal2));
    //Check it runs
    CuAssertDblEquals(testCase, totalProb, dpDiagonal_dotProduct(dpDiagonal, dpDiagonal2), 0.001);

    dpDiagonal_destruct(dpDiagonal);
    dpDiagonal_destruct(dpDiagonal2);

}

static void test_dpMatrix(CuTest *testCase) {
    int64_t lX = 3, lY = 2;
    char *s = "ATGXAATT";
    //         ATGCAA   0
    //          TGCAAT  1
    //           GCAATT 2
    Sequence *sX = makeKmerSequence(s);

    DpMatrix *dpMatrix = dpMatrix_construct(lX + lY, 5);

    // check initialization
    CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), 0);

    // make sure there aren't any fantom diagonals
    for (int64_t i = -1; i <= lX + lY + 10; i++) {
        CuAssertTrue(testCase, dpMatrix_getDiagonal(dpMatrix, i) == NULL);
    }

    // make some diagonals in the dpMatrix, and check them, then make sure that
    // the number of active diagonals is correct.
    for (int64_t i = 0; i <= lX + lY; i++) {
        DpDiagonal *dpDiagonal = dpMatrix_createDiagonal(dpMatrix, diagonal_construct(i, -i, i), sX);
        CuAssertTrue(testCase, dpDiagonal == dpMatrix_getDiagonal(dpMatrix, i));
        CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), i + 1);
    }

    // test for destroying diagonals
    for (int64_t i = lX + lY; i >= 0; i--) {
        dpMatrix_deleteDiagonal(dpMatrix, i);
        CuAssertTrue(testCase, dpMatrix_getDiagonal(dpMatrix, i) == NULL);
        CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), i);
    }

    // double check that they are gone
    CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), 0);

    dpMatrix_destruct(dpMatrix);
}

static void test_sm3_diagonalDPCalculations(CuTest *testCase) {
    // make some DNA sequences and fake nanopore read data
    char *sX = "ACGATACGGACAT";
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
    int64_t lX = sequence_correctSeqLength(strlen(sX), event);
    int64_t lY = 7;

    // make Sequence objects
    //Sequence *SsX = sequence_construct(lX, sX, sequence_getKmer3, kmer);
    Sequence *SsX = makeKmerSequence(sX);
    Sequence *SsY = sequence_construct(lY, sY, sequence_getEvent, event);

    // make stateMachine, forward and reverse DP matrices and banding stuff
    char *modelFile = stString_print("../../signalAlign/models/testModel_template.model");
    StateMachine *sM = getStrawManStateMachine3(modelFile);

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
    CuAssertIntEquals(testCase, 8, (int) stList_length(alignedPairs));

    // clean up
    stateMachine_destruct(sM);
    sequence_sequenceDestroy(SsX);
    sequence_sequenceDestroy(SsY);
}

static void print_stateMachine_transitions(StateMachine *sM) {
    StateMachine3_HDP *sMhdp = (StateMachine3_HDP *)sM;
    st_uglyf("Match continue %f\n", exp(sMhdp->TRANSITION_MATCH_CONTINUE));
    st_uglyf("Match to gap x %f\n", exp(sMhdp->TRANSITION_GAP_OPEN_X));
    st_uglyf("Match to gap y %f\n", exp(sMhdp->TRANSITION_GAP_OPEN_Y));

    st_uglyf("gap X continue %f\n", exp(sMhdp->TRANSITION_GAP_EXTEND_X));
    st_uglyf("gap X to match %f\n", exp(sMhdp->TRANSITION_MATCH_FROM_GAP_X));
    st_uglyf("gap X from gap Y %f\n", exp(sMhdp->TRANSITION_GAP_SWITCH_TO_X));

    st_uglyf("gap Y continue %f\n", exp(sMhdp->TRANSITION_GAP_EXTEND_Y));
    st_uglyf("gap Y to match %f\n", exp(sMhdp->TRANSITION_MATCH_FROM_GAP_Y));
    st_uglyf("gap Y from gap X %f\n", exp(sMhdp->TRANSITION_GAP_SWITCH_TO_Y));
}
static void test_sm3Hdp_diagonalDPCalculations(CuTest *testCase) {
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
    int64_t lX = sequence_correctSeqLength(strlen(sX), event);
    int64_t lY = 7;

    // make Sequence objects
    Sequence *SsX = sequence_construct(lX, sX, sequence_getKmer3, kmer);
    Sequence *SsY = sequence_construct(lY, sY, sequence_getEvent, event);
    char *modelFile = stString_print("../../signalAlign/models/testModel_template.model");
    NanoporeHDP *nHdp = deserialize_nhdp("../../signalAlign/models/templateSingleLevelFixed.nhdp");
    StateMachine *sM = getHdpStateMachine3(nHdp, modelFile);

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
    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.01);

    ;
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
    sequence_sequenceDestroy(SsX);
    sequence_sequenceDestroy(SsY);
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

static void checkAlignedPairsForEchelon(CuTest *testCase, stList *blastPairs, int64_t lX, int64_t lY) {
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

static void test_strawMan_getAlignedPairsWithBanding(CuTest *testCase) {
    // load the reference sequence and the nanopore read
    char *ZymoReference = stString_print("../../signalAlign/tests/test_npReads/ZymoRef.txt");
    FILE *fH = fopen(ZymoReference, "r");
    char *ZymoReferenceSeq = stFile_getLineFromFile(fH);
    char *npReadFile = stString_print("../../signalAlign/tests/test_npReads/ZymoC_ch_1_file1.npRead");
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);

    // get sequence lengths
    int64_t lX = sequence_correctSeqLength(strlen(ZymoReferenceSeq), event);
    int64_t lY = npRead->nbTemplateEvents;

    // load stateMachine from model file
    char *templateModelFile = stString_print("../../signalAlign/models/testModel_template.model");
    StateMachine *sMt = getStrawManStateMachine3(templateModelFile);

    // scale model
    emissions_signal_scaleModel(sMt, npRead->templateParams.scale, npRead->templateParams.shift,
                                npRead->templateParams.var, npRead->templateParams.scale_sd,
                                npRead->templateParams.var_sd); // clunky

    // parameters for pairwise alignment using defaults
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

    // get anchors using lastz
    stList *anchorPairs = getBlastPairsForPairwiseAlignmentParameters(ZymoReferenceSeq, npRead->twoDread, p);

    // remap and filter
    stList *remappedAnchors = nanopore_remapAnchorPairs(anchorPairs, npRead->templateEventMap);
    stList *filteredRemappedAnchors = filterToRemoveOverlap(remappedAnchors);

    // make Sequences for reference and template events
    Sequence *refSeq = sequence_constructReferenceSequence(lX, ZymoReferenceSeq,
                                                           sequence_getKmer3, sequence_sliceNucleotideSequence2,
                                                           "CEO", 3, kmer);

    Sequence *templateSeq = sequence_construct2(lY, npRead->templateEvents, sequence_getEvent,
                                                sequence_sliceEventSequence2, event);

    // do alignment of template events
    stList *alignedPairs = getAlignedPairsUsingAnchors(sMt, refSeq, templateSeq, filteredRemappedAnchors, p,
                                                       diagonalCalculationPosteriorMatchProbs,
                                                       0, 0);
    checkAlignedPairs(testCase, alignedPairs, lX, lY);


    // for ch1_file1 template there should be this many aligned pairs with banding
    //st_uglyf("got %lld alignedPairs with anchors\n", stList_length(alignedPairs));
    // used to be 1001 before I added prior skip prob
    CuAssertTrue(testCase, stList_length(alignedPairs) == 1069);

    // check against alignment without banding

    stList *alignedPairs2 = getAlignedPairsWithoutBanding(sMt, ZymoReferenceSeq, npRead->templateEvents, lX,
                                                          lY, p,
                                                          sequence_getKmer, sequence_getEvent,
                                                          diagonalCalculationPosteriorMatchProbs,
                                                          0, 0);

    checkAlignedPairs(testCase, alignedPairs2, lX, lY);
    //st_uglyf("got %lld alignedPairs with anchors\n", stList_length(alignedPairs2));
    CuAssertTrue(testCase, stList_length(alignedPairs2) == 1069);

    // clean
    pairwiseAlignmentBandingParameters_destruct(p);
    nanopore_nanoporeReadDestruct(npRead);
    sequence_sequenceDestroy(refSeq);
    sequence_sequenceDestroy(templateSeq);
    stList_destruct(alignedPairs);
    //stList_destruct(alignedPairs2);
    stateMachine_destruct(sMt);
}

static void test_strawMan_getDescaledAlignedPairsWithBanding(CuTest *testCase) {
    // load the reference sequence and the nanopore read
    char *ZymoReference = stString_print("../../signalAlign/tests/test_npReads/ZymoRef.txt");
    FILE *fH = fopen(ZymoReference, "r");
    char *ZymoReferenceSeq = stFile_getLineFromFile(fH);
    char *npReadFile = stString_print("../../signalAlign/tests/test_npReads/ZymoC_ch_1_file1.npRead");
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);

    // get sequence lengths
    int64_t lX = sequence_correctSeqLength(strlen(ZymoReferenceSeq), event);
    int64_t lY = npRead->nbTemplateEvents;

    // load stateMachine from model file
    char *templateModelFile = stString_print("../../signalAlign/models/testModel_template.model");
    StateMachine *sMt = getSM3_descaled(templateModelFile, npRead->templateParams);

    // parameters for pairwise alignment using defaults
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

    // get anchors using lastz
    stList *anchorPairs = getBlastPairsForPairwiseAlignmentParameters(ZymoReferenceSeq, npRead->twoDread, p);

    // remap and filter
    stList *remappedAnchors = nanopore_remapAnchorPairs(anchorPairs, npRead->templateEventMap);
    stList *filteredRemappedAnchors = filterToRemoveOverlap(remappedAnchors);

    // make Sequences for reference and template events
    Sequence *refSeq = sequence_constructReferenceSequence(lX, ZymoReferenceSeq,
                                                           sequence_getKmer3, sequence_sliceNucleotideSequence2,
                                                           THREE_CYTOSINES, NB_CYTOSINE_OPTIONS, kmer);
    Sequence *templateSeq = sequence_construct2(lY, npRead->templateEvents, sequence_getEvent,
                                                sequence_sliceEventSequence2, event);

    // do alignment of template events
    stList *alignedPairs = getAlignedPairsUsingAnchors(sMt, refSeq, templateSeq, filteredRemappedAnchors, p,
                                                       diagonalCalculationPosteriorMatchProbs,
                                                       0, 0);
    checkAlignedPairs(testCase, alignedPairs, lX, lY);


    // for ch1_file1 template there should be this many aligned pairs with banding
    //st_uglyf("got %lld alignedPairs with anchors\n", stList_length(alignedPairs));
    // used to be 1001 before I added prior skip prob
    CuAssertTrue(testCase, stList_length(alignedPairs) == 1069);

    // check against alignment without banding

    stList *alignedPairs2 = getAlignedPairsWithoutBanding(sMt, ZymoReferenceSeq, npRead->templateEvents, lX,
                                                          lY, p,
                                                          sequence_getKmer, sequence_getEvent,
                                                          diagonalCalculationPosteriorMatchProbs,
                                                          0, 0);

    checkAlignedPairs(testCase, alignedPairs2, lX, lY);
    //st_uglyf("got %lld alignedPairs without anchors\n", stList_length(alignedPairs2));
    CuAssertTrue(testCase, stList_length(alignedPairs2) == 1069);

    // clean
    pairwiseAlignmentBandingParameters_destruct(p);
    nanopore_nanoporeReadDestruct(npRead);
    sequence_sequenceDestroy(refSeq);
    sequence_sequenceDestroy(templateSeq);
    stList_destruct(alignedPairs);
    stList_destruct(alignedPairs2);
    stateMachine_destruct(sMt);
}

static void test_sm3Hdp_getAlignedPairsWithBanding(CuTest *testCase) {
    // load the reference sequence and the nanopore read
    char *ZymoReference = stString_print("../../signalAlign/tests/test_npReads/ZymoRef.txt");
    FILE *fH = fopen(ZymoReference, "r");
    char *ZymoReferenceSeq = stFile_getLineFromFile(fH);
    char *npReadFile = stString_print("../../signalAlign/tests/test_npReads/ZymoC_ch_1_file1.npRead");
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);
    nanopore_descaleNanoporeRead(npRead);

    // get sequence lengths
    int64_t lX = sequence_correctSeqLength(strlen(ZymoReferenceSeq), event);
    int64_t lY = npRead->nbTemplateEvents;
    char *modelFile = stString_print("../../signalAlign/models/testModel_template.model");
    NanoporeHDP *nHdp = deserialize_nhdp("../../signalAlign/models/templateSingleLevelFixed.nhdp");
    StateMachine *sMt = getHdpMachine(nHdp, modelFile, npRead->templateParams);
    // parameters for pairwise alignment using defaults
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
    p->threshold = 0.1;

    // get anchors using lastz
    stList *anchorPairs = getBlastPairsForPairwiseAlignmentParameters(ZymoReferenceSeq, npRead->twoDread, p);

    // remap and filter
    stList *remappedAnchors = nanopore_remapAnchorPairs(anchorPairs, npRead->templateEventMap);
    stList *filteredRemappedAnchors = filterToRemoveOverlap(remappedAnchors);

    // make Sequences for reference and template events
    Sequence *refSeq = sequence_constructReferenceSequence(lX, ZymoReferenceSeq,
                                                           sequence_getKmer3, sequence_sliceNucleotideSequence2,
                                                           THREE_CYTOSINES, NB_CYTOSINE_OPTIONS, kmer);
    Sequence *templateSeq = sequence_construct2(lY, npRead->templateEvents, sequence_getEvent,
                                                sequence_sliceEventSequence2, event);

    // do alignment of template events
    stList *alignedPairs = getAlignedPairsUsingAnchors(sMt, refSeq, templateSeq, filteredRemappedAnchors, p,
                                                       diagonalCalculationPosteriorMatchProbs,
                                                       0, 0);
    checkAlignedPairs(testCase, alignedPairs, lX, lY);

    // for ch1_file1 template there should be this many aligned pairs with banding
    //st_uglyf("got %lld alignedPairs with anchors\n", stList_length(alignedPairs));
    CuAssertTrue(testCase, stList_length(alignedPairs) == 1213);

    // clean
    pairwiseAlignmentBandingParameters_destruct(p);
    nanopore_nanoporeReadDestruct(npRead);
    sequence_sequenceDestroy(refSeq);
    sequence_sequenceDestroy(templateSeq);
    stList_destruct(alignedPairs);
    stateMachine_destruct(sMt);
}

static void test_sm3Hdp_getAlignedPairsWithBanding_withReplacement(CuTest *testCase) {
    // load the reference sequence and the nanopore read
    char *ZymoReference = stString_print("../../signalAlign/tests/test_npReads/ZymoRef.txt");
    FILE *fH = fopen(ZymoReference, "r");
    char *ZymoReferenceSeq = stFile_getLineFromFile(fH);
    char *CtoM_referenceSeq = stString_replace(ZymoReferenceSeq, "C", "X");
    // npRead
    char *npReadFile = stString_print("../../signalAlign/tests/test_npReads/ZymoC_ch_1_file1.npRead");
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);
    nanopore_descaleNanoporeRead(npRead);

    // get sequence lengths
    int64_t lX = sequence_correctSeqLength(strlen(ZymoReferenceSeq), event);
    int64_t lY = npRead->nbTemplateEvents;
    char *modelFile = stString_print("../../signalAlign/models/testModel_template.model");
    NanoporeHDP *nHdp = deserialize_nhdp("../../signalAlign/models/templateSingleLevelFixed.nhdp");
    StateMachine *sMt = getHdpMachine(nHdp, modelFile, npRead->templateParams);

    // parameters for pairwise alignment using defaults
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
    p->threshold = 0.1;

    // get anchors using lastz
    stList *anchorPairs = getBlastPairsForPairwiseAlignmentParameters(ZymoReferenceSeq, npRead->twoDread, p);

    // remap and filter
    stList *remappedAnchors = nanopore_remapAnchorPairs(anchorPairs, npRead->templateEventMap);
    stList *filteredRemappedAnchors = filterToRemoveOverlap(remappedAnchors);

    // make Sequences for reference and template events
    Sequence *refSeq = sequence_construct2(lX, ZymoReferenceSeq, sequence_getKmer3,
                                           sequence_sliceNucleotideSequence2, kmer);
    Sequence *CtoM_refSeq = sequence_construct2(lX, CtoM_referenceSeq, sequence_getKmer3,
                                                sequence_sliceNucleotideSequence2, kmer);
    Sequence *templateSeq = sequence_construct2(lY, npRead->templateEvents, sequence_getEvent,
                                                sequence_sliceEventSequence2, event);

    // do alignment of template events
    stList *alignedPairs = getAlignedPairsUsingAnchors(sMt, refSeq, templateSeq, filteredRemappedAnchors, p,
                                                       diagonalCalculationPosteriorMatchProbs,
                                                       0, 0);
    p->threshold = 0.01;
    stList *alignedPairs2 = getAlignedPairsUsingAnchors(sMt, CtoM_refSeq, templateSeq, filteredRemappedAnchors, p,
                                                        diagonalCalculationPosteriorMatchProbs,
                                                        0, 0);


    checkAlignedPairs(testCase, alignedPairs, lX, lY);
    checkAlignedPairsForEchelon(testCase, alignedPairs2, lX, lY);

    // for ch1_file1 template there should be this many aligned pairs with banding
    st_uglyf("got %lld alignedPairs on the normal sequence\n", stList_length(alignedPairs));
    st_uglyf("got %lld alignedPairs on the methyl sequence\n", stList_length(alignedPairs2));
    //CuAssertTrue(testCase, stList_length(alignedPairs) == 2876);
    //CuAssertTrue(testCase, stList_length(alignedPairs2) == 2887);

    // clean
    pairwiseAlignmentBandingParameters_destruct(p);
    nanopore_nanoporeReadDestruct(npRead);
    sequence_sequenceDestroy(refSeq);
    sequence_sequenceDestroy(templateSeq);
    //stList_destruct(alignedPairs);
    stList_destruct(alignedPairs2);
    stateMachine_destruct(sMt);
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
    StateMachine *sM = getStrawManStateMachine3(templateLookupTableFile);

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
    StateMachine *sM = getStrawManStateMachine3(model);
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

static void test_continuousPairHmm_em(CuTest *testCase) {
    char *referencePath = stString_print("../../cPecan/tests/test_npReads/ZymoRef.txt");
    FILE *fH = fopen(referencePath, "r");
    char *ZymoReferenceSeq = stFile_getLineFromFile(fH);

    // load the npRead
    char *npReadFile = stString_print("../../cPecan/tests/test_npReads/ZymoC_ch_1_file1.npRead");
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);

    // get sequence lengths
    int64_t lX = sequence_correctSeqLength(strlen(ZymoReferenceSeq), event);
    int64_t lY = npRead->nbTemplateEvents;

    // start with random model
    double pLikelihood = -INFINITY;

    // load stateMachine and model
    const char *model = "../../signalAlign/models/testModel_template.model";
    if (!stFile_exists(model)) {
        st_errAbort("Didn't find model file %s\n", model);
    }

    // load it into the stateMachine
    StateMachine *sM = getSM3_descaled(model, npRead->templateParams);
    //StateMachine *sM = getStrawManStateMachine3(model);

    //emissions_signal_scaleModel(sM, npRead->templateParams.scale, npRead->templateParams.shift,
    //                            npRead->templateParams.var, npRead->templateParams.scale_sd,
    //                            npRead->templateParams.var_sd);

    // parameters for pairwise alignment using defaults
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

    for (int64_t iter = 0; iter < 10; iter++) {
        Hmm *hmm = continuousPairHmm_construct(0.001, 0.001, 3, NUM_OF_KMERS, threeState,
                                               npRead->templateParams.scale,
                                               npRead->templateParams.shift,
                                               npRead->templateParams.var);

        //print_stateMachine_transitions(sM);
        ContinuousPairHmm *cpHmm = (ContinuousPairHmm *)hmm;

        implantModelFromStateMachineIntoHmm(sM, cpHmm);
        // E step
        // get anchors using lastz
        stList *anchorPairs = getBlastPairsForPairwiseAlignmentParameters(ZymoReferenceSeq, npRead->twoDread, p);

        // remap and filter
        stList *remappedAnchors = nanopore_remapAnchorPairs(anchorPairs, npRead->templateEventMap);
        stList *filteredRemappedAnchors = filterToRemoveOverlap(remappedAnchors);

        // make Sequences for reference and template events
        //Sequence *refSeq = sequence_construct2(lX, ZymoReferenceSeq, sequence_getKmer3,
        //                                       sequence_sliceNucleotideSequence2, kmer);
        Sequence *refSeq = sequence_constructReferenceSequence(lX, ZymoReferenceSeq,
                                                               sequence_getKmer3, sequence_sliceNucleotideSequence2,
                                                               THREE_CYTOSINES, NB_CYTOSINE_OPTIONS, kmer);
        Sequence *templateSeq = sequence_construct2(lY, npRead->templateEvents, sequence_getEvent,
                                                    sequence_sliceEventSequence2, event);

        // get expectations
        getExpectationsUsingAnchors(sM, (Hmm *)cpHmm, refSeq, templateSeq, filteredRemappedAnchors,
                                    p, diagonalCalculation_Expectations, 0, 0);

        // normalize
        continuousPairHmm_normalize((Hmm *)cpHmm);

        //st_uglyf("->->-> Got expected likelihood %f for iteration %" PRIi64 "\n", cpHmm->baseHmm.likelihood, iter);
        // M step
        continuousPairHmm_loadTransitionsIntoStateMachine(sM, (Hmm *)cpHmm);
        continuousPairHmm_loadEmissionsIntoStateMachine(sM, (Hmm *)cpHmm);
        //emissions_signal_scaleEmissions(sM, npRead->templateParams.scale, npRead->templateParams.shift,
        //                                npRead->templateParams.var);

        //continuousPairHmm_dumpToFile((Hmm *)cpHmm, priorHmm);

        // Tests
        if (iter > 1) {
            assert(pLikelihood <= cpHmm->baseHmm.likelihood * 0.95);
            CuAssertTrue(testCase, pLikelihood <= cpHmm->baseHmm.likelihood * 0.85);
        }

        // update
        pLikelihood = cpHmm->baseHmm.likelihood;

        // per iteration clean up
        continuousPairHmm_destruct((Hmm *)cpHmm);
        sequence_sequenceDestroy(refSeq);
        sequence_sequenceDestroy(templateSeq);
        stList_destruct(filteredRemappedAnchors);
    }
    //stFile_rmrf(priorHmm);
    nanopore_nanoporeReadDestruct(npRead);
    pairwiseAlignmentBandingParameters_destruct(p);
    stateMachine_destruct(sM);
}

static void test_hdpHmm_em(CuTest *testCase) {
    // load the reference sequence
    char *referencePath = stString_print("../../signalAlign/tests/test_npReads/ZymoRef.txt");
    FILE *fH = fopen(referencePath, "r");
    char *ZymoReferenceSeq = stFile_getLineFromFile(fH);

    // parameters for pairwise alignment using defaults
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

    // load the npRead
    char *npReadFile = stString_print("../../signalAlign/tests/test_npReads/ZymoC_ch_1_file1.npRead");
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);

    nanopore_descaleNanoporeRead(npRead);

    // get sequence lengths
    int64_t lX = sequence_correctSeqLength(strlen(ZymoReferenceSeq), event);
    int64_t lY = npRead->nbTemplateEvents;

    double pLikelihood = -INFINITY;

    Hmm *hdpHmm = hdpHmm_constructEmpty(0.0001, 3, threeStateHdp, p->threshold);

    //continuousPairHmm_randomize(hdpHmm);
    hmmDiscrete_randomizeTransitions(hdpHmm);

    // make initial NanoporeHdp from alignment
    char *modelFile = stString_print("../../signalAlign/models/testModel_template.model");
    char *alignmentFile = stString_print("../../signalAlign/tests/test_alignments/simple_alignment.tsv");
    char *strand = "t";
    NanoporeHDP *nHdp = deserialize_nhdp("../../signalAlign/models/templateSingleLevelFixed.nhdp");
    execute_nhdp_gibbs_sampling(nHdp, 500, 10000, 100, TRUE);
    //finalize_nhdp_distributions(nHdp);
    //NanoporeHDP *nHdp = flat_hdp_model("ACEGOT", SYMBOL_NUMBER_NO_N, KMER_LENGTH,
    //                                   5.0, 0.5, // base_gamma, leaf_gamma
    //                                   0.0, 100.0, 100, modelFile);
    //update_nhdp_from_alignment_with_filter(nHdp, alignmentFile, FALSE, strand);

    // make statemachine
    //StateMachine *sMt = getHdpStateMachine3(nHdp, modelFile);
    StateMachine *sMt = getHdpMachine(nHdp, modelFile, npRead->templateParams);

    // load (random) transitions into stateMachine
    continuousPairHmm_loadEmissionsIntoStateMachine(sMt, hdpHmm);

    // close hmm
    hdpHmm_destruct(hdpHmm);

    for (int64_t iter = 0; iter < 10; iter++) {
        st_uglyf("at start of iteration %lld\n", iter);
        Hmm *hmmExpectations = hdpHmm_constructEmpty(0.0001, 3, threeStateHdp, p->threshold);
        // E step
        // get anchors using lastz
        stList *anchorPairs = getBlastPairsForPairwiseAlignmentParameters(ZymoReferenceSeq, npRead->twoDread, p);

        // remap and filter
        stList *remappedAnchors = nanopore_remapAnchorPairs(anchorPairs, npRead->templateEventMap);
        stList *filteredRemappedAnchors = filterToRemoveOverlap(remappedAnchors);

        // make Sequences for reference and template events
        Sequence *refSeq = sequence_construct2(lX, ZymoReferenceSeq, sequence_getKmer3,
                                               sequence_sliceNucleotideSequence2, kmer);
        Sequence *templateSeq = sequence_construct2(lY, npRead->templateEvents, sequence_getEvent,
                                                    sequence_sliceEventSequence2, event);

        getExpectationsUsingAnchors(sMt, hmmExpectations, refSeq, templateSeq, filteredRemappedAnchors,
                                    p, diagonalCalculation_Expectations, 0, 0);

        // norm
        //continuousPairHmm_normalize(hmmExpectations);
        hmmDiscrete_normalizeTransitions(hmmExpectations);
        st_uglyf("normalized transitions\n");
        st_uglyf("->->-> Got expected likelihood %f for iteration %" PRIi64 "\n", hmmExpectations->likelihood, iter);

        // M step
        // dump hmm (expectations) to disk
        char *tempHdpHmmFile = stString_print("../../signalAlign/tests/test_hdp/tempHdpHmm.hmm");
        FILE *fH = fopen(tempHdpHmmFile, "w");
        hdpHmm_writeToFile(hmmExpectations, fH);
        fclose(fH);

        // load into nHdp
        updateStateMachineHDP(tempHdpHmmFile, sMt);
        remove(tempHdpHmmFile);

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
        sequence_sequenceDestroy(refSeq);
        sequence_sequenceDestroy(templateSeq);
        stList_destruct(filteredRemappedAnchors);
    }
    nanopore_nanoporeReadDestruct(npRead);
    pairwiseAlignmentBandingParameters_destruct(p);
    stateMachine_destruct(sMt);
}

CuSuite *highOrderPairwiseAlignerTestSuite(void) {
    CuSuite *suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_findPotentialMethylation);
    SUITE_ADD_TEST(suite, test_pathLegalTransitions);
    SUITE_ADD_TEST(suite, test_methylPermutations);
    SUITE_ADD_TEST(suite, test_substitutedKmers);
    SUITE_ADD_TEST(suite, test_hdCellConstruct);
    SUITE_ADD_TEST(suite, test_hdCellConstructWorstCase);
    SUITE_ADD_TEST(suite, test_dpDiagonal);
    SUITE_ADD_TEST(suite, test_dpMatrix);
    SUITE_ADD_TEST(suite, test_getKmer4);
    SUITE_ADD_TEST(suite, test_sm3_diagonalDPCalculations);
    //SUITE_ADD_TEST(suite, test_sm3Hdp_diagonalDPCalculations); // known fail
    SUITE_ADD_TEST(suite, test_strawMan_getAlignedPairsWithBanding);
    SUITE_ADD_TEST(suite, test_strawMan_getDescaledAlignedPairsWithBanding);
    SUITE_ADD_TEST(suite, test_sm3Hdp_getAlignedPairsWithBanding);
    //SUITE_ADD_TEST(suite, test_sm3Hdp_getAlignedPairsWithBanding_withReplacement);
    SUITE_ADD_TEST(suite, test_hdpHmmWithoutAssignments);
    SUITE_ADD_TEST(suite, test_continuousPairHmm);
    //SUITE_ADD_TEST(suite, test_continuousPairHmm_em); // TODO make 'sample' method
    //SUITE_ADD_TEST(suite, test_hdpHmm_em);

    return suite;
}
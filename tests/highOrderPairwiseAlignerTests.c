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

// legal path test
static void test_pathLegalTransitions(CuTest *testCase) {
    Path *path1 = path_construct("AAETTT", 3);
    Path *path2 = path_construct("AETTTC", 3);
    Path *path3 = path_construct("AETTTC", 2);
    Path *path4 = path_construct("ACTTTC", 3);
    CuAssertTrue(testCase, path_checkLegal(path1, path2) == TRUE);
    CuAssertTrue(testCase, path_checkLegal(path1, path3) == FALSE);
    CuAssertTrue(testCase, path_checkLegal(path1, path4) == FALSE);
    path_destruct(path1);
    path_destruct(path2);
    path_destruct(path3);
    path_destruct(path4);
}


// test make potentially methylated kmers
static void test_methylPermutations(CuTest *testCase) {
    CuAssertTrue(testCase, path_getMehtylPermutations(0) == NULL);
    for (int64_t i = 1; i < 6; i++) {
        stList *patterns = path_getMehtylPermutations(i);
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
    HDCell *cell = hdCell_construct(ambigKmer, 3);
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
    HDCell *cell = hdCell_construct(ambigKmer, 3);
    Path *path = hdCell_getPath(cell, 0);
    Path *path2 = hdCell_getPath(cell, 728);
    CuAssertIntEquals(testCase, 729, (int )cell->numberOfPaths);
    CuAssertStrEquals(testCase, path->kmer, "CCCCCC");
    CuAssertStrEquals(testCase, path2->kmer, "OOOOOO");
    hdCell_destruct(cell);
}

Sequence *makeTestKmerSequence() {
    char *s = "ATGXAXA"; // has 2 6mers
    int64_t lX = sequence_correctSeqLength(strlen(s), kmer);
    Sequence *seq = sequence_construct(lX, s, sequence_getKmer3, kmer);
    return seq;
}

Sequence *makeKmerSequence(char *nucleotides) {
    int64_t lX = sequence_correctSeqLength(strlen(nucleotides), kmer);
    Sequence *seq = sequence_construct(lX, nucleotides, sequence_getKmer3, kmer);
    return seq;
}

static void test_dpDiagonal(CuTest *testCase) {
    // load model and make stateMachine
    char *modelFile = stString_print("../../signalAlign/models/template_median68pA.model");
    char *canonicalAlphabet = "ACGT\0";

    NanoporeHDP *nHdp = flat_hdp_model(canonicalAlphabet, 4, KMER_LENGTH, 1.0, 1.0, 40.0, 100.0, 100, modelFile);
    StateMachine *sM = getHdpStateMachine3(nHdp);
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
            totalProb = logAdd(totalProb, (1/c1->numberOfPaths) * path1->cells[s]);
            totalProb = logAdd(totalProb, (1/c1->numberOfPaths) * path2->cells[s]);
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
    Sequence *SsX = sequence_construct(lX, sX, sequence_getKmer3, kmer);
    Sequence *SsY = sequence_construct(lY, sY, sequence_getEvent, event);

    // make stateMachine, forward and reverse DP matrices and banding stuff
    char *modelFile = stString_print("../../signalAlign/models/template_median68pA.model");
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
    //dpDoagonal_setValues(dpMatrix_getDiagonal(dpMatrixForward, 0), sM, sM->startStateProb);
    //dpDoagonal_setValues(dpMatrix_getDiagonal(dpMatrixBackward, lX + lY), sM, sM->endStateProb);


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
        st_uglyf("path %lld's kmer is %s\n", p, path->kmer);
        double *cells = path_getCell(path);
        for (int64_t s = 0; s < path->stateNumber; s++) {
            st_uglyf("s:%lld - cells:%f\n", s, cells[s]);
        }
        totalProbForward = logAdd(totalProbForward, cell_dotProduct2(path->cells, sM, sM->endStateProb));
    }
    st_uglyf("totalprobforward: %f\n", totalProbForward);
    breakpoint("HERE!");

    double totalProbBackward = LOG_ZERO;
    DpDiagonal *dpDiagonalB = dpMatrix_getDiagonal(dpMatrixBackward, 0);
    HDCell *cellB = dpDiagonal_getCell(dpDiagonalB, 0);
    for (int64_t p = 0; p < cellB->numberOfPaths; p++) {
        Path *path = hdCell_getPath(cellB, p);
        for (int64_t s = 0; s < path->stateNumber; s++) {
            st_uglyf("s:%lld - cells:%f\n", s, path->cells[s]);
        }
        totalProbBackward = logAdd(totalProbBackward, cell_dotProduct2(path->cells, sM, sM->startStateProb));
    }
    st_uglyf("totalprobBackward: %f\n", totalProbBackward);
    //double totalProbBackward = cell_dotProduct2(
    //        dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixBackward, 0), 0), sM, sM->startStateProb);
    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.001);

    //st_logInfo("Total forward and backward prob %f %f\n", (float) totalProbForward, (float) totalProbBackward);
    CuAssertTrue(testCase, FALSE);
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
        st_uglyf("Pair %f %" PRIi64 " %" PRIi64 "\n", (float) stIntTuple_get(pair, 0) / PAIR_ALIGNMENT_PROB_1, x, y);
        //CuAssertTrue(testCase, stSortedSet_search(alignedPairsSet, stIntTuple_construct2(x, y)) != NULL);
    }
    //CuAssertIntEquals(testCase, 8, (int) stList_length(alignedPairs));

    // clean up
    stateMachine_destruct(sM);
    sequence_sequenceDestroy(SsX);
    sequence_sequenceDestroy(SsY);
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
    //SUITE_ADD_TEST(suite, test_sm3_diagonalDPCalculations);

    return suite;
}
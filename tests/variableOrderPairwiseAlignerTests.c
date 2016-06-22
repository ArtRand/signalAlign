// Tests for high order HMM

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include "stateMachine.h"
#include "CuTest.h"
#include "pairwiseAligner.h"
//#include "discreteHmm.h"
//#include "continuousHmm.h"
#include "randomSequences.h"

// helper functions
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

static inline char *homopolymer(char base, int64_t repeat) {
    char *homopolymer = (char *)st_malloc(sizeof(char) * repeat);
    for (int64_t i = 0; i < repeat; i++) {
        homopolymer[i] = base;
    }
    return homopolymer;
}

//static void breakpoint(char *msg) {
//    st_errAbort("Breakpoint %s\n", msg);
//}


//////////////////////////////////////////////Function Tests/////////////////////////////////////////////////////////
static void test_getKmerWithBoundsCheck(CuTest *testCase) {
    char *nucleotides = stString_print("ATGCATAGC");
    Sequence *sX = sequence_constructKmerSequence(sequence_correctSeqLength(strlen(nucleotides), kmer),nucleotides,
                                                  sequence_getKmer, sequence_sliceNucleotideSequence,
                                                  THREE_CYTOSINES, NB_CYTOSINE_OPTIONS, kmer);
    //ATGCAT
    // TGCATA
    //  GCATAG
    //   CATAGC
    void *checkMinusOne = sequence_getKmerWithBoundsCheck(sX, -1);
    CuAssertTrue(testCase, checkMinusOne == NULL);
    void *checkGreaterThanLength = sequence_getKmerWithBoundsCheck(sX, (sX->length + 1));
    CuAssertTrue(testCase, checkGreaterThanLength == NULL);
    void *firstKmer = sequence_getKmerWithBoundsCheck(sX, 0);
    CuAssert(testCase, "", firstKmer != NULL);
    CuAssertStrEquals(testCase, firstKmer, nucleotides);
}

static void test_findDegeneratePositions(CuTest *testCase) {
    stList *noMethyls = path_findDegeneratePositions("ATGCAC");
    CuAssertTrue(testCase, stList_length(noMethyls) == 0);
    stList *methyls = path_findDegeneratePositions("ATGXAX");
    CuAssertTrue(testCase, stList_length(methyls) == 2);
    int64_t methyl_1 = *(int64_t *)stList_get(methyls, 0);
    int64_t methyl_2 = *(int64_t *)stList_get(methyls, 1);
    CuAssertTrue(testCase, methyl_1 == 3);
    CuAssertTrue(testCase, methyl_2 == 5);
    stList_destruct(noMethyls);
    stList_destruct(methyls);
}

static void test_checkPathConstruct(CuTest *testCase) {
    Path *path = path_construct(NULL, 3);
    CuAssert(testCase, "",path->kmer == NULL);
    path_destruct(path);
    char *tS = getRandomSequence(KMER_LENGTH);
    path = path_construct(tS, 3);
    CuAssertStrEquals(testCase, tS, path->kmer);
    CuAssertIntEquals(testCase, path->stateNumber, 3);
    double *cells = path_getCell(path);
    for (int64_t i = 0; i < 3; i++) {
        CuAssertDblEquals(testCase, cells[i], 0.0, 0.0);
    }
    // put some numbers in the cells
    for (int64_t i = 0; i < 3; i++) {
        cells[i] = (double )i;
        CuAssertDblEquals(testCase, cells[i], (double )i, 0.0);
    }

    Path *copy = path_construct(tS, 3);
    path_copyCells(path, copy);
    double *copyCells = path_getCell(copy);
    for (int64_t i = 0; i < 3; i++) {
        CuAssertDblEquals(testCase, copyCells[i], (double )i, 0.0);
    }

    path_destruct(path);
    path_destruct(copy);
}

static void test_pathLegalTransitions(CuTest *testCase) {
    Path *path0 = path_construct(NULL, 3);
    Path *path1 = path_construct("AAETTT", 3);
    Path *path2 = path_construct("AETTTC", 3);
    Path *path3 = path_construct("AETTTC", 2);
    Path *path4 = path_construct("ACTTTC", 3);
    CuAssertTrue(testCase, path_checkLegal(path0, path1) == TRUE);
    CuAssertTrue(testCase, path_checkLegal(path1, path0) == TRUE);
    CuAssertTrue(testCase, path_checkLegal(path1, path2) == TRUE);
    CuAssertTrue(testCase, path_checkLegal(path1, path3) == FALSE);
    CuAssertTrue(testCase, path_checkLegal(path1, path4) == FALSE);
    path_destruct(path0);
    path_destruct(path1);
    path_destruct(path2);
    path_destruct(path3);
    path_destruct(path4);
}

static void test_checkKmerLengths(CuTest *testCase, stList *kmers, int64_t correctLength) {
    for (int64_t i = 0; i < stList_length(kmers); i++) {
        char *kmer = stList_get(kmers, i);
        CuAssert(testCase, "", strlen(kmer) == correctLength);
    }
}

static void test_firstAndLastKmerOfPermutation(CuTest *testCase, stList *permutations, int64_t length,
                                               char *options) {
    char *firstKmer = stList_get(permutations, 0);
    char *lastKmer = stList_get(permutations, (stList_length(permutations) - 1));
    CuAssertStrEquals(testCase, firstKmer, homopolymer(options[0], length));
    CuAssertStrEquals(testCase, lastKmer, homopolymer(options[strlen(options) - 1], length));
}

static void test_checkPermutationsP(CuTest *testCase, char *options) {
    CuAssertTrue(testCase, path_listPotentialKmers(0, strlen(options), options) == NULL);
    for (int64_t i = 1; i < KMER_LENGTH + 1; i++) {
        stList *patterns = path_listPotentialKmers(i, strlen(options), options);
        test_checkKmerLengths(testCase, patterns, i);
        test_firstAndLastKmerOfPermutation(testCase, patterns, i, options);
        int64_t correctLength = intPow(strlen(options), i);
        CuAssertTrue(testCase, stList_length(patterns) == correctLength);
        stList_destruct(patterns);
    }
}

static void test_Permutations(CuTest *testCase) {
    test_checkPermutationsP(testCase, THREE_CYTOSINES);
    test_checkPermutationsP(testCase, TWO_CYTOSINES);
    test_checkPermutationsP(testCase, CANONICAL_NUCLEOTIDES);
}

static void test_getKmerIndex(CuTest *testCase) {
    stList *kmers = path_listPotentialKmers(KMER_LENGTH, strlen(ALL_BASES), ALL_BASES);
    for (int64_t i = 0; i < stList_length(kmers); i++) {
        char *kmer = stList_get(kmers, i);
        int64_t index = emissions_discrete_getKmerIndexFromPtr(kmer);
        CuAssertIntEquals(testCase, index, i);
    }
}

static void test_substitutedKmers(CuTest *testCase) {
    char *ambigKmer = "ATGXAX";
    stList *positions = path_findDegeneratePositions(ambigKmer);
    char *pattern = "CE";
    char *kmer = hdCell_getSubstitutedKmer(positions, 2, pattern, ambigKmer);
    CuAssertTrue(testCase, stList_length(positions) == 2);
    CuAssertStrEquals(testCase, kmer, "ATGCAE");
    stList_destruct(positions);
    free(kmer);
}



CuSuite *variableOrderPairwiseAlignerTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_findDegeneratePositions);
    SUITE_ADD_TEST(suite, test_checkPathConstruct);
    SUITE_ADD_TEST(suite, test_pathLegalTransitions);
    SUITE_ADD_TEST(suite, test_Permutations);
    SUITE_ADD_TEST(suite, test_substitutedKmers);
    SUITE_ADD_TEST(suite, test_getKmerIndex);
    SUITE_ADD_TEST(suite, test_getKmerWithBoundsCheck);
    return suite;
}
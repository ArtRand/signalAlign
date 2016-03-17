/*
 * pairwiseAligner.c
 *
 *  Created on: 1 Mar 2012
 *      Author: benedictpaten
 */
//This is being included to make popen work!
#define _XOPEN_SOURCE 500

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <stdbool.h>
#include <inttypes.h>
#include <stdint.h>
#include <discreteHmm.h>
#include "nanopore.h"
#include "sonLib.h"
#include "bioioC.h"
#include "pairwiseAlignment.h"
#include "pairwiseAligner.h"
#include "continuousHmm.h"
#include "stateMachine.h"
#include "emissionMatrix.h"




/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Diagonal
//Structure for working with x-y diagonal of dp matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////

const char *PAIRWISE_ALIGNMENT_EXCEPTION_ID = "PAIRWISE_ALIGNMENT_EXCEPTION";

Diagonal diagonal_construct(int64_t xay, int64_t xmyL, int64_t xmyR) {
    if ((xay + xmyL) % 2 != 0 || (xay + xmyR) % 2 != 0 || xmyL > xmyR) {
        stThrowNew(PAIRWISE_ALIGNMENT_EXCEPTION_ID,
                "Attempt to create diagonal with invalid coordinates: xay %" PRIi64 " xmyL %" PRIi64 " xmyR %" PRIi64 "",
                xay, xmyL, xmyR);
    }
    Diagonal diagonal;
    diagonal.xay = xay;
    diagonal.xmyL = xmyL;
    diagonal.xmyR = xmyR;
    assert(xmyL <= xmyR);
    assert(xay >= 0);
    return diagonal;
}

inline int64_t diagonal_getXay(Diagonal diagonal) {
    return diagonal.xay;
}

inline int64_t diagonal_getMinXmy(Diagonal diagonal) {
    return diagonal.xmyL;
}

inline int64_t diagonal_getMaxXmy(Diagonal diagonal) {
    return diagonal.xmyR;
}

inline int64_t diagonal_getWidth(Diagonal diagonal) {
    return (diagonal.xmyR - diagonal.xmyL) / 2 + 1;
}

inline int64_t diagonal_getXCoordinate(int64_t xay, int64_t xmy) {
    assert((xay + xmy) % 2 == 0);
    return (xay + xmy) / 2;
}

inline int64_t diagonal_equals(Diagonal diagonal1, Diagonal diagonal2) {
    return diagonal1.xay == diagonal2.xay && diagonal1.xmyL == diagonal2.xmyL && diagonal1.xmyR == diagonal2.xmyR;
}

inline int64_t diagonal_getYCoordinate(int64_t xay, int64_t xmy) {
    assert((xay - xmy) % 2 == 0);
    return (xay - xmy) / 2;
}

inline char *diagonal_getString(Diagonal diagonal) {
    return stString_print("Diagonal, xay: %" PRIi64 " xmyL %" PRIi64 ", xmyR: %" PRIi64 "", diagonal_getXay(diagonal),
            diagonal_getMinXmy(diagonal), diagonal_getMaxXmy(diagonal));
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Band Iterator
//Iterator for walking along x+y diagonals in banded fashion
//(using a set of anchor constraints)
/////////////////////////////////////////////////////////////////////////////////////////////////////////

struct _band {
    Diagonal *diagonals;
    int64_t lXalY;
};

static int64_t band_avoidOffByOne(int64_t xay, int64_t xmy) {
    return (xay + xmy) % 2 == 0 ? xmy : xmy + 1;
}

static void band_setCurrentDiagonalP(int64_t *xmy, int64_t i, int64_t j, int64_t k) {
    if (i < j) {
        *xmy += (int64_t) (2 * ((int64_t) j - i) * (int64_t) k);
    }
}

static Diagonal band_setCurrentDiagonal(int64_t xay, int64_t xL, int64_t yL, int64_t xU, int64_t yU) {
    int64_t xmyL = xL - yL;
    int64_t xmyR = xU - yU;

    assert(xay >= xL + yU);
    assert(xay <= xU + yL);

    //Avoid in-between undefined x,y coordinate positions when intersecting xay and xmy.
    xmyL = band_avoidOffByOne(xay, xmyL);
    xmyR = band_avoidOffByOne(xay, xmyR);

    //Bound the xmy coordinates by the xL, yL and xU, yU band boundaries
    band_setCurrentDiagonalP(&xmyL, diagonal_getXCoordinate(xay, xmyL), xL, 1);
    band_setCurrentDiagonalP(&xmyL, yL, diagonal_getYCoordinate(xay, xmyL), 1);
    band_setCurrentDiagonalP(&xmyR, xU, diagonal_getXCoordinate(xay, xmyR), -1);
    band_setCurrentDiagonalP(&xmyR, diagonal_getYCoordinate(xay, xmyR), yU, -1);

    return diagonal_construct(xay, xmyL, xmyR);
}

static int64_t band_boundCoordinate(int64_t z, int64_t lZ) {
    return z < 0 ? 0 : (z > lZ ? lZ : z);
}

Band *band_construct(stList *anchorPairs, int64_t lX, int64_t lY, int64_t expansion) {
    //Prerequisities
    assert(lX >= 0);
    assert(lY >= 0);
    assert(expansion % 2 == 0);

    Band *band = st_malloc(sizeof(Band));
    band->diagonals = st_malloc(sizeof(Diagonal) * (lX + lY + 1));
    band->lXalY = lX + lY;

    //Now initialise the diagonals
    int64_t anchorPairIndex = 0;
    int64_t xay = 0;
    int64_t pxay = 0, pxmy = 0;
    int64_t nxay = 0, nxmy = 0;
    int64_t xL = 0, yL = 0, xU = 0, yU = 0;

    while (xay <= band->lXalY) {
        band->diagonals[xay] = band_setCurrentDiagonal(xay, xL, yL, xU, yU);
        if (nxay == xay++) {
            //The previous diagonals become the next
            pxay = nxay;
            pxmy = nxmy;

            int64_t x = lX, y = lY;
            if (anchorPairIndex < stList_length(anchorPairs)) {
                stIntTuple *anchorPair = stList_get(anchorPairs, anchorPairIndex++);
                x = stIntTuple_get(anchorPair, 0) + 1; //Plus ones, because matrix coordinates are +1 the sequence ones
                y = stIntTuple_get(anchorPair, 1) + 1;

                //Check the anchor pairs
                assert(x > diagonal_getXCoordinate(pxay, pxmy));
                assert(y > diagonal_getYCoordinate(pxay, pxmy));
                assert(x <= lX);
                assert(y <= lY);
                assert(x > 0);
                assert(y > 0);
            }

            nxay = x + y;
            nxmy = x - y;

            //Now call to set the lower and upper x,y coordinates
            xL = band_boundCoordinate(diagonal_getXCoordinate(pxay, pxmy - expansion), lX);
            yL = band_boundCoordinate(diagonal_getYCoordinate(nxay, nxmy - expansion), lY);
            xU = band_boundCoordinate(diagonal_getXCoordinate(nxay, nxmy + expansion), lX);
            yU = band_boundCoordinate(diagonal_getYCoordinate(pxay, pxmy + expansion), lY);
        }
    }

    return band;
}

void band_destruct(Band *band) {
    free(band->diagonals);
    free(band);
}

struct _bandIterator {
    Band *band;
    int64_t index;
};

BandIterator *bandIterator_construct(Band *band) {
    BandIterator *bandIterator = st_malloc(sizeof(BandIterator));
    bandIterator->band = band;
    bandIterator->index = 0;
    return bandIterator;
}

BandIterator *bandIterator_clone(BandIterator *bandIterator) {
    BandIterator *bandIterator2 = st_malloc(sizeof(BandIterator));
    memcpy(bandIterator2, bandIterator, sizeof(BandIterator));
    return bandIterator2;
}

void bandIterator_destruct(BandIterator *bandIterator) {
    free(bandIterator);
}

Diagonal bandIterator_getNext(BandIterator *bandIterator) {
    Diagonal diagonal = bandIterator->band->diagonals[
            bandIterator->index > bandIterator->band->lXalY ? bandIterator->band->lXalY : bandIterator->index];
    if (bandIterator->index <= bandIterator->band->lXalY) {
        bandIterator->index++;
    }
    return diagonal;
}

Diagonal bandIterator_getPrevious(BandIterator *bandIterator) {
    if (bandIterator->index > 0) {
        bandIterator->index--;
    }
    return bandIterator->band->diagonals[bandIterator->index];
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Log Add functions
//Interpolation function for doing log add
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#define logUnderflowThreshold 7.5
#define posteriorMatchThreshold 0.01

static inline double lookup(double x) {
    //return log (exp (x) + 1);
    assert(x >= 0.00f);
    assert(x <= logUnderflowThreshold);
    if (x <= 1.00f)
        return ((-0.009350833524763f * x + 0.130659527668286f) * x + 0.498799810682272f) * x + 0.693203116424741f;
    if (x <= 2.50f)
        return ((-0.014532321752540f * x + 0.139942324101744f) * x + 0.495635523139337f) * x + 0.692140569840976f;
    if (x <= 4.50f)
        return ((-0.004605031767994f * x + 0.063427417320019f) * x + 0.695956496475118f) * x + 0.514272634594009f;
    return ((-0.000458661602210f * x + 0.009695946122598f) * x + 0.930734667215156f) * x + 0.168037164329057f;
}

double logAdd(double x, double y) {
    if (x < y)
        return (x == LOG_ZERO || y - x >= logUnderflowThreshold) ? y : lookup(y - x) + x;
    return (y == LOG_ZERO || x - y >= logUnderflowThreshold) ? x : lookup(x - y) + y;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Sequence Object generalized way to represent a sequence of symbols or measurements (events)
/////////////////////////////////////////////////////////////////////////////////////////////////////////

static double _NULLEVENT[] = {LOG_ZERO, 0};
static double *NULLEVENT = _NULLEVENT;

Sequence *sequence_construct(int64_t length, void *elements, void *(*getFcn)(void *, int64_t), SequenceType type) {
    Sequence *self = malloc(sizeof(Sequence));
    self->type = type;
    self->length = length;
    self->elements = elements;
    self->get = getFcn;
    return self;
}

Sequence *sequence_construct2(int64_t length, void *elements, void *(*getFcn)(void *, int64_t),
                              Sequence *(*sliceFcn)(Sequence *, int64_t, int64_t), SequenceType type) {
    Sequence *self = malloc(sizeof(Sequence));
    self->type = type;
    self->length = length;
    self->elements = elements;
    self->get = getFcn;
    self->sliceFcn = sliceFcn;
    return self;
}

void sequence_padSequence(Sequence *sequence) {
    char *endPadding = "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnn";
    sequence->elements = stString_print("%s%s", sequence->elements, endPadding);
}

/*
Sequence *sequence_sliceNucleotideSequence(Sequence *inputSequence, int64_t start, int64_t sliceLength,
                                           void *(*getFcn)(void *, int64_t)) {
    size_t elementSize = sizeof(char);
    void *elementSlice = (char *)inputSequence->elements + (start * elementSize);
    Sequence *newSequence = sequence_construct(sliceLength, elementSlice, getFcn);
    return newSequence;
}
*/
Sequence *sequence_sliceNucleotideSequence2(Sequence *inputSequence, int64_t start, int64_t sliceLength) {
    size_t elementSize = sizeof(char);
    void *elementSlice = (char *)inputSequence->elements + (start * elementSize);
    Sequence *newSequence = sequence_construct2(sliceLength, elementSlice,
                                                inputSequence->get, inputSequence->sliceFcn, inputSequence->type);
    return newSequence;
}
/*
Sequence *sequence_sliceEventSequence(Sequence *inputSequence, int64_t start, int64_t sliceLength,
                                      void *(*getFcn)(void *, int64_t)) {
    size_t elementSize = sizeof(double);
    void *elementSlice = (char *)inputSequence->elements + ((start * NB_EVENT_PARAMS) * elementSize);
    Sequence *newSequence = sequence_construct(sliceLength, elementSlice, getFcn);
    return newSequence;
}
*/
Sequence *sequence_sliceEventSequence2(Sequence *inputSequence, int64_t start, int64_t sliceLength) {
    size_t elementSize = sizeof(double);
    void *elementSlice = (char *)inputSequence->elements + ((start * NB_EVENT_PARAMS) * elementSize);
    Sequence *newSequence = sequence_construct2(sliceLength, elementSlice,
                                                inputSequence->get, inputSequence->sliceFcn, inputSequence->type);
    return newSequence;
}

void sequence_sequenceDestroy(Sequence *seq) {
    //assert(seq != NULL);
    free(seq);
}

void *sequence_getBase(void *elements, int64_t index) {
    char* n;
    n = "n";
    return index >= 0 ? &(((char *)elements)[index]) : n;
}

void *sequence_getKmer(void *elements, int64_t index) {
    char *n = "n";
    //int64_t i = index;
    return index >= 0 ? &(((char *) elements)[index]) : n;
}

void *sequence_getKmer2(void *elements, int64_t index) {
    if (index < 0){
        return &(((char *) elements)[0]);
    }
    return index > 0 ? &(((char *) elements)[index - 1]) : &(((char *) elements)[index]);
}

void *sequence_getKmer3(void *elements, int64_t index) {
    return index >= 0 ? &(((char *) elements)[index]) : &(((char *) elements)[0]);
}

void *sequence_getKmer4(Sequence *sequence, int64_t index) {
    if ((index >= sequence->length) || (index < 0)) {
        return NULL;
    } else {
        return &(((char *) sequence->elements)[index]);
    }
}


void *sequence_getEvent(void *elements, int64_t index) {
    index = index * NB_EVENT_PARAMS;
    //return index >= 0 ? &(((double *)elements)[index]) : NULL;
    return index >= 0 ? &(((double *)elements)[index]) : NULLEVENT;
}

int64_t sequence_correctSeqLength(int64_t length, SequenceType type) {
    // for trivial case
    if (length == 0) {
        return 0;
    }
    if (length > 0) {
        switch (type) {
            case nucleotide: // nucleotide sequence
                return length;
            case kmer: // event and kmer sequence
            case event:
                return length - (KMER_LENGTH - 1);
        }
    }
    return 0;
}
// Path
// from SO: http://stackoverflow.com/questions/213042/how-do-you-do-exponentiation-in-c
int64_t intPow(int64_t base, int64_t exp) {
    if (exp == 0) {
        return 1;
    } else if (exp % 2) {
        return base * intPow(base, exp - 1);
    } else {
        int64_t tmp = intPow(base, exp / 2);
        return tmp * tmp;
    }
}

Path *path_construct(char *kmer, int64_t stateNumber) {
    Path *path = st_malloc(sizeof(Path));
    path->kmer = kmer;
    path->stateNumber = stateNumber;
    //path->cells = st_malloc(sizeof(double) * stateNumber);
    path->cells = st_calloc(stateNumber, sizeof(double));
    return path;
}

stList *path_findPotentialMethylation(char *kmer) {
    stList *methyls = stList_construct3(0, &free);
    if (kmer == NULL) {
        return methyls;
    }
    for (int64_t i = 0; i < KMER_LENGTH; i++) {
        char n = *(kmer + i);
        if (strchr(CYTOSINE_METHYL_AMBIG, n)) {
            int64_t *methylPosition = (int64_t *)st_malloc(sizeof(int64_t));
            *methylPosition = i;
            stList_append(methyls, methylPosition);
        }
    }
    return methyls;
}

static bool path_checkKmerLegalTransition(char *kmerOne, char *kmerTwo) {
    // checks for xATAGA
    //             ATAGAy
    for (int64_t x = 0; x < (KMER_LENGTH - 1); x++) {
        char x_i = *(kmerOne + (x + 1));
        char x_j = *(kmerTwo + x);
        if (x_i == x_j) {
            continue;
        } else {
            return FALSE;
        }
    }
    return TRUE;
}

bool path_checkLegal(Path *path1, Path *path2) {
    if (path1->stateNumber != path2->stateNumber) {
        return FALSE;
    } else if (path1->kmer == NULL) {
        return TRUE;
    } else {
        return path_checkKmerLegalTransition(path1->kmer, path2->kmer);
    }
}

double *path_getCell(Path *path) {
    return path->cells;
}

void path_copyCells(Path *master, Path *copy) {
    if (master->stateNumber != copy->stateNumber) {
        st_errAbort("path_copyCells: master and copy have different stateNumbers\n");
    }
    for (int64_t i = 0; i < master->stateNumber; i++) {
        copy->cells[i] = master->cells[i];
    }
}


void path_destruct(Path *path) {
    free(path->cells);
    free(path);
}

// adapted from SO:
// http://stackoverflow.com/questions/17506848/print-all-possible-strings-of-length-p-that-can-be-formed-from-the-given-set
static void path_permutePattern(stList *methylList, int64_t currentLength, int64_t methylPositions, char *arr) {
    if (currentLength == methylPositions) {
        arr[methylPositions] = '\0';
        char *kmer = stString_copy(arr);
        stList_append(methylList, kmer);
        return;
    } else {
        for(int64_t i = 0; i < NB_CYTOSINE_OPTIONS; i++) {
            arr[currentLength] = CYTOSINES[i];
            path_permutePattern(methylList, currentLength + 1, methylPositions, arr);
        }
    }
}

stList *path_getMehtylPermutations(int64_t methylPositions) {
    if (methylPositions == 0) {
        return NULL;
    } else {
        stList *methylPatterns = stList_construct3(0, &free);
        char *arr = (char *) malloc(sizeof(char) * methylPositions);
        path_permutePattern(methylPatterns, 0, methylPositions, arr);
        free(arr);
        return methylPatterns;
    }
}

// TODO comment
HDCell *hdCell_construct(void *nucleotideSequence, int64_t stateNumber) {
    char *kmer_i;
    if (nucleotideSequence != NULL) {
        kmer_i = (char *)st_malloc((KMER_LENGTH) * sizeof(char));
        for (int64_t x = 0; x < KMER_LENGTH; x++) {
            kmer_i[x] = *((char *)nucleotideSequence + x);
        }
        kmer_i[KMER_LENGTH] = '\0';
    } else {
        kmer_i = NULL;
    }

    HDCell *cell = st_malloc(sizeof(HDCell));

    stList *methylPositions = path_findPotentialMethylation(kmer_i);  // will return empty list for NULL kmer
    int64_t nbMethylPositions = stList_length(methylPositions);  // NULL kmer will be 0

    cell->numberOfPaths = intPow(NB_CYTOSINE_OPTIONS, nbMethylPositions);

    if (cell->numberOfPaths < 1) {
        st_errAbort("hdCell_construct: got illegal number of paths: %lld", cell->numberOfPaths);
    }

    cell->paths = (Path **)st_malloc(sizeof(Path *) * cell->numberOfPaths); // todo make this calloc?

    if (cell->numberOfPaths > 1) {
        stList *methylationPatterns = path_getMehtylPermutations(nbMethylPositions);
        for (int64_t i = 0; i < cell->numberOfPaths; i++) {
            char *pattern = stList_get(methylationPatterns, i);
            Path *path = path_construct(hdCell_getSubstitutedKmer(methylPositions, nbMethylPositions, pattern,
                                                                  kmer_i), stateNumber);
            cell->paths[i] = path;
        }
    } else {
        char *onePathKmer = stString_copy(kmer_i);
        Path *path = path_construct(onePathKmer, stateNumber);
        cell->paths[0] = path;
    }

    cell->init = TRUE;
    free(kmer_i);

    return cell;
}

Path *hdCell_getPath(HDCell *cell, int64_t pathNumber) {
    if (pathNumber >= cell->numberOfPaths) {
        return NULL;
    }
    return cell->paths[pathNumber];
}

char *hdCell_getSubstitutedKmer(stList *methylPositions, int64_t nbPositions, char *pattern, char *ambigKmer) {
    char *kmer = stString_copy(ambigKmer);

    for (int64_t i = 0; i < nbPositions; i++) {
        int64_t position = *(int64_t *)stList_get(methylPositions, i);
        kmer[position] = pattern[i];
    }
    return kmer;
}

void hdCell_destruct(HDCell *cell) {
    for (int64_t i = 0; i < cell->numberOfPaths; i++) {
        path_destruct(cell->paths[i]);
    }
    free(cell);
}

double hdCell_totalProbability(HDCell *cell1, HDCell *cell2) {
    double totalProbability = LOG_ZERO;
    for (int64_t x = 0; x < cell1->numberOfPaths; x++) {
        Path *path1 = hdCell_getPath(cell1, x);
        for (int64_t y = 0; y < cell2->numberOfPaths; y++) {
            Path *path2 = hdCell_getPath(cell2, y);
            if (stString_eq(path1->kmer, path2->kmer)) {
                totalProbability = logAdd(totalProbability,
                                          cell_dotProduct(path1->cells,
                                                          path2->cells,
                                                          path1->stateNumber));
            }
        }
    }
    return totalProbability;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Cell calculations
//A cell is a set of states associated with an x, y coordinate.
//These functions do the forward/backward calculations for the pairwise
//alignment model.
/////////////////////////////////////////////////////////////////////////////////////////////////////////

static inline void doTransitionForward(double *fromCells, double *toCells,
                                       int64_t from, int64_t to,
                                       double eP, double tP,
                                       void *extraArgs) {
    //st_uglyf("TRANSITION FORWARD %f-->%f eP: %f tP: %f\n", fromCells[from], toCells[to], eP, tP);
    toCells[to] = logAdd(toCells[to], fromCells[from] + (eP + tP));
}

void cell_calculateForward(StateMachine *sM,
                           void *current, void *lower, void *middle, void *upper,
                           void* cX, void* cY, void *extraArgs) {
    sM->cellCalculate(sM, current, lower, middle, upper, cX, cY, doTransitionForward, extraArgs);
}

static inline void doTransitionBackward(double *fromCells, double *toCells,
                                        int64_t from, int64_t to,
                                        double eP, double tP,
                                        void *extraArgs) {
    fromCells[from] = logAdd(fromCells[from], toCells[to] + (eP + tP));
}

void cell_calculateBackward(StateMachine *sM,
                            void *current, void *lower, void *middle, void *upper,
                            void* cX, void* cY, void *extraArgs) {
    sM->cellCalculate(sM, current, lower, middle, upper, cX, cY, doTransitionBackward, extraArgs);
}

double cell_dotProduct(double *cell1, double *cell2, int64_t stateNumber) {
    double totalProb = cell1[0] + cell2[0];
    for (int64_t i = 1; i < stateNumber; i++) {
        totalProb = logAdd(totalProb, cell1[i] + cell2[i]);
    }
    return totalProb;
}

double cell_dotProduct2(double *cell, StateMachine *sM, double (*getStateValue)(StateMachine *, int64_t)) {
    double totalProb = cell[0] + getStateValue(sM, 0);
    for (int64_t i = 1; i < sM->stateNumber; i++) {
        totalProb = logAdd(totalProb, cell[i] + getStateValue(sM, i));
    }
    return totalProb;
}

void cell_updateExpectations(double *fromCells, double *toCells, int64_t from, int64_t to, double eP, double tP,
                             void *extraArgs) {

    //void *extraArgs2[2] = { &totalProbability, hmmExpectations };
    double totalProbability = *((double *) ((void **) extraArgs)[0]);
    HmmDiscrete *hmmExpectations = ((void **) extraArgs)[1];

    int64_t x = hmmExpectations->getElementIndexFcn(((void **) extraArgs)[2]); // this gives you the base/kmer index
    int64_t y = hmmExpectations->getElementIndexFcn(((void **) extraArgs)[3]);

    //Calculate posterior probability of the transition/emission pair
    double p = exp(fromCells[from] + toCells[to] + (eP + tP) - totalProbability);
    hmmExpectations->baseHmm.addToTransitionExpectationFcn((Hmm *)hmmExpectations, from, to, p);

    if(x < hmmExpectations->baseHmm.symbolSetSize && y < hmmExpectations->baseHmm.symbolSetSize) { //Ignore gaps involving Ns.
        hmmExpectations->addToEmissionExpectationFcn((Hmm *)hmmExpectations, to, x, y, p);
    }
}

void cell_signal_updateExpectations(double *fromCells, double *toCells, int64_t from, int64_t to,
                                    double eP, double tP, void *extraArgs) {
    //void *extraArgs2[2] = { &totalProbability, hmmExpectations };
    double totalProbability = *((double *) ((void **) extraArgs)[0]);
    ContinuousPairHmm *hmmExpectations = ((void **) extraArgs)[1];
    if (!hmmExpectations->hasModel) {
        st_errAbort("cell_signal_updateExpectations: HMM needs to have model\n");
    }

    int64_t kmerIndex = hmmExpectations->getElementIndexFcn(((void **) extraArgs)[2]); // this gives you the kmer index
    double eventMean = *(double *)((void **) extraArgs)[3];
    double descaledEventMean = hmmExpectations->getDescaledEvent(
            hmmExpectations, eventMean, *(hmmExpectations->getEventModelEntry((Hmm *)hmmExpectations, kmerIndex)));

    // Calculate posterior probability of the transition/emission pair
    double p = exp(fromCells[from] + toCells[to] + (eP + tP) - totalProbability);

    // update transitions expectation
    hmmExpectations->baseHmm.addToTransitionExpectationFcn((Hmm *)hmmExpectations, from, to, p);
    //
    if (to == match) {
        hmmExpectations->addToEmissionExpectationFcn((Hmm *)hmmExpectations, kmerIndex, descaledEventMean, p);
    }
}

void cell_signal_updateTransAndKmerSkipExpectations2(double *fromCells, double *toCells, int64_t from, int64_t to,
                                                    double eP, double tP, void *extraArgs) {
    //void *extraArgs2[2] = { &totalProbability, hmmExpectations };
    // unpack the extraArgs thing
    double totalProbability = *((double *) ((void **) extraArgs)[0]);
    HdpHmm *hmmExpectations = ((void **) extraArgs)[1];

    char *kmer = (char *)((void **) extraArgs)[2];  // pointer to the position in the sequence (kmer)

    double *event = (double *)((void **) extraArgs)[3];  // pointer to the event mean

    // Calculate posterior probability of the transition/emission pair
    double p = exp(fromCells[from] + toCells[to] + (eP + tP) - totalProbability);

    // update transitions expectation
    hmmExpectations->baseHmm.addToTransitionExpectationFcn((Hmm *)hmmExpectations, from, to, p);

    if ((to == match) & (p >= hmmExpectations->threshold)) {
        //st_uglyf("SENTINAL - adding to expectations kmer %s\n", kmer);
        // add to emissions expectations function here
        hmmExpectations->addToAssignments(hmmExpectations, kmer, event);
    }
}

static void cell_calculateUpdateExpectation(StateMachine *sM,
                                            void *current, void *lower, void *middle, void *upper,
                                            void *cX, void *cY,
                                            void *extraArgs) {
    void *extraArgs2[4] = { ((void **)extraArgs)[0], // &totalProbabability
                            ((void **)extraArgs)[1], // hmmExpectations
                            cX,   // pointer to char in sequence
                            cY }; // pointer to event array, can remove..?
    sM->cellCalculate(sM, current, lower, middle, upper, cX, cY,
                      sM->cellCalculateUpdateExpectations,
                      extraArgs2);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//DpDiagonal
/////////////////////////////////////////////////////////////////////////////////////////////////////////

DpDiagonal *dpDiagonal_construct(Diagonal diagonal, int64_t stateNumber, Sequence *nucleotideSequence) {
    if (nucleotideSequence->type != kmer) {
        st_errAbort("dpDiagonal_construct: got illegal sequence type %i", nucleotideSequence->type);
    }

    DpDiagonal *dpDiagonal = st_malloc(sizeof(DpDiagonal));
    dpDiagonal->diagonal = diagonal;
    dpDiagonal->stateNumber = stateNumber;
    dpDiagonal->sequence = nucleotideSequence;
    dpDiagonal->totalPaths = 0;
    assert(diagonal_getWidth(diagonal) >= 0);

    int64_t nCells = diagonal_getWidth(diagonal);
    dpDiagonal->cells = (HDCell **) st_calloc(nCells, sizeof(HDCell *));

    int64_t xmy = diagonal_getMinXmy(diagonal);
    int64_t maxXmy = diagonal_getMaxXmy(diagonal);
    int64_t xay = diagonal_getXay(diagonal);

    while (xmy <= maxXmy) {
        int64_t x = diagonal_getXCoordinate(xay, xmy);

        //void *k = nucleotideSequence->get(nucleotideSequence->elements, (x - 1));
        void *k = sequence_getKmer4(nucleotideSequence, (x -1));

        HDCell *hdCell = hdCell_construct(k, stateNumber);

        dpDiagonal->totalPaths += hdCell->numberOfPaths;
        dpDiagonal_assignCell(dpDiagonal, hdCell, xmy);

        xmy += 2;
    }
    return dpDiagonal;
}

DpDiagonal *dpDiagonal_clone(DpDiagonal *diagonal) { // TODO make a function that clones without needing the sequence?
    // make empty dpDiagonal
    DpDiagonal *diagonal2 = dpDiagonal_construct(diagonal->diagonal, diagonal->stateNumber, diagonal->sequence);

    // copy the cells (doubles) from the paths
    int64_t nCells = diagonal_getWidth(diagonal2->diagonal);

    if (nCells != diagonal_getWidth(diagonal->diagonal)) {
        st_errAbort("dpDiagonal_clone: error with diagonal width\n");
    }

    for (int64_t c = 0; c < nCells; c++) {
        HDCell *hdCell1 = diagonal->cells[c];
        HDCell *hdCell2 = diagonal2->cells[c];

        if ((hdCell1->init != TRUE) || (hdCell2->init != TRUE)) {
            st_errAbort("dpDiagonal_clone: got uninitialized cells\n");
        }

        for (int64_t p = 0; p < hdCell1->numberOfPaths; p++) {
            Path *path1 = hdCell_getPath(hdCell1, p);
            Path *path2 = hdCell_getPath(hdCell2, p);
            if (path1->stateNumber != path2->stateNumber) {
                st_errAbort("dpDiagonal_clone: paths have different stateNumbers\n");
            }
            path_copyCells(path1, path2);
        }
    }
    return diagonal2;
}

bool dpDiagonal_equals(DpDiagonal *diagonal1, DpDiagonal *diagonal2) {
    if (!diagonal_equals(diagonal1->diagonal, diagonal2->diagonal)) {
        return 0;
    }

    if(diagonal1->stateNumber != diagonal2->stateNumber) {
        return 0;
    }
    // iterate over the cells
    for (int64_t i = 0; i < diagonal_getWidth(diagonal1->diagonal); i++) {
        HDCell *hdCell1 = diagonal1->cells[i];
        HDCell *hdCell2 = diagonal2->cells[i];
        if (hdCell1->numberOfPaths != hdCell2->numberOfPaths) {
            fprintf(stderr, "dpDiagonal_equals: different numberOfPaths\n");
            return 0;
        }
        for (int64_t p = 0; p < hdCell1->numberOfPaths; p++) {
            Path *path1 = hdCell_getPath(hdCell1, p);
            Path *path2 = hdCell_getPath(hdCell2, p);
            if (path1->stateNumber != path2->stateNumber) {
                fprintf(stderr, "dpDiagonal_equals: different stateNumber in paths %lld\n", p);
                return 0;
            }
            for (int64_t s = 0; s < path1->stateNumber; s++) {
                if (path1->cells[s] != path2->cells[s]) {
                    return 0;
                }
            }
        }
    }
    return 1;
}

void dpDiagonal_assignCell(DpDiagonal *dpDiagonal, HDCell *hdCell, int64_t xmy) {
    if (xmy < dpDiagonal->diagonal.xmyL || xmy > dpDiagonal->diagonal.xmyR) {
        st_errAbort("dpDiagonal_assignCell: cannot put a diagonal at those coordinates\n");
    }
    assert((diagonal_getXay(dpDiagonal->diagonal) + xmy) % 2 == 0);

    if (((xmy - dpDiagonal->diagonal.xmyL) / 2) >= diagonal_getWidth(dpDiagonal->diagonal)) {
        st_errAbort("dpDiagonal_assignCell: cannot put a diagonal at those coordinates OUT OF BOUNDS\n");
    }

    dpDiagonal->cells[(xmy - dpDiagonal->diagonal.xmyL) / 2] = hdCell;
}

HDCell *dpDiagonal_getCell(DpDiagonal *dpDiagonal, int64_t xmy) {
    if (xmy < dpDiagonal->diagonal.xmyL || xmy > dpDiagonal->diagonal.xmyR) {
        return NULL;
    }
    assert((diagonal_getXay(dpDiagonal->diagonal) + xmy) % 2 == 0);
    return &(*dpDiagonal->cells[(xmy - dpDiagonal->diagonal.xmyL) / 2]);
}

void dpDiagonal_zeroValues(DpDiagonal *diagonal) {
    for (int64_t i = 0; i < diagonal_getWidth(diagonal->diagonal); i++) {
        HDCell *hdCell = diagonal->cells[i];
        for (int64_t p = 0; p < hdCell->numberOfPaths; p++) {
            Path *path = hdCell_getPath(hdCell, p);
            for (int64_t s = 0; s < path->stateNumber; s++) {
                path->cells[s] = LOG_ZERO;
            }
        }
    }
}

void dpDoagonal_setValues(DpDiagonal *diagonal, StateMachine *sM,
                          double (*getStateValue)(StateMachine *, int64_t)) {
    for (int64_t i = 0; i < diagonal_getWidth(diagonal->diagonal); i++) {
        HDCell *hdCell = diagonal->cells[i];
        for (int64_t p = 0; p < hdCell->numberOfPaths; p++) {
            Path *path = hdCell_getPath(hdCell, p);
            double *pathCells = path_getCell(path);
            for (int64_t s = 0; s < path->stateNumber; s++) {
                pathCells[s] = getStateValue(sM, s);
            }
        }
    }
}

void dpDiagonal_initialiseValues(DpDiagonal *diagonal, StateMachine *sM,
                                 double (*getStateValue)(StateMachine *, int64_t)) {
    for (int64_t i = diagonal_getMinXmy(diagonal->diagonal); i <= diagonal_getMaxXmy(diagonal->diagonal); i += 2) {
        HDCell *hdCell = dpDiagonal_getCell(diagonal, i);

        if (hdCell->init == FALSE) {
            st_errAbort("dpDiagonal_initializeValues: got uninitialized cell at %lld\n", i);
        }

        for (int64_t p = 0; p < hdCell->numberOfPaths; p++) {
            Path *path = hdCell_getPath(hdCell, p);
            double *cell = path_getCell(path);
            for (int64_t j = 0; j < path->stateNumber; j++) {
                cell[j] = getStateValue(sM, j);
            }
        }
    }
}

double dpDiagonal_dotProduct(DpDiagonal *diagonal1, DpDiagonal *diagonal2) {
    double totalProbability = LOG_ZERO;
    Diagonal diagonal = diagonal1->diagonal;
    int64_t xmy = diagonal_getMinXmy(diagonal);

    while (xmy <= diagonal_getMaxXmy(diagonal)) {
        HDCell *hdCell1 = dpDiagonal_getCell(diagonal1, xmy);
        HDCell *hdCell2 = dpDiagonal_getCell(diagonal2, xmy);
        totalProbability = logAdd(totalProbability, hdCell_totalProbability(hdCell1, hdCell2));
        xmy += 2;
    }

    return totalProbability;
}

void dpDiagonal_destruct(DpDiagonal *dpDiagonal) {
    int64_t xmy = diagonal_getMinXmy(dpDiagonal->diagonal);
    int64_t maxXmy = diagonal_getMaxXmy(dpDiagonal->diagonal);

    while (xmy < maxXmy) {
        HDCell *hdCell = dpDiagonal_getCell(dpDiagonal, xmy);
        hdCell_destruct(hdCell);
        xmy +=2;
    }
    //for (int64_t i = 0; i < diagonal_getWidth(dpDiagonal->diagonal); i++) {
    //    HDCell *cell = dpDiagonal->cells[i];
    //    hdCell_destruct(cell);
    //}

    free(dpDiagonal->cells);
    free(dpDiagonal);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//DpMatrix
//
//Structure for storing dp-matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////

struct _dpMatrix {
    DpDiagonal **diagonals;
    int64_t diagonalNumber;
    int64_t activeDiagonals;
    int64_t stateNumber;
};

DpMatrix *dpMatrix_construct(int64_t diagonalNumber, int64_t stateNumber) {
    assert(diagonalNumber >= 0);
    DpMatrix *dpMatrix = st_malloc(sizeof(DpMatrix));
    dpMatrix->diagonalNumber = diagonalNumber;
    dpMatrix->diagonals = st_calloc(dpMatrix->diagonalNumber + 1, sizeof(DpDiagonal *));
    dpMatrix->activeDiagonals = 0;
    dpMatrix->stateNumber = stateNumber;
    return dpMatrix;
}

void dpMatrix_destruct(DpMatrix *dpMatrix) {
    assert(dpMatrix->activeDiagonals == 0);
    free(dpMatrix->diagonals);
    free(dpMatrix);
}

DpDiagonal *dpMatrix_getDiagonal(DpMatrix *dpMatrix, int64_t xay) {
    if (xay < 0 || xay > dpMatrix->diagonalNumber) {
        return NULL;
    }
    return dpMatrix->diagonals[xay];
}

int64_t dpMatrix_getActiveDiagonalNumber(DpMatrix *dpMatrix) {
    return dpMatrix->activeDiagonals;
}

DpDiagonal *dpMatrix_createDiagonal(DpMatrix *dpMatrix, Diagonal diagonal, Sequence *sX) {
    if (sX->type != kmer) {
        st_errAbort("dpMatrix_createDiagonal: wrong king of sequence got: %lld\n", sX->type);
    }

    assert(diagonal.xay >= 0);
    assert(diagonal.xay <= dpMatrix->diagonalNumber);
    assert(dpMatrix_getDiagonal(dpMatrix, diagonal.xay) == NULL);
    DpDiagonal *dpDiagonal = dpDiagonal_construct(diagonal, dpMatrix->stateNumber, sX);
    dpMatrix->diagonals[diagonal_getXay(diagonal)] = dpDiagonal;
    dpMatrix->activeDiagonals++;
    return dpDiagonal;
}

void dpMatrix_deleteDiagonal(DpMatrix *dpMatrix, int64_t xay) {
    assert(xay >= 0);
    assert(xay <= dpMatrix->diagonalNumber);
    if (dpMatrix->diagonals[xay] != NULL) {
        dpMatrix->activeDiagonals--;
        assert(dpMatrix->activeDiagonals >= 0);
        dpDiagonal_destruct(dpMatrix->diagonals[xay]);
        dpMatrix->diagonals[xay] = NULL;
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Diagonal DP Calculations
//
//Functions which do forward/backward/posterior calculations
//between diagonal rows of a dp-matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////

static int64_t getXposition(Sequence *sX, int64_t xay, int64_t xmy) {
    int64_t x = diagonal_getXCoordinate(xay, xmy);
    assert(x >= 0 && x <= sX->length);
    return x;
}

static int64_t getYposition(Sequence *sY, int64_t xay, int64_t xmy) {
    int64_t y = diagonal_getYCoordinate(xay, xmy);
    assert(y >= 0 && y <= sY->length);
    return y;
}

static void diagonalCalculation(StateMachine *sM,
                                DpDiagonal *dpDiagonal, DpDiagonal *dpDiagonalM1, DpDiagonal *dpDiagonalM2,
                                Sequence* sX, Sequence* sY,
                                void (*cellCalculation)(StateMachine *, void *, void *, void *, void *,
                                                        void *, void *, void *),
                                void *extraArgs) {
    Diagonal diagonal = dpDiagonal->diagonal;

    // get the smallest x - y coordinate
    int64_t xmy = diagonal_getMinXmy(diagonal);

    // work from smallest to largest
    while (xmy <= diagonal_getMaxXmy(diagonal)) {
        // get the position in the sequence based on the diagonals
        int64_t indexX = getXposition(sX, diagonal_getXay(diagonal), xmy) - 1;
        int64_t indexY = getYposition(sY, diagonal_getXay(diagonal), xmy) - 1;

        // get the element from the sequence at that position. At this point the element is still a void, so
        // it could be anything (base, kmer, event, etc.)
        void* x = sX->get(sX->elements, indexX);
        void* y = sY->get(sY->elements, indexY);

        // do the calculations

        HDCell *current = dpDiagonal_getCell(dpDiagonal, xmy);
        HDCell *lower = dpDiagonalM1 == NULL ? NULL : dpDiagonal_getCell(dpDiagonalM1, xmy - 1);
        HDCell *middle = dpDiagonalM2 == NULL ? NULL : dpDiagonal_getCell(dpDiagonalM2, xmy);
        HDCell *upper = dpDiagonalM1 == NULL ? NULL : dpDiagonal_getCell(dpDiagonalM1, xmy + 1);
        cellCalculation(sM, current, lower, middle, upper, x, y, extraArgs);
        xmy += 2;

    }
}

void diagonalCalculationForward(StateMachine *sM, int64_t xay, DpMatrix *dpMatrix,
                                Sequence* sX, Sequence* sY) {
    diagonalCalculation(sM,
                        dpMatrix_getDiagonal(dpMatrix, xay),
                        dpMatrix_getDiagonal(dpMatrix, xay - 1),
                        dpMatrix_getDiagonal(dpMatrix, xay - 2),
                        sX, sY,
                        cell_calculateForward,
                        NULL);
}

void diagonalCalculationBackward(StateMachine *sM, int64_t xay, DpMatrix *dpMatrix,
                                 Sequence* sX, Sequence* sY) {
    diagonalCalculation(sM,
                        dpMatrix_getDiagonal(dpMatrix, xay),
                        dpMatrix_getDiagonal(dpMatrix, xay - 1),
                        dpMatrix_getDiagonal(dpMatrix, xay - 2),
                        sX, sY,
                        cell_calculateBackward,
                        NULL);
}

double diagonalCalculationTotalProbability(StateMachine *sM, int64_t xay, DpMatrix *forwardDpMatrix,
                                           DpMatrix *backwardDpMatrix, Sequence* sX, Sequence* sY) {
    //Get the forward and backward diagonals
    DpDiagonal *forwardDiagonal = dpMatrix_getDiagonal(forwardDpMatrix, xay);
    DpDiagonal *backDiagonal = dpMatrix_getDiagonal(backwardDpMatrix, xay);
    double totalProbability = dpDiagonal_dotProduct(forwardDiagonal, backDiagonal);

    //Now calculate the contribution of matches through xay.
    forwardDiagonal = dpMatrix_getDiagonal(forwardDpMatrix, xay - 1);
    backDiagonal = dpMatrix_getDiagonal(backwardDpMatrix, xay + 1);
    if (backDiagonal != NULL && forwardDiagonal != NULL) {
        DpDiagonal *matchDiagonal = dpDiagonal_clone(backDiagonal);
        dpDiagonal_zeroValues(matchDiagonal);
        diagonalCalculation(sM, matchDiagonal, NULL, forwardDiagonal, sX, sY, cell_calculateForward, NULL);
        totalProbability = logAdd(totalProbability, dpDiagonal_dotProduct(matchDiagonal, backDiagonal));
        dpDiagonal_destruct(matchDiagonal);
    }
    return totalProbability;
}

void diagonalCalculationPosteriorMatchProbs(StateMachine *sM, int64_t xay, DpMatrix *forwardDpMatrix,
                                            DpMatrix *backwardDpMatrix, Sequence* sX, Sequence* sY,
                                            double totalProbability, PairwiseAlignmentParameters *p, void *extraArgs) {
    assert(p->threshold >= 0.0);
    assert(p->threshold <= 1.0);
    stList *alignedPairs = ((void **) extraArgs)[0];
    DpDiagonal *forwardDiagonal = dpMatrix_getDiagonal(forwardDpMatrix, xay);
    DpDiagonal *backDiagonal = dpMatrix_getDiagonal(backwardDpMatrix, xay);
    Diagonal diagonal = forwardDiagonal->diagonal;
    int64_t xmy = diagonal_getMinXmy(diagonal);

    //Walk over the cells computing the posteriors
    while (xmy <= diagonal_getMaxXmy(diagonal)) {
        int64_t x = diagonal_getXCoordinate(diagonal_getXay(diagonal), xmy);
        int64_t y = diagonal_getYCoordinate(diagonal_getXay(diagonal), xmy);
        if (x > 0 && y > 0) {
            HDCell *cellForward = dpDiagonal_getCell(forwardDiagonal, xmy);
            //st_uglyf("X: %lld Y: %lld cellForward->MatchState: %f ", x, y, cellForward[sM->matchState]);
            HDCell *cellBackward = dpDiagonal_getCell(backDiagonal, xmy);
            //st_uglyf("cellBackward->MatchState: %f ", cellBackward[sM->matchState]);
            for (int64_t fp = 0; fp < cellForward->numberOfPaths; fp++) {
                Path *pathForward = cellForward->paths[fp];
                for (int64_t bp = 0; bp < cellBackward->numberOfPaths; bp++) {
                    Path *pathBackwards = cellBackward->paths[bp];
                    if (stString_eq(pathForward->kmer, pathBackwards->kmer)) {
                        double *forwardCells = path_getCell(pathForward);
                        double *backwardCells = path_getCell(pathBackwards);
                        double posteriorProbability = exp(
                                (forwardCells[sM->matchState] + backwardCells[sM->matchState]) - totalProbability);
                        if (posteriorProbability >= p->threshold) {
                            if (posteriorProbability > 1.0) {
                                posteriorProbability = 1.0;
                            }
                            //st_uglyf("Adding to alignedPairs! posteriorProb: %f, X: %lld (%s), Y: %lld (%f)\n", posteriorProbability, x - 1, sX->get(sX->elements, x-1), y - 1, *(double *)sY->get(sY->elements, y-1));
                            //st_uglyf("Adding to alignedPairs! posteriorProb: %f, X: %lld (%s), Y: %lld (%f)\n", posteriorProbability, x - 1, pathForward->kmer, y - 1, *(double *)sY->get(sY->elements, y-1));
                            //st_uglyf("Adding to alignedPairs! posteriorProb: %f, X: %lld, Y: %lld (%f)\n", posteriorProbability, x - 1, y - 1, *(double *)sY->get(sY->elements, y-1));
                            posteriorProbability = floor(posteriorProbability * PAIR_ALIGNMENT_PROB_1);
                            stList_append(alignedPairs, stIntTuple_construct4((int64_t) posteriorProbability,
                                                                              x - 1, y - 1,
                                                                              (int64_t )pathForward->kmer));
                        }
                    }
                }
            }
            //if (posteriorProbability <= p->threshold) {
            //    //st_uglyf("NOT Adding to alignedPairs! posteriorProb: %f, X: %lld, Y: %lld (%f)\n", posteriorProbability, x - 1, y - 1, *(double *)sY->get(sY->elements, y-1));
            //}
        }
        xmy += 2;
    }
    //st_uglyf("final length for alignedPairs: %lld\n", stList_length(alignedPairs));
}

void diagonalCalculationMultiPosteriorMatchProbs(StateMachine *sM, int64_t xay, DpMatrix *forwardDpMatrix,
                                                 DpMatrix *backwardDpMatrix, Sequence* sX, Sequence* sY,
                                                 double totalProbability, PairwiseAlignmentParameters *p,
                                                 void *extraArgs) {
    assert(p->threshold >= 0.0);
    assert(p->threshold <= 1.0);
    stList *alignedPairs = ((void **) extraArgs)[0];
    DpDiagonal *forwardDiagonal = dpMatrix_getDiagonal(forwardDpMatrix, xay);
    DpDiagonal *backDiagonal = dpMatrix_getDiagonal(backwardDpMatrix, xay);
    Diagonal diagonal = forwardDiagonal->diagonal;
    int64_t xmy = diagonal_getMinXmy(diagonal);

    //Walk over the cells computing the posteriors
    while (xmy <= diagonal_getMaxXmy(diagonal)) {
        int64_t x = diagonal_getXCoordinate(diagonal_getXay(diagonal), xmy);
        int64_t y = diagonal_getYCoordinate(diagonal_getXay(diagonal), xmy);
        if (x > 0 && y > 0) {
            double *cellForward = dpDiagonal_getCell(forwardDiagonal, xmy);
            //st_uglyf("X: %lld Y: %lld cellForward->MatchState: %f ", x, y, cellForward[sM->matchState]);
            double *cellBackward = dpDiagonal_getCell(backDiagonal, xmy);
            //st_uglyf("cellBackward->MatchState: %f ", cellBackward[sM->matchState]);
            for (int64_t s = sM->matchState; s < 6; s++) {
                double posteriorProbability = exp(
                        (cellForward[s] + cellBackward[s]) - totalProbability);
                if (posteriorProbability >= p->threshold) {
                    if (posteriorProbability > 1.0) {
                        posteriorProbability = 1.0;
                    }
                    posteriorProbability = floor(posteriorProbability * PAIR_ALIGNMENT_PROB_1);
                    for (int64_t n = 0; n < s; n++) {
                        //st_uglyf("Adding to alignedPairs! posteriorProb: %f, X: %lld, Y: %lld (%f), state:%lld \n", posteriorProbability/PAIR_ALIGNMENT_PROB_1, (x + n) - 1, y - 1, *(double *)sY->get(sY->elements, y-1), s);
                        stList_append(alignedPairs, stIntTuple_construct3((int64_t) posteriorProbability, (x + n) - 1, y - 1));
                    }
                }
                //if (posteriorProbability <= p->threshold) {
                //    //st_uglyf("NOT adding to alignedPairs! posteriorProb: %f, X: %lld, Y: %lld (%f), s:%lld \n", posteriorProbability, x - 1, y - 1, *(double *)sY->get(sY->elements, y-1), s);
                //}
            }
        }
        xmy += 2;
    }
    //st_uglyf("final length for alignedPairs: %lld\n", stList_length(alignedPairs));
}

void diagonalCalculation_Expectations(StateMachine *sM, int64_t xay,
                                      DpMatrix *forwardDpMatrix, DpMatrix *backwardDpMatrix,
                                      Sequence *sX, Sequence *sY,
                                      double totalProbability,
                                      PairwiseAlignmentParameters *p, void *extraArgs) {
    /*
     * Updates the expectations of the transitions/emissions for the given diagonal.
     */
    Hmm *hmmExpectations = extraArgs;
    void *extraArgs2[2] = { &totalProbability, hmmExpectations }; // this is where you pack in totalprob

    // update likelihood
    hmmExpectations->likelihood += totalProbability;

    // We do this once per diagonal, which is a hack, rather than for the
    // whole matrix. The correction factor is approximately 1/number of
    // diagonals.
    diagonalCalculation(sM,
                        dpMatrix_getDiagonal(backwardDpMatrix, xay),
                        dpMatrix_getDiagonal(forwardDpMatrix, xay - 1),
                        dpMatrix_getDiagonal(forwardDpMatrix, xay - 2),
                        sX, sY, cell_calculateUpdateExpectation, extraArgs2);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Banded alignment routine to calculate posterior match probs
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void getPosteriorProbsWithBanding(StateMachine *sM,
                                  stList *anchorPairs,
                                  Sequence *sX, Sequence *sY,
                                  PairwiseAlignmentParameters *p,
                                  bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd,
                                  void (*diagonalPosteriorProbFn)(StateMachine *, int64_t, DpMatrix *, DpMatrix *,
                                                                  Sequence*, Sequence*,
                                                                  double, PairwiseAlignmentParameters *, void *),
                                  void *extraArgs) {
    //Prerequisites
    assert(p->traceBackDiagonals >= 1);
    assert(p->diagonalExpansion >= 0);
    assert(p->diagonalExpansion % 2 == 0);
    assert(p->minDiagsBetweenTraceBack >= 2);
    assert(p->traceBackDiagonals + 1 < p->minDiagsBetweenTraceBack);

    int64_t diagonalNumber = sX->length + sY->length;
    if (diagonalNumber == 0) { //Deal with trivial case
        return;
    }

    //Primitives for the forward matrix recursion
    Band *band = band_construct(anchorPairs, sX->length, sY->length, p->diagonalExpansion);

    BandIterator *forwardBandIterator = bandIterator_construct(band);
    DpMatrix *forwardDpMatrix = dpMatrix_construct(diagonalNumber, sM->stateNumber);
    //Initialise forward matrix.
    dpDiagonal_initialiseValues(dpMatrix_createDiagonal(forwardDpMatrix, bandIterator_getNext(forwardBandIterator), sX), // added sX here
                                sM, alignmentHasRaggedLeftEnd ? sM->raggedStartStateProb : sM->startStateProb);

    //Backward matrix.
    DpMatrix *backwardDpMatrix = dpMatrix_construct(diagonalNumber, sM->stateNumber);

    int64_t tracedBackTo = 0;
    int64_t totalPosteriorCalculations = 0;

    while (1) { //Loop that moves through the matrix forward

        Diagonal diagonal = bandIterator_getNext(forwardBandIterator);

        //Forward calculation
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(forwardDpMatrix, diagonal, sX)); // added sX here
        diagonalCalculationForward(sM, diagonal_getXay(diagonal), forwardDpMatrix, sX, sY);

        //Condition true at the end of the matrix
        bool atEnd = diagonal_getXay(diagonal) == diagonalNumber;
        //Condition true when we want to do an intermediate traceback.
        bool tracebackPoint = diagonal_getXay(diagonal) >= tracedBackTo + p->minDiagsBetweenTraceBack
                              && diagonal_getWidth(diagonal) <= p->diagonalExpansion * 2 + 1;

        //Traceback
        if (atEnd || tracebackPoint) {
            //Initialise the last row (until now) of the backward matrix to represent an end point
            dpDiagonal_initialiseValues(dpMatrix_createDiagonal(backwardDpMatrix, diagonal, sX), sM, // added sX here
                                        (atEnd && alignmentHasRaggedRightEnd) ? sM->raggedEndStateProb : sM->endStateProb);
            //This is a diagonal between the place we trace back to and where we trace back from
            if (diagonal_getXay(diagonal) > tracedBackTo + 1) {
                DpDiagonal *j = dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal) - 1);
                assert(j != NULL);
                dpDiagonal_zeroValues(dpMatrix_createDiagonal(backwardDpMatrix, j->diagonal, sX)); // added sX here
            }

            //Do walk back
            BandIterator *backwardBandIterator = bandIterator_clone(forwardBandIterator);
            Diagonal diagonal2 = bandIterator_getPrevious(backwardBandIterator);
            assert(diagonal_getXay(diagonal2) == diagonal_getXay(diagonal));
            int64_t tracedBackFrom = diagonal_getXay(diagonal) - (atEnd ? 0 : p->traceBackDiagonals + 1);
            double totalProbability = LOG_ZERO;
            int64_t totalPosteriorCalculationsThisTraceback = 0;
            while (diagonal_getXay(diagonal2) > tracedBackTo) {
                //Create the earlier diagonal
                if (diagonal_getXay(diagonal2) > tracedBackTo + 2) {
                    DpDiagonal *j = dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2) - 2);
                    assert(j != NULL);
                    dpDiagonal_zeroValues(dpMatrix_createDiagonal(backwardDpMatrix, j->diagonal, sX)); // added sX here
                }
                if (diagonal_getXay(diagonal2) > tracedBackTo + 1) {
                    diagonalCalculationBackward(sM, diagonal_getXay(diagonal2), backwardDpMatrix, sX, sY);
                }
                if (diagonal_getXay(diagonal2) <= tracedBackFrom) {
                    assert(dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2)) != NULL);
                    assert(dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2)-1) != NULL);
                    assert(dpMatrix_getDiagonal(backwardDpMatrix, diagonal_getXay(diagonal2)) != NULL);
                    if (diagonal_getXay(diagonal2) != diagonalNumber) {
                        assert(dpMatrix_getDiagonal(backwardDpMatrix, diagonal_getXay(diagonal2)+1) != NULL);
                    }
                    if (totalPosteriorCalculationsThisTraceback++ % 10 == 0) {
                        double newTotalProbability = diagonalCalculationTotalProbability(
                                sM, diagonal_getXay(diagonal2),
                                forwardDpMatrix, backwardDpMatrix, sX, sY
                        );
                        if (totalPosteriorCalculationsThisTraceback != 1) {
                            assert(totalProbability + 1.0 > newTotalProbability);
                            assert(newTotalProbability + 1.0 > newTotalProbability);
                        }
                        totalProbability = newTotalProbability;
                    }
                    diagonalPosteriorProbFn(sM, diagonal_getXay(diagonal2),
                                            forwardDpMatrix, backwardDpMatrix,
                                            sX, sY,
                                            totalProbability, p, extraArgs);
                    //Delete forward diagonal after last access in posterior calculation
                    if (diagonal_getXay(diagonal2) < tracedBackFrom || atEnd) {
                        dpMatrix_deleteDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2));
                    }
                }
                //Delete backward diagonal after last access in backward calculation
                if (diagonal_getXay(diagonal2) + 1 <= diagonalNumber) {
                    dpMatrix_deleteDiagonal(backwardDpMatrix, diagonal_getXay(diagonal2) + 1);
                }
                diagonal2 = bandIterator_getPrevious(backwardBandIterator);
            }
            tracedBackTo = tracedBackFrom;
            bandIterator_destruct(backwardBandIterator);
            dpMatrix_deleteDiagonal(backwardDpMatrix, diagonal_getXay(diagonal2) + 1);
            dpMatrix_deleteDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2));
            //Check memory state.
            assert(dpMatrix_getActiveDiagonalNumber(backwardDpMatrix) == 0);
            totalPosteriorCalculations += totalPosteriorCalculationsThisTraceback;
            if (!atEnd) {
                assert(dpMatrix_getActiveDiagonalNumber(forwardDpMatrix) == p->traceBackDiagonals + 2);
            }
        }
        if (atEnd) {
            break;
        }
    }
    assert(totalPosteriorCalculations == diagonalNumber);
    assert(tracedBackTo == diagonalNumber);
    assert(dpMatrix_getActiveDiagonalNumber(backwardDpMatrix) == 0);
    assert(dpMatrix_getActiveDiagonalNumber(forwardDpMatrix) == 0);
    //Cleanup
    dpMatrix_destruct(forwardDpMatrix);
    dpMatrix_destruct(backwardDpMatrix);
    bandIterator_destruct(forwardBandIterator);
    band_destruct(band);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Blast anchoring functions
//Use lastz to get sets of anchors
/////////////////////////////////////////////////////////////////////////////////////////////////////////

int sortByXPlusYCoordinate(const void *i, const void *j) {
    int64_t k = stIntTuple_get((stIntTuple *) i, 0) + stIntTuple_get((stIntTuple *) i, 1);
    int64_t l = stIntTuple_get((stIntTuple *) j, 0) + stIntTuple_get((stIntTuple *) j, 1);
    return k > l ? 1 : (k < l ? -1 : 0);
}

int sortByXPlusYCoordinate2(const void *i, const void *j) {
    int64_t k = stIntTuple_get((stIntTuple *)i, 1) + stIntTuple_get((stIntTuple *)i, 2);
    int64_t l = stIntTuple_get((stIntTuple *)j, 1) + stIntTuple_get((stIntTuple *)j, 2);
    return k > l ? 1 : (k < l ? -1 : 0);
}

static char *makeUpperCase(const char *s, int64_t l) {
    char *s2 = stString_copy(s);
    for (int64_t i = 0; i < l; i++) {
        s2[i] = toupper(s[i]);
    }
    return s2;
}

static void writeSequenceToFile(char *file, const char *name, const char *sequence) {
    FILE *fileHandle = fopen(file, "w");
    fastaWrite((char *) sequence, (char *) name, fileHandle);
    fclose(fileHandle);
}

stList *convertPairwiseForwardStrandAlignmentToAnchorPairs(struct PairwiseAlignment *pA, int64_t trim) {
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct); //the list to put the output in
    int64_t j = pA->start1;
    int64_t k = pA->start2;
    assert(pA->strand1);
    assert(pA->strand2);
    for (int64_t i = 0; i < pA->operationList->length; i++) {
        struct AlignmentOperation *op = pA->operationList->list[i];
        if (op->opType == PAIRWISE_MATCH) {
            for (int64_t l = trim; l < op->length - trim; l++) {
                stList_append(alignedPairs, stIntTuple_construct2(j + l, k + l));
            }
        }
        if (op->opType != PAIRWISE_INDEL_Y) {
            j += op->length;
        }
        if (op->opType != PAIRWISE_INDEL_X) {
            k += op->length;
        }
    }

    assert(j == pA->end1);
    assert(k == pA->end2);
    return alignedPairs;
}

stList *getBlastPairs(const char *sX, const char *sY, int64_t trim, bool repeatMask) {
    /*
     * Uses lastz to compute a bunch of monotonically increasing pairs such that for any pair of consecutive
     * pairs in the list (x1, y1) (x2, y2) in the set of aligned pairs x1 appears before x2 in X and y1
     * appears before y2 in Y.
     */
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    int64_t lX = strlen(sX);
    int64_t lY = strlen(sY);

    if (lX == 0 || lY == 0) {
        return alignedPairs;
    }

    if (!repeatMask) {
        sX = makeUpperCase(sX, lX);
        sY = makeUpperCase(sY, lY);
    }

    //Write one sequence to file..
    char *tempFile1 = getTempFile();
    char *tempFile2 = NULL;

    writeSequenceToFile(tempFile1, "a", sX);

    char *command;

    if (lY > 1000) {
        tempFile2 = getTempFile();
        writeSequenceToFile(tempFile2, "b", sY);
        command =
                stString_print(
                        //"./cPecanLastz --hspthresh=800 --chain --strand=plus --gapped --format=cigar --ambiguous=iupac,100,100 %s %s",
                        "./cPecanLastz --hspthresh=1800 --chain --strand=plus --gapped --format=cigar --gap=100,100 --ambiguous=iupac,100,100 %s %s",
                        tempFile1, tempFile2);
    } else {
        command =
                stString_print(
                        //"echo '>b\n%s\n' | ./cPecanLastz --hspthresh=800 --chain --strand=plus --gapped --format=cigar --ambiguous=iupac,100,100 %s",
                        "echo '>b\n%s\n' | ./cPecanLastz --hspthresh=1800 --chain --strand=plus --gapped --gap=100,100 --format=cigar --ambiguous=iupac,100,100 %s",
                        sY, tempFile1);
    }
    FILE *fileHandle = popen(command, "r");
    if (fileHandle == NULL) {
        st_errAbort("Problems with lastz pipe");
    }
    //Read from stream
    struct PairwiseAlignment *pA;
    while ((pA = cigarRead(fileHandle)) != NULL) {
        assert(strcmp(pA->contig1, "a") == 0);
        assert(strcmp(pA->contig2, "b") == 0);
        stList *alignedPairsForCigar = convertPairwiseForwardStrandAlignmentToAnchorPairs(pA, trim);
        stList_appendAll(alignedPairs, alignedPairsForCigar);
        stList_setDestructor(alignedPairsForCigar, NULL);
        stList_destruct(alignedPairsForCigar);
        destructPairwiseAlignment(pA);
    }
    int64_t status = pclose(fileHandle);
    if (status != 0) {
        st_errAbort("pclose failed when getting rid of lastz pipe with value %" PRIi64 " and command %s", status,
                command);
    }
    free(command);

    stList_sort(alignedPairs, sortByXPlusYCoordinate); //Ensure the coordinates are increasing
    //Remove old files
    st_system("rm %s", tempFile1);
    free(tempFile1);
    if (tempFile2 != NULL) {
        st_system("rm %s", tempFile2);
        free(tempFile2);
    }

    if (!repeatMask) {
        free((char *) sX);
        free((char *) sY);
    }

    return alignedPairs;
}

static void convertBlastPairs(stList *alignedPairs2, int64_t offsetX, int64_t offsetY) {
    /*
     * Convert the coordinates of the computed pairs.
     */
    for (int64_t k = 0; k < stList_length(alignedPairs2); k++) {
        stIntTuple *i = stList_get(alignedPairs2, k);
        assert(stIntTuple_length(i) == 2);
        stList_set(alignedPairs2, k,
                stIntTuple_construct2(stIntTuple_get(i, 0) + offsetX, stIntTuple_get(i, 1) + offsetY));
        stIntTuple_destruct(i);
    }
}

stList *filterToRemoveOverlap(stList *sortedOverlappingPairs) {
    stList *nonOverlappingPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    //Traverse backwards
    stSortedSet *set = stSortedSet_construct3((int (*)(const void *, const void *)) stIntTuple_cmpFn, NULL);
    int64_t pX = INT64_MAX, pY = INT64_MAX;
    for (int64_t i = stList_length(sortedOverlappingPairs) - 1; i >= 0; i--) {
        stIntTuple *pair = stList_get(sortedOverlappingPairs, i);
        int64_t x = stIntTuple_get(pair, 0);
        int64_t y = stIntTuple_get(pair, 1);
        if (x < pX && y < pY) {
            stSortedSet_insert(set, pair);
        }
        pX = x < pX ? x : pX;
        pY = y < pY ? y : pY;
    }

    //Traverse forwards to final set of pairs
    pX = INT64_MIN;
    pY = INT64_MIN;
    int64_t pY2 = INT64_MIN;
    for (int64_t i = 0; i < stList_length(sortedOverlappingPairs); i++) {
        stIntTuple *pair = stList_get(sortedOverlappingPairs, i);
        int64_t x = stIntTuple_get(pair, 0);
        int64_t y = stIntTuple_get(pair, 1);
        if (x > pX && y > pY && stSortedSet_search(set, pair) != NULL) {
            stList_append(nonOverlappingPairs, stIntTuple_construct2(x, y));
        }
        //Check things are sorted in the input
        assert(x >= pX);
        if (x == pX) {
            assert(y >= pY2);
        }
        pY2 = y;
        pX = x > pX ? x : pX;
        pY = y > pY ? y : pY;
    }
    stSortedSet_destruct(set);

    return nonOverlappingPairs;
}

static void getBlastPairsForPairwiseAlignmentParametersP(
                        const char *sX, const char *sY, int64_t pX, int64_t pY,
                        int64_t x, int64_t y, PairwiseAlignmentParameters *p,
                        stList *combinedAnchorPairs) {

    int64_t lX2 = x - pX;
    assert(lX2 >= 0);
    int64_t lY2 = y - pY;
    assert(lY2 >= 0);
    int64_t matrixSize = (int64_t) lX2 * lY2;
    if (matrixSize > p->repeatMaskMatrixBiggerThanThis) {
        char *sX2 = stString_getSubString(sX, pX, lX2);
        char *sY2 = stString_getSubString(sY, pY, lY2);
        stList *unfilteredBottomLevelAnchorPairs = getBlastPairs(sX2, sY2, p->constraintDiagonalTrim, 0);
        stList_sort(unfilteredBottomLevelAnchorPairs, (int (*)(const void *, const void *)) stIntTuple_cmpFn);
        stList *bottomLevelAnchorPairs = filterToRemoveOverlap(unfilteredBottomLevelAnchorPairs);
        st_logDebug("Got %" PRIi64 " bottom level anchor pairs, which reduced to %" PRIi64 " after filtering \n",
                stList_length(unfilteredBottomLevelAnchorPairs), stList_length(bottomLevelAnchorPairs));
        stList_destruct(unfilteredBottomLevelAnchorPairs);
        convertBlastPairs(bottomLevelAnchorPairs, pX, pY);
        free(sX2);
        free(sY2);
        stList_appendAll(combinedAnchorPairs, bottomLevelAnchorPairs);
        stList_setDestructor(bottomLevelAnchorPairs, NULL);
        stList_destruct(bottomLevelAnchorPairs);
    }
}

stList *getBlastPairsForPairwiseAlignmentParameters(void *sX, void *sY, PairwiseAlignmentParameters *p) {

    // cast to char arrays for lastz
    char *cX = (char *) sX;
    char *cY = (char *) sY;
    int64_t lX = strlen(cX);
    int64_t lY = strlen(cY);

    if ((int64_t) lX * lY <= p->anchorMatrixBiggerThanThis) {
        return stList_construct();
    }
    // anchorPairs
    // Get anchors
    stList *unfilteredTopLevelAnchorPairs = getBlastPairs(cX, cY, p->constraintDiagonalTrim, 1);
    // sort them
    stList_sort(unfilteredTopLevelAnchorPairs, (int (*)(const void *, const void *)) stIntTuple_cmpFn);
    // filter
    stList *topLevelAnchorPairs = filterToRemoveOverlap(unfilteredTopLevelAnchorPairs);
    st_logDebug("Got %" PRIi64 " top level anchor pairs, which reduced to %" PRIi64 " after filtering \n",
            stList_length(unfilteredTopLevelAnchorPairs), stList_length(topLevelAnchorPairs));
    // intermediate cleanup
    stList_destruct(unfilteredTopLevelAnchorPairs);

    // go though topLevelAnchorPairs and combine them based on certain conditions
    int64_t pX = 0;
    int64_t pY = 0;
    stList *combinedAnchorPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(topLevelAnchorPairs); i++) {
        stIntTuple *anchorPair = stList_get(topLevelAnchorPairs, i);
        int64_t x = stIntTuple_get(anchorPair, 0);
        int64_t y = stIntTuple_get(anchorPair, 1);
        // make sure x and y are within the length of the sequence
        assert(x >= 0 && x < lX);
        assert(y >= 0 && y < lY);
        // make sure x and y are 'in front of' the last pair
        assert(x >= pX);
        assert(y >= pY);
        // see if we want to split the matrix into two
        getBlastPairsForPairwiseAlignmentParametersP(cX, cY, pX, pY, x, y, p, combinedAnchorPairs);
        // finally append
        stList_append(combinedAnchorPairs, anchorPair);
        // increment for next iteration
        pX = x + 1;
        pY = y + 1;
    }
    // one final check
    getBlastPairsForPairwiseAlignmentParametersP(cX, cY, pX, pY, lX, lY, p, combinedAnchorPairs);
    stList_setDestructor(topLevelAnchorPairs, NULL);
    stList_destruct(topLevelAnchorPairs);
    st_logDebug("Got %" PRIi64 " combined anchor pairs\n", stList_length(combinedAnchorPairs));
    return combinedAnchorPairs;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Split large gap functions
//Functions to split up alignment around gaps in the anchors that are too large.
/////////////////////////////////////////////////////////////////////////////////////////////////////////


static bool getSplitPointsP(int64_t *x1, int64_t *y1, int64_t x2, int64_t y2, int64_t x3, int64_t y3,
        stList *splitPoints, int64_t splitMatrixBiggerThanThis, bool skipBlock) {
    /*
     * x2/y2 are the previous anchor point, x3/y3 are the next anchor point. Gaps greater than (x3-x2)*(y3-y2) are split up.
     */
    int64_t lX2 = x3 - x2;
    int64_t lY2 = y3 - y2;
    int64_t matrixSize = lX2 * lY2;
    if (matrixSize > splitMatrixBiggerThanThis) {
        st_logDebug("Split point found at x1: %" PRIi64 " x2: %" PRIi64 " y1: %" PRIi64 " y2: %" PRIi64 "\n", x2, x3,
                y2, y3);
        int64_t maxSequenceLength = sqrt(splitMatrixBiggerThanThis);
        int64_t hX = lX2 / 2 > maxSequenceLength ? maxSequenceLength : lX2 / 2;
        int64_t hY = lY2 / 2 > maxSequenceLength ? maxSequenceLength : lY2 / 2;
        if(!skipBlock) {
            stList_append(splitPoints, stIntTuple_construct4(*x1, *y1, x2 + hX, y2 + hY));
        }
        *x1 = x3 - hX;
        *y1 = y3 - hY;
        return 1;
    }
    return 0;
}

stList *getSplitPoints(stList *anchorPairs, int64_t lX, int64_t lY, int64_t splitMatrixBiggerThanThis,
                       bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {
    int64_t x1 = 0, y1 = 0, x2 = 0, y2 = 0;
    assert(lX >= 0);
    assert(lY >= 0);
    stList *splitPoints = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(anchorPairs); i++) {
        stIntTuple *anchorPair = stList_get(anchorPairs, i);
        int64_t x3 = stIntTuple_get(anchorPair, 0), y3 = stIntTuple_get(anchorPair, 1);
        getSplitPointsP(&x1, &y1, x2, y2, x3, y3, splitPoints, splitMatrixBiggerThanThis, alignmentHasRaggedLeftEnd && i == 0);
        assert(x3 >= x2);
        assert(y3 >= y2);
        assert(x3 < lX);
        assert(y3 < lY);
        x2 = x3 + 1;
        y2 = y3 + 1;
    }
    if(!getSplitPointsP(&x1, &y1, x2, y2, lX, lY, splitPoints, splitMatrixBiggerThanThis,
            alignmentHasRaggedLeftEnd && stList_length(anchorPairs) == 0) || !alignmentHasRaggedRightEnd) {
        stList_append(splitPoints, stIntTuple_construct4(x1, y1, lX, lY));
    }

    if (stList_length(splitPoints) > 1) {
        st_logDebug("For sequences of length %" PRIi64 " and %" PRIi64 " we got %" PRIi64 " splits\n", lX, lY,
                stList_length(splitPoints));
    }
    return splitPoints;
}

static void convertAlignedPairs(stList *alignedPairs2, int64_t offsetX, int64_t offsetY) {
    /*
     * Convert the coordinates of the computed pairs.
     */
    for (int64_t k = 0; k < stList_length(alignedPairs2); k++) {
        stIntTuple *i = stList_get(alignedPairs2, k);
        assert(stIntTuple_length(i) == 4);
        stList_set(alignedPairs2, k,
                stIntTuple_construct4(stIntTuple_get(i, 0), stIntTuple_get(i, 1) + offsetX, //TODO make this work with kmers
                        stIntTuple_get(i, 2) + offsetY, stIntTuple_get(i, 3)));
        stIntTuple_destruct(i);
    }
}

void getPosteriorProbsWithBandingSplittingAlignmentsByLargeGaps(
        StateMachine *sM, stList *anchorPairs, Sequence *SsX, Sequence *SsY,
        PairwiseAlignmentParameters *p,
        bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd,
        void (*diagonalPosteriorProbFn)(StateMachine *, int64_t, DpMatrix *,
                                        DpMatrix *, Sequence*, Sequence*, double,
                                        PairwiseAlignmentParameters *, void *),
        void (*coordinateCorrectionFn)(), void *extraArgs) {
    // you are going to cut the sequences into subSequences anyways, so not having the correct
    // number of elements in length, ie having it reflect the number of nucleotides might be ok?
    int64_t lX = SsX->length; // so here you want the total number of elements
    int64_t lY = SsY->length;

    stList *splitPoints = getSplitPoints(anchorPairs, lX, lY,
                                         p->splitMatrixBiggerThanThis,
                                         alignmentHasRaggedLeftEnd,
                                         alignmentHasRaggedRightEnd);
    int64_t j = 0;

    //Now to the actual alignments
    for (int64_t i = 0; i < stList_length(splitPoints); i++) {
        stIntTuple *subRegion = stList_get(splitPoints, i);
        int64_t x1 = stIntTuple_get(subRegion, 0);
        int64_t y1 = stIntTuple_get(subRegion, 1);
        int64_t x2 = stIntTuple_get(subRegion, 2);
        int64_t y2 = stIntTuple_get(subRegion, 3);

        Sequence *sX3 = SsX->sliceFcn(SsX, x1, x2 - x1);
        Sequence *sY3 = SsY->sliceFcn(SsY, y1, y2 - y1);

        //List of anchor pairs
        stList *subListOfAnchorPoints = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

        while (j < stList_length(anchorPairs)) {
            stIntTuple *anchorPair = stList_get(anchorPairs, j);
            int64_t x = stIntTuple_get(anchorPair, 0);
            int64_t y = stIntTuple_get(anchorPair, 1);

            assert(x + y >= x1 + y1);
            if (x + y >= x2 + y2) {
                break;
            }
            assert(x >= x1 && x < x2);
            assert(y >= y1 && y < y2);
            stList_append(subListOfAnchorPoints, stIntTuple_construct2(x - x1, y - y1));
            j++;
        }

        //Make the alignments

        getPosteriorProbsWithBanding(sM, subListOfAnchorPoints, sX3, sY3, p,
                                     (alignmentHasRaggedLeftEnd || i > 0),
                                     (alignmentHasRaggedRightEnd || i < stList_length(splitPoints) - 1),
                                     diagonalPosteriorProbFn, extraArgs);

        if (coordinateCorrectionFn != NULL) {
            coordinateCorrectionFn(x1, y1, extraArgs);
        }

        //Clean up
        stList_destruct(subListOfAnchorPoints);
        sequence_sequenceDestroy(sX3);
        sequence_sequenceDestroy(sY3);
    }
    assert(j == stList_length(anchorPairs));
    stList_destruct(splitPoints);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Core public functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////

PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters_construct() {
    PairwiseAlignmentParameters *p = st_malloc(sizeof(PairwiseAlignmentParameters));
    p->threshold = 0.01;
    p->minDiagsBetweenTraceBack = 1000;
    p->traceBackDiagonals = 40;
    p->diagonalExpansion = 20;
    p->constraintDiagonalTrim = 14;
    p->anchorMatrixBiggerThanThis = 500 * 500;
    p->repeatMaskMatrixBiggerThanThis = 500 * 500;
    p->splitMatrixBiggerThanThis = (int64_t) 3000 * 3000;
    p->alignAmbiguityCharacters = 0;
    p->gapGamma = 0.5;
    return p;
}

void pairwiseAlignmentBandingParameters_destruct(PairwiseAlignmentParameters *p) {
    free(p);
}

static void alignedPairCoordinateCorrectionFn(int64_t offsetX, int64_t offsetY, void *extraArgs) {
    stList *subListOfAlignedPairs = ((void **) extraArgs)[0];
    stList *alignedPairs = ((void **) extraArgs)[1];
    convertAlignedPairs(subListOfAlignedPairs, offsetX, offsetY); //Shift back the aligned pairs to the appropriate coordinates
    while (stList_length(subListOfAlignedPairs) > 0) {
        stList_append(alignedPairs, stList_pop(subListOfAlignedPairs));
    }
}

stList *getAlignedPairsUsingAnchors(StateMachine *sM,
                                    Sequence *SsX, Sequence *SsY,
                                    stList *anchorPairs,
                                    PairwiseAlignmentParameters *p,
                                    void (*diagonalPosteriorProbFn)(StateMachine *, int64_t, DpMatrix *,
                                                                    DpMatrix *, Sequence *, Sequence *, double,
                                                                    PairwiseAlignmentParameters *, void *),
                                    bool alignmentHasRaggedLeftEnd,
                                    bool alignmentHasRaggedRightEnd) {

    //This list of pairs to be returned. Not in any order, but points must be unique
    stList *subListOfAlignedPairs = stList_construct();
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    void *extraArgs[2] = { subListOfAlignedPairs, alignedPairs };

    getPosteriorProbsWithBandingSplittingAlignmentsByLargeGaps(sM, anchorPairs,
                                                               SsX, SsY,
                                                               p,
                                                               alignmentHasRaggedLeftEnd,
                                                               alignmentHasRaggedRightEnd,
                                                               diagonalPosteriorProbFn,
                                                               alignedPairCoordinateCorrectionFn,
                                                               extraArgs);

    assert(stList_length(subListOfAlignedPairs) == 0);
    stList_destruct(subListOfAlignedPairs);

    return alignedPairs;
}

stList *getAlignedPairs(StateMachine *sM, void *cX, void *cY, int64_t lX, int64_t lY,
                        PairwiseAlignmentParameters *p,
                        void *(*getXFcn)(void *, int64_t),
                        void *(*getYFcn)(void *, int64_t),
                        stList *(*getAnchorPairFcn)(void *, void *, PairwiseAlignmentParameters *),
                        bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {

    //stList *anchorPairs = getBlastPairsForPairwiseAlignmentParameters(cX, cY, p);
    stList *anchorPairs = getAnchorPairFcn(cX, cY, p);

    //Sequence *SsX = sequence_construct(lX, cX, getXFcn);
    //Sequence *SsY = sequence_construct(lY, cY, getYFcn);
    Sequence *SsX = sequence_construct2(lX, cX, getXFcn, sequence_sliceNucleotideSequence2, nucleotide);
    Sequence *SsY = sequence_construct2(lY, cY, getYFcn, sequence_sliceNucleotideSequence2, nucleotide);

    stList *alignedPairs = getAlignedPairsUsingAnchors(sM, SsX, SsY,
                                                       anchorPairs, p,
                                                       diagonalCalculationPosteriorMatchProbs,
                                                       alignmentHasRaggedLeftEnd,
                                                       alignmentHasRaggedRightEnd);
    sequence_sequenceDestroy(SsX);
    sequence_sequenceDestroy(SsY);
    stList_destruct(anchorPairs);
    return alignedPairs;
}

stList *getAlignedPairsWithoutBanding(StateMachine *sM, void *cX, void *cY, int64_t lX, int64_t lY,
                                      PairwiseAlignmentParameters *p,
                                      void *(*getXFcn)(void *, int64_t),
                                      void *(*getYFcn)(void *, int64_t),
                                      void (*diagonalPosteriorProbFn)(StateMachine *, int64_t, DpMatrix *,
                                                                      DpMatrix *, Sequence *, Sequence *, double,
                                                                      PairwiseAlignmentParameters *, void *),
                                      bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {
    // make sequence objects
    Sequence *ScX = sequence_construct(lX, cX, getXFcn, kmer);
    if (sM->type == echelon) {
        sequence_padSequence(ScX);
    }
    Sequence *ScY = sequence_construct(lY, cY, getYFcn, event);

    // make matrices and bands
    int64_t diagonalNumber = ScX->length + ScY->length;
    DpMatrix *forwardDpMatrix = dpMatrix_construct(diagonalNumber, sM->stateNumber);
    DpMatrix *backwardDpMatrix = dpMatrix_construct(diagonalNumber, sM->stateNumber);
    stList *emptyList = stList_construct(); // place holder for band_construct
    Band *band = band_construct(emptyList, ScX->length, ScY->length, 2); // why 2?
    BandIterator *bandIt = bandIterator_construct(band);

    for (int64_t i = 0; i <= diagonalNumber; i++) {
        Diagonal d = bandIterator_getNext(bandIt);
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(backwardDpMatrix, d, ScX)); // added sX here
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(forwardDpMatrix, d, ScX)); // added sX here
    }

    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(forwardDpMatrix, 0), sM,
                                alignmentHasRaggedLeftEnd ? sM->raggedStartStateProb : sM->startStateProb);
    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(backwardDpMatrix, diagonalNumber), sM,
                                alignmentHasRaggedRightEnd ? sM->raggedEndStateProb : sM->endStateProb);

    // perform forward algorithm
    for (int64_t i = 0; i <= diagonalNumber; i++) {
        diagonalCalculationForward(sM, i, forwardDpMatrix, ScX, ScY);
    }
    // perform backward algorithm
    for (int64_t i = diagonalNumber; i > 0; i--) {
        diagonalCalculationBackward(sM, i, backwardDpMatrix, ScX, ScY);
    }
    // calculate total probability, make a place for the aligned pairs to go then run
    // diagoinalCalculationPosteriorMatchProbs
    double totalProbability = diagonalCalculationTotalProbability(sM, diagonalNumber, forwardDpMatrix,
                                                                  backwardDpMatrix, ScX, ScY);
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    void *extraArgs[1] = { alignedPairs };
    for (int64_t i = 0; i <= diagonalNumber; i++) {
        diagonalPosteriorProbFn(sM, i, forwardDpMatrix, backwardDpMatrix, ScX, ScY, totalProbability, p, extraArgs);
    }

    // cleanup
    sequence_sequenceDestroy(ScX);
    sequence_sequenceDestroy(ScY);

    return alignedPairs;
}

void getExpectationsUsingAnchors(StateMachine *sM, Hmm *hmmExpectations,
                                 Sequence *SsX, Sequence *SsY,
                                 stList *anchorPairs,
                                 PairwiseAlignmentParameters *p,
                                 void (*diagonalCalcExpectationFcn)(StateMachine *sM, int64_t xay,
                                                                    DpMatrix *forwardDpMatrix,
                                                                    DpMatrix *backwardDpMatrix,
                                                                    Sequence* sX, Sequence* sY,
                                                                    double totalProbability,
                                                                    PairwiseAlignmentParameters *p,
                                                                    void *extraArgs),
                                 bool alignmentHasRaggedLeftEnd,
                                 bool alignmentHasRaggedRightEnd) {
    getPosteriorProbsWithBandingSplittingAlignmentsByLargeGaps(sM, anchorPairs,
                                                               SsX, SsY,
                                                               p,
                                                               alignmentHasRaggedLeftEnd,
                                                               alignmentHasRaggedRightEnd,
                                                               diagonalCalcExpectationFcn,
                                                               NULL, hmmExpectations);
}

void getExpectations(StateMachine *sM, Hmm *hmmExpectations,
                     void *sX, void *sY, int64_t lX, int64_t lY, PairwiseAlignmentParameters *p,
                     void *(getFcn)(void *, int64_t),
                     stList *(*getAnchorPairFcn)(void *, void *, PairwiseAlignmentParameters *),
                     bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {
    // get anchors
    stList *anchorPairs = getAnchorPairFcn(sX, sY, p);
    // make Sequence objects
    Sequence *SsX = sequence_construct2(lX, sX, getFcn, sequence_sliceNucleotideSequence2, nucleotide);
    Sequence *SsY = sequence_construct2(lY, sY, getFcn, sequence_sliceNucleotideSequence2, nucleotide);

    getExpectationsUsingAnchors(sM, hmmExpectations, SsX, SsY,
                                anchorPairs, p,
                                diagonalCalculation_Expectations,
                                alignmentHasRaggedLeftEnd,
                                alignmentHasRaggedRightEnd);

    sequence_sequenceDestroy(SsX);
    sequence_sequenceDestroy(SsY);
    stList_destruct(anchorPairs);
}

/*
 * Functions for adjusting weights to account for probability of alignment to a gap. These were unchanged.
 */
int64_t *getIndelProbabilities(stList *alignedPairs, int64_t seqLength, bool xIfTrueElseY) {
    int64_t *indelProbs = st_malloc(seqLength * sizeof(int64_t));
    for(int64_t i=0; i<seqLength; i++) {
        indelProbs[i] = PAIR_ALIGNMENT_PROB_1;
    }
    for(int64_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *j = stList_get(alignedPairs, i);
        indelProbs[stIntTuple_get(j, xIfTrueElseY ? 1 : 2)] -= stIntTuple_get(j, 0);
    }
    for(int64_t i=0; i<seqLength; i++) {
        if(indelProbs[i] < 0) {
            indelProbs[i] = 0;
        }
    }
    return indelProbs;
}

stList *reweightAlignedPairs(stList *alignedPairs,
                             int64_t *indelProbsX, int64_t *indelProbsY,
                             double gapGamma) {
    stList *reweightedAlignedPairs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for(int64_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(aPair, 1);
        int64_t y = stIntTuple_get(aPair, 2);
        int64_t updatedWeight = stIntTuple_get(aPair, 0) - gapGamma * (indelProbsX[x] + indelProbsY[y]);
        stList_append(reweightedAlignedPairs, stIntTuple_construct3(updatedWeight, x, y));
    }
    stList_destruct(alignedPairs);
    return reweightedAlignedPairs;
}

stList *reweightAlignedPairs2(stList *alignedPairs,
                              int64_t seqLengthX, int64_t seqLengthY,
                              double gapGamma) {
    if(gapGamma <= 0.0) {
        return alignedPairs;
    }
    int64_t *indelProbsX = getIndelProbabilities(alignedPairs, seqLengthX, 1);
    int64_t *indelProbsY = getIndelProbabilities(alignedPairs, seqLengthY, 0);
    alignedPairs = reweightAlignedPairs(alignedPairs, indelProbsX, indelProbsY, gapGamma);
    free(indelProbsX);
    free(indelProbsY);
    return alignedPairs;
}

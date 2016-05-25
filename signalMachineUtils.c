#include <stdio.h>
#include "signalMachineUtils.h"

// Borrowed from Bob Stout http://stjarnhimlen.se/snippets/strrev.c
static char *strrev(char *str) {
    char *p1, *p2;

    if (! str || ! *str)
        return str;
    for (p1 = str, p2 = str + strlen(str) - 1; p2 > p1; ++p1, --p2) {
        *p1 ^= *p2;
        *p2 ^= *p1;
        *p1 ^= *p2;
    }
    return str;
}

char *signalUtils_getSubSequence(char *seq, int64_t start, int64_t end, bool strand) {
    if (strand) {
        seq = stString_getSubString(seq, start, end - start);
        return seq;
    }
    seq = stString_getSubString(seq, end, start - end);
    return seq;
}

static inline void referenceSequence_loadReference(ReferenceSequence *self,
                                                   char *forwardReferencePath, char *backwardReferencePath) {
    if ((!stFile_exists(forwardReferencePath)) || (!stFile_exists(backwardReferencePath))) {
        st_errAbort("Couldn't find forward or backward reference files, %s, %s",
                    forwardReferencePath, backwardReferencePath);
    }

    self->reference = stFile_getLineFromFile(fopen(forwardReferencePath, "r"));
    self->reverseComplementOfReference = stFile_getLineFromFile(fopen(backwardReferencePath, "r"));
}

static inline struct PairwiseAlignment *referenceSequence_copyPairwiseAlignment(struct PairwiseAlignment *pA) {
    struct PairwiseAlignment *A = constructPairwiseAlignment(pA->contig1, pA->start1, pA->end1, pA->strand1,
                                                              pA->contig2, pA->start2, pA->end2, pA->strand2,
                                                              pA->score, NULL);
    return A;
}

static inline void referenceSequence_setTrimmedSeqeuences(ReferenceSequence *self) {
    self->trimmedForwardSequence = signalUtils_getSubSequence(self->reference, self->A->start1,
                                                              self->A->end1, self->A->strand1);
    self->trimmedBackwardSequence = strrev(signalUtils_getSubSequence(self->reverseComplementOfReference,
                                                                      self->A->start1, self->A->end1,
                                                                      self->A->strand1));

}

static inline void referenceSequence_reset(ReferenceSequence *self) {
    free(self->reference);
    free(self->reverseComplementOfReference);
    free(self->trimmedForwardSequence);
    free(self->trimmedBackwardSequence);
    self->initialized = FALSE;
}

char *referenceSequence_getTemplateTarget(ReferenceSequence *self) {
    //return self->forward ? self->trimmedForwardSequence : self->trimmedBackwardSequence;
    return self->A->strand1 ? self->trimmedForwardSequence : self->trimmedBackwardSequence;
}

char *referenceSequence_getComplementTarget(ReferenceSequence *self) {
    //return self->forward ? self->trimmedBackwardSequence : self->trimmedForwardSequence;
    return self->A->strand1 ? self->trimmedBackwardSequence : self->trimmedForwardSequence;
}

ReferenceSequence *signalUtils_ReferenceSequenceConstructFull(char *forwardReferencePath, char *backwardReferencePath,
                                                              struct PairwiseAlignment *pA) {
    ReferenceSequence *R = st_malloc(sizeof(ReferenceSequence));
    referenceSequence_loadReference(R, forwardReferencePath, backwardReferencePath);

    R->A = referenceSequence_copyPairwiseAlignment(pA);

    //R->trimmedForwardSequence = signalUtils_getSubSequence(R->reference, pA->start1, pA->end1, pA->strand1);
    //R->trimmedBackwardSequence = strrev(signalUtils_getSubSequence(R->reverseComplementOfReference,
    //                                                        pA->start1, pA->end1, pA->strand1));

    referenceSequence_setTrimmedSeqeuences(R);
    //R->trimmedForwardSequence = signalUtils_getSubSequence(R->reference, R->A->start1, R->A->end1, R->A->strand1);
    //R->trimmedBackwardSequence = strrev(signalUtils_getSubSequence(R->reverseComplementOfReference,
    //                                                               R->A->start1, R->A->end1, R->A->strand1));

    R->getTemplateTargetSequence = referenceSequence_getTemplateTarget;
    R->getComplementTargetSequence = referenceSequence_getComplementTarget;

    R->initialized = TRUE;

    return R;
}

ReferenceSequence *signalUtils_ReferenceSequenceConstructEmpty(struct PairwiseAlignment *pA) {
    ReferenceSequence *R = st_malloc(sizeof(ReferenceSequence));
    R->A = referenceSequence_copyPairwiseAlignment(pA);

    R->reference = NULL;
    R->reverseComplementOfReference = NULL;
    R->trimmedForwardSequence = NULL;
    R->trimmedBackwardSequence = NULL;

    R->getTemplateTargetSequence = referenceSequence_getTemplateTarget;
    R->getComplementTargetSequence = referenceSequence_getComplementTarget;

    R->initialized = FALSE;

    return R;
}

void signalUtils_ReferenceSequenceSet(ReferenceSequence *self,
                                      char *forwardReferencePath, char *backwardReferencePath) {
    if (self->initialized) {
        referenceSequence_reset(self);
    }

    referenceSequence_loadReference(self, forwardReferencePath, backwardReferencePath);
    referenceSequence_setTrimmedSeqeuences(self);

    self->initialized = TRUE;
}

void signalUtils_ReferenceSequenceDestruct(ReferenceSequence *self) {
    //destructPairwiseAlignment(self->A);
    free(self->A->contig1);
    free(self->A->contig2);
    free(self->A);
    free(self->reference);
    free(self->reverseComplementOfReference);
    free(self->trimmedForwardSequence);
    free(self->trimmedBackwardSequence);
    free(self);
}

void printFoo() {
    fprintf(stderr, "HELLO FOO!\n");
}


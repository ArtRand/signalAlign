#ifndef NANOPORE
#define NANOPORE
#include "sonLibTypes.h"
#define NB_EVENT_PARAMS 3

typedef struct _nanoporeReadAdjustmentParameters {
    double scale;
    double shift;
    double var;
    double scale_sd;
    double var_sd;
} NanoporeReadAdjustmentParameters;

typedef struct _nanoporeRead {
    int64_t readLength; // 2D read length in nucleotides
    int64_t nbTemplateEvents; // number of events in the array
    int64_t nbComplementEvents; // same
    NanoporeReadAdjustmentParameters templateParams;
    NanoporeReadAdjustmentParameters complementParams;

    char *twoDread; // dread indeed

    int64_t *templateEventMap;
    double *templateEvents; // mean, stdev, length

    int64_t *complementEventMap;
    double *complementEvents;
    bool scaled;
} NanoporeRead;

NanoporeRead *nanopore_loadNanoporeReadFromFile(const char *nanoporeReadFile);

stList *nanopore_remapAnchorPairs(stList *anchorPairs, int64_t *eventMap);

stList *nanopore_remapAnchorPairsWithOffset(stList *unmappedPairs, int64_t *eventMap, int64_t mapOffset);

void nanopore_descaleNanoporeRead(NanoporeRead *npRead);

void nanopore_nanoporeReadDestruct(NanoporeRead *npRead);

#endif

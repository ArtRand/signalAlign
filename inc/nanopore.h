#ifndef NANOPORE
#define NANOPORE
#include "sonLibTypes.h"
#define NB_EVENT_PARAMS 3

#ifndef MACHEP
#define MACHEP 1.11022302462515654042E-16
#endif


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

typedef struct _eventKmerTuple {
    double eventMean;
    double eventSd;
    double eventDuration;
    int64_t kmerIndex;
} EventKmerTuple;


// loads a nanopore read (.npRead) from a file
// TODO refactor format so that it can handle 1D reads also
NanoporeRead *nanopore_loadNanoporeReadFromFile(const char *nanoporeReadFile);

EventKmerTuple *nanopore_eventKmerTupleConstruct(double mean, double sd, double duration, int64_t kmerIndex);

NanoporeReadAdjustmentParameters *nanopore_readAdjustmentParametersConstruct();

stList *nanopore_remapAnchorPairs(stList *anchorPairs, int64_t *eventMap);

stList *nanopore_remapAnchorPairsWithOffset(stList *unmappedPairs, int64_t *eventMap, int64_t mapOffset);

void nanopore_descaleNanoporeRead(NanoporeRead *npRead);

stList *nanopore_getAnchorKmersToEventsMap(stList *anchorPairs, double *eventSequence, char *nucleotideSequence);

void nanopore_nanoporeReadDestruct(NanoporeRead *npRead);

void nanopore_compute_scale_params(double *model, stList *kmerToEventMap, NanoporeReadAdjustmentParameters *params,
                                   bool drift_out, bool var_out);

void nanopore_lineq_solve(double *A, double *b, double *x_out, int64_t n);

double rand_uniform2(double a, double b);

#endif

#ifndef NANOPORE
#define NANOPORE
#include "sonLibTypes.h"
#define NB_EVENT_PARAMS 4

#ifndef MACHEP
#define MACHEP 1.11022302462515654042E-16
#endif


typedef struct _nanoporeReadAdjustmentParameters {
    double scale;
    double shift;
    double var;
    double shift_sd;
    double scale_sd;
    double var_sd;
    double drift;
} NanoporeReadAdjustmentParameters;

typedef struct _nanoporeRead {
    int64_t readLength; // 2D read length in nucleotides
    int64_t templateReadLength;
    int64_t complementReadLength;
    int64_t nbTemplateEvents; // number of events in the array
    int64_t nbComplementEvents; // same
    NanoporeReadAdjustmentParameters templateParams;
    NanoporeReadAdjustmentParameters complementParams;

    char *twoDread; // dread indeed
    char *templateRead;
    char *complementRead;

    int64_t *templateEventMap;
    int64_t *templateStrandEventMap;
    double *templateEvents; // mean, stdev, length, start_time

    int64_t *complementEventMap;
    int64_t *complementStrandEventMap;
    double *complementEvents;

    // these are temporary until I come up with something better
    int64_t *templateModelState;
    double *templatePModel;

    int64_t *complementModelState;
    double *complementPModel;

    bool scaled;
} NanoporeRead;

typedef struct _eventKmerTuple {
    double eventMean;
    double eventSd;
    double deltaTime;
    int64_t kmerIndex;
} EventKmerTuple;


// loads a nanopore read (.npRead) from a file
// TODO refactor format so that it can handle 1D reads also
NanoporeRead *nanopore_loadNanoporeReadFromFile(const char *nanoporeReadFile);

EventKmerTuple *nanopore_eventKmerTupleConstruct(double mean, double sd, double deltaTime, int64_t kmerIndex);

NanoporeReadAdjustmentParameters *nanopore_readAdjustmentParametersConstruct();

stList *nanopore_remapAnchorPairs(stList *anchorPairs, int64_t *eventMap);

stList *nanopore_remapAnchorPairsWithOffset(stList *unmappedPairs, int64_t *eventMap, int64_t mapOffset);

void nanopore_descaleNanoporeRead(NanoporeRead *npRead);

stList *nanopore_getAnchorKmersToEventsMap(stList *anchorPairs, double *eventSequence, char *nucleotideSequence);

stList *nanopore_getTemplateOneDAssignments(NanoporeRead *npRead, double threshold);

stList *nanopore_getComplementOneDAssignments(NanoporeRead *npRead, double threshold);

stList *nanopore_getOneDAssignmentsFromRead(int64_t *strandEventMap, double *events, char *strandRead,
                                            int64_t readLength);

stList *nanopore_complementOneDAssignmentsFromRead(NanoporeRead *npRead, double ignore);

stList *nanopore_templateOneDAssignmentsFromRead(NanoporeRead *npRead, double ignore);

void nanopore_adjustEventsForDrift(NanoporeRead *npRead);

void nanopore_adjustTemplateEventsForDrift(NanoporeRead *npRead);

void nanopore_adjustComplementEventsForDrift(NanoporeRead *npRead);

// this is just a placeholder function that doesn't do anything. it's used in estimateNanoporeParams so that
// the program returns 'raw' un-adjusted event means
void nanopore_dontAdjustEvents(NanoporeRead *npRead);

void nanopore_nanoporeReadDestruct(NanoporeRead *npRead);

void nanopore_compute_mean_scale_params(double *model, stList *kmerToEventMap, NanoporeReadAdjustmentParameters *params, bool drift_out, bool var_out);

void nanopore_compute_noise_scale_params(double *model, stList *kmerToEventMap, NanoporeReadAdjustmentParameters *params);

void nanopore_convert_to_lognormal_params(int64_t alphabet_size, int64_t kmer_length, double *model,
                                          stList *kmerToEventMap);

void nanopore_lineq_solve(double *A, double *b, double *x_out, int64_t n);

double rand_uniform2(double a, double b);

#endif

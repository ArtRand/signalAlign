//
// Created by Arthur Rand on 6/24/16.
//

#include <getopt.h>
#include <string.h>
#include "pairwiseAligner.h"
#include "signalMachineUtils.h"

void usage() {
    fprintf(stderr, "estimateNanoporeParams: need more arguments\n");
}

double absPercentDiff(double obs, double exp) {
    double percentDiff = ((obs - exp) / exp) * 100;
    return percentDiff > 0 ? percentDiff : -percentDiff;
}

Sequence *initializeSequenceFromNpReadFile(NanoporeRead *npRead, bool templateStrand) {
    return sequence_construct2(
            (templateStrand ? npRead->nbTemplateEvents : npRead->nbComplementEvents),
            (templateStrand ? npRead->templateEvents : npRead->complementEvents),
            sequence_getEvent, sequence_sliceEventSequence, event);
}

//stList *getRemappedAnchorPairs(stList *unmappedAnchors, int64_t *eventMap, int64_t mapOffset) {
//    stList *remapedAnchors = nanopore_remapAnchorPairsWithOffset(unmappedAnchors, eventMap, mapOffset);
//    stList *filteredRemappedAnchors = filterToRemoveOverlap(remapedAnchors);
//    return filteredRemappedAnchors;
//}

Sequence *makeEventSequenceFromPairwiseAlignment(double *events, int64_t queryStart, int64_t queryEnd,
                                                 int64_t *eventMap) {
    // find the event mapped to the start and end of the 2D read alignment
    int64_t startIdx = eventMap[queryStart];
    int64_t endIdx = eventMap[queryEnd];

    // move the event pointer to the first event
    size_t elementSize = sizeof(double);
    void *elements = (char *)events + ((startIdx * NB_EVENT_PARAMS) * elementSize);

    // make the eventSequence
    Sequence *eventS = sequence_constructEventSequence(endIdx - startIdx, elements);

    return eventS;
}

void printEstimateOfParams(NanoporeRead *npRead, StateMachine *sM, double threshold,
                           NanoporeReadAdjustmentParameters *trueParams, char *readLabel, char *strand,
                           stList *(*assignmentFunction)(NanoporeRead *, double)) {
    NanoporeReadAdjustmentParameters *params = nanopore_readAdjustmentParametersConstruct();
    stList *map = assignmentFunction(npRead, threshold);
    fprintf(stderr, "Got %"PRId64" assignments\n", stList_length(map));
    nanopore_compute_scale_params(sM->EMISSION_MATCH_MATRIX, map, params, FALSE, TRUE);

    double scale_err = absPercentDiff(params->scale, trueParams->scale);
    double shift_err = absPercentDiff(params->shift, trueParams->shift);
    double var_err = absPercentDiff(params->var, trueParams->var);

    fprintf(stdout, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\n",
            params->scale, params->shift, params->var,
            trueParams->scale, trueParams->shift, trueParams->var,
            scale_err, shift_err, var_err,
            threshold, strand, readLabel);
}

int main(int argc, char *argv[]) {
    int64_t j = 0;
    char *npReadFile = NULL;
    char *templateModelFile = stString_print("../models/testModel_template.model");
    char *complementModelFile = stString_print("../models/testModel_complement.model");
    double threshold = 0.8;

    int key;
    while (1) {
        static struct option long_options[] = {
                {"help",                    no_argument,        0,  'h'},
                {"templateModel",           required_argument,  0,  'T'},
                {"complementModel",         required_argument,  0,  'C'},
                {"npRead",                  required_argument,  0,  'q'},
                {"threshold",               required_argument,  0,  'D'},
                {0, 0, 0, 0} };

        int option_index = 0;

        key = getopt_long(argc, argv, "h:T:C:q:f:b:D:m:",
                          long_options, &option_index);

        if (key == -1) {
            //usage();
            break;
        }
        switch (key) {
            case 'h':
                usage();
                return 1;
            case 'T':
                templateModelFile = stString_copy(optarg);
                break;
            case 'C':
                complementModelFile = stString_copy(optarg);
                break;
            case 'q':
                npReadFile = stString_copy(optarg);
                break;
            case 'D':
                j = sscanf(optarg, "%lf", &threshold);
                assert (j == 1);
                assert (threshold >= 0);
                break;
            default:
                usage();
                return 1;
        }
    }

    if (!stFile_exists(npReadFile)) {
        st_errAbort("Could not find npRead here: %s\n", npReadFile);
    }

    double thresholds[5] = {
            0.0,
            0.2,
            0.4,
            0.6,
            0.8,
    };

    for (int64_t i = 0; i < 5; i++) {
        NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);
        StateMachine *sMt = getStateMachine3(templateModelFile);
        StateMachine *sMc = getStateMachine3(complementModelFile);
        printEstimateOfParams(npRead, sMt, thresholds[i], &npRead->templateParams, npReadFile, "t",
                              nanopore_getTemplateOneDAssignments);
        printEstimateOfParams(npRead, sMc, thresholds[i], &npRead->complementParams, npReadFile, "c",
                              nanopore_getComplementOneDAssignments);
    }

    //fprintf(stdout, "\n");
    (void) j;  // silence unused variable warning.
    return 0;
}
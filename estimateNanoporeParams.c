//
// Created by Arthur Rand on 6/24/16.
//

#include <getopt.h>
#include <string.h>
#include "pairwiseAligner.h"
#include "signalMachineUtils.h"

#define ASSIGNMENT_THRESHOLD 0.0

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

void printEventsAndParams(NanoporeRead *npRead) {
    // kmer | mean | stDev | prob | scale | shift | var | (drift)
    for (int64_t i = 0; i < npRead->nbTemplateEvents; i++) {
        int64_t index = i * NB_EVENT_PARAMS;
        fprintf(stdout, "%"PRId64"\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n",
                npRead->templateModelState[i],
                npRead->templateEvents[index], // mean
                npRead->templateEvents[index + 1], // st_dev
                npRead->templateEvents[index + 3], // start
                npRead->templatePModel[i],
                npRead->templateParams.scale,
                npRead->templateParams.shift,
                npRead->templateParams.var,
                npRead->templateParams.drift,
                "t");
    }
    for (int64_t i = 0; i < npRead->nbComplementEvents; i++) {
        int64_t index = i * NB_EVENT_PARAMS;
        fprintf(stdout, "%"PRId64"\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n",
                npRead->complementModelState[i],
                npRead->complementEvents[index], // mean
                npRead->complementEvents[index + 1], // st_dev
                npRead->complementEvents[index + 3], // start
                npRead->complementPModel[i],
                npRead->complementParams.scale,
                npRead->complementParams.shift,
                npRead->complementParams.var,
                npRead->complementParams.drift,
                "c");
    }
}

int main(int argc, char *argv[]) {
    int64_t j = 0;
    char *npReadFile = NULL;
    char *templateModelFile = stString_print("../models/testModelR9_template.model");
    char *complementModelFile = stString_print("../models/testModelR9_complement_pop2.model");
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

    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);
    StateMachine *sMt = getStateMachine3(templateModelFile);
    StateMachine *sMc = getStateMachine3(complementModelFile);
    signalUtils_estimateNanoporeParams(sMt, npRead, &npRead->templateParams, ASSIGNMENT_THRESHOLD,
                                       nanopore_templateOneDAssignmentsFromRead);
    signalUtils_estimateNanoporeParams(sMc, npRead, &npRead->complementParams, ASSIGNMENT_THRESHOLD,
                                       nanopore_complementOneDAssignmentsFromRead);
    printEventsAndParams(npRead);

    (void) j;  // silence unused variable warning.
    return 0;
}
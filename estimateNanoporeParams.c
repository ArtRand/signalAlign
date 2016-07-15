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

stList *lineTokensFromFile(const char *filePath, int64_t getLine) {
    FILE *fH = fopen(filePath, "r");
    int64_t lineCount = 1;
    stList *tokens;
    char *string = stFile_getLineFromFile(fH);
    while (string != NULL) {
        string = stFile_getLineFromFile(fH);
        if (lineCount == getLine) {
            tokens = stString_split(string);
            return tokens;
        } else {
            lineCount++;
            free(string);
        }
    }
    fclose(fH);
    return stList_construct3(0, &free);
}

void printEstimateOfParams(NanoporeRead *npRead, StateMachine *sM, double threshold,
                           NanoporeReadAdjustmentParameters *trueParams, char *readLabel, char *strand,
                           stList *(*assignmentFunction)(NanoporeRead *, double)) {
    NanoporeReadAdjustmentParameters *params = nanopore_readAdjustmentParametersConstruct();
    stList *map = assignmentFunction(npRead, threshold);
    fprintf(stderr, "Got %"PRId64" assignments\n", stList_length(map));
    nanopore_compute_mean_scale_params(sM->EMISSION_MATCH_MATRIX, map, params, FALSE, TRUE);

    double scale_err = absPercentDiff(params->scale, trueParams->scale);
    double shift_err = absPercentDiff(params->shift, trueParams->shift);
    double var_err = absPercentDiff(params->var, trueParams->var);

    fprintf(stdout, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\n",
            params->scale, params->shift, params->var,
            trueParams->scale, trueParams->shift, trueParams->var,
            scale_err, shift_err, var_err,
            threshold, strand, readLabel);
}

void printEventMeansAndParams(NanoporeRead *npRead, stList *templateKmers, stList *complementKmers) {
    char *t_modelState, *c_modelState;
    // kmer | mean | stDev | start | prob | scale | shift | var | drift | strand
    for (int64_t i = 0; i < npRead->nbTemplateEvents; i++) {
        int64_t index = i * NB_EVENT_PARAMS;
        t_modelState = (char *)stList_get(templateKmers, i);
        fprintf(stdout, "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n",
                t_modelState,
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
        c_modelState = (char *)stList_get(complementKmers, i);
        fprintf(stdout, "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n",
                c_modelState,
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

void printEventNoisesAndParams(NanoporeRead *npRead, stList *templateKmers, stList *complementKmers) {
    char *t_modelState, *c_modelState;
    // kmer | stDev | scale_sd | shift_sd | var_sd | strand
    for (int64_t i = 0; i < npRead->nbTemplateEvents; i++) {
        int64_t index = i * NB_EVENT_PARAMS;
        t_modelState = (char *)stList_get(templateKmers, i);
        fprintf(stdout, "%s\t%f\t%f\t%f\t%f\t%s\n",
                t_modelState,
                npRead->templateEvents[index + 1], // st_dev
                npRead->templateParams.scale_sd,
                npRead->templateParams.shift_sd,
                npRead->templateParams.var_sd,
                "t");
    }
    for (int64_t i = 0; i < npRead->nbComplementEvents; i++) {
        int64_t index = i * NB_EVENT_PARAMS;
        c_modelState = (char *)stList_get(complementKmers, i);
        fprintf(stdout, "%s\t%f\t%f\t%f\t%f\t%s\n",
                c_modelState,
                npRead->complementEvents[index + 1], // st_dev
                npRead->complementParams.scale_sd,
                npRead->complementParams.shift_sd,
                npRead->complementParams.var_sd,
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
    // read in the .npRead file
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);

    // build state machines (to use the look up table)
    StateMachine *sMt = getStateMachine3(templateModelFile);
    StateMachine *sMc = getStateMachine3(complementModelFile);

    // make 1D map of events (mean, noise) to kmers
    stList *templateMap = signalUtils_templateOneDAssignmentsFromRead(npRead, sMt, ASSIGNMENT_THRESHOLD);
    stList *complementMap = signalUtils_complementOneDAssignmentsFromRead(npRead, sMc, ASSIGNMENT_THRESHOLD);

    // convert template to log normal
    nanopore_convert_to_lognormal_params(sMt->alphabetSize, sMt->kmerLength, sMt->EMISSION_MATCH_MATRIX, templateMap);
    // convert complement to log normal
    nanopore_convert_to_lognormal_params(sMc->alphabetSize, sMc->kmerLength, sMc->EMISSION_MATCH_MATRIX, complementMap);

    // error log report
    st_uglyf("SENTINEL - Before: shift_sd: %f scale_sd: %f var_sd: %f [template]\n",
             npRead->templateParams.shift_sd, npRead->templateParams.scale_sd, npRead->templateParams.var_sd);
    st_uglyf("SENTINEL - Before: shift_sd: %f scale_sd: %f var_sd: %f [template]\n",
             npRead->complementParams.shift_sd, npRead->complementParams.scale_sd, npRead->complementParams.var_sd);

    // compute template params
    nanopore_compute_noise_scale_params(sMt->EMISSION_MATCH_MATRIX, templateMap, &npRead->templateParams);
    // compute complement params
    nanopore_compute_noise_scale_params(sMc->EMISSION_MATCH_MATRIX, complementMap, &npRead->complementParams);

    // error log report
    st_uglyf("SENTINEL - After: shift_sd: %f scale_sd: %f var_sd: %f [template]\n",
             npRead->templateParams.shift_sd, npRead->templateParams.scale_sd, npRead->templateParams.var_sd);
    st_uglyf("SENTINEL - After: shift_sd: %f scale_sd: %f var_sd: %f [template]\n",
             npRead->complementParams.shift_sd, npRead->complementParams.scale_sd, npRead->complementParams.var_sd);


    //signalUtils_estimateNanoporeParams(sMt, npRead, &npRead->templateParams, ASSIGNMENT_THRESHOLD,
    //                                   signalUtils_templateOneDAssignmentsFromRead, nanopore_dontAdjustEvents);
    //signalUtils_estimateNanoporeParams(sMc, npRead, &npRead->complementParams, ASSIGNMENT_THRESHOLD,
    //                                   signalUtils_complementOneDAssignmentsFromRead, nanopore_dontAdjustEvents);

    stList *templateKmers = lineTokensFromFile(npReadFile, 10);
    stList *complementKmers = lineTokensFromFile(npReadFile, 12);
    printEventNoisesAndParams(npRead, templateKmers, complementKmers);
    //printEventMeansAndParams(npRead, templateKmers, complementKmers);

    stList_destruct(templateKmers);
    stList_destruct(complementKmers);
    stList_destruct(templateMap);
    stList_destruct(complementMap);
    nanopore_nanoporeReadDestruct(npRead);
    stateMachine_destruct(sMt);
    stateMachine_destruct(sMc);

    (void) j;  // silence unused variable warning.
    return 0;
}
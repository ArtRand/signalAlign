// Utility that builds a Nanopore HDP from an alignment, allows for experimenting with hyperparameters
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "pairwiseAligner.h"
#include "continuousHmm.h"

void usage() {
    fprintf(stderr, "USAGE: \n");
    exit(1);
}

void printStartMessage(int64_t hdpType, char *alignmentsFile, char *templateHdpOutfile, char *complementHdpOutfile) {
    fprintf(stderr, "Building Nanopore HDP\n");
    fprintf(stderr, "Making HDP type %lld\n", hdpType);
    if (alignmentsFile != NULL) {
        fprintf(stderr, "Using alignment from %s\n", alignmentsFile);
    }
    fprintf(stderr, "Putting template here: %s\n", templateHdpOutfile);
    fprintf(stderr, "Putting complement here: %s\n", complementHdpOutfile);
}

void updateHdpFromAssignments(const char *nHdpFile, const char *expectationsFile, const char *nHdpOutFile,
                              int64_t nbSamples, int64_t burnIn, int64_t thinning, bool verbose) {
    NanoporeHDP *nHdp = deserialize_nhdp(nHdpFile);
    Hmm *hdpHmm = hdpHmm_loadFromFile(expectationsFile, nHdp);
    hmmContinuous_destruct(hdpHmm, hdpHmm->type);

    fprintf(stderr, "signalAlign - Running Gibbs on HDP doing %lld samples %lld burn in %lld thinning\n",
            nbSamples, burnIn, thinning);
    execute_nhdp_gibbs_sampling(nHdp, nbSamples, burnIn, thinning, verbose);
    finalize_nhdp_distributions(nHdp);

    fprintf(stderr, "signalAlign - Serializing HDP to %s\n", nHdpOutFile);
    serialize_nhdp(nHdp, nHdpOutFile);
    destroy_nanopore_hdp(nHdp);
}

int main(int argc, char *argv[]) {
    int64_t hdpType = NULL;
    char *templateLookupTable = stString_print("../../signalAlign/models/template_median68pA.model");
    char *complementLookupTable = stString_print("../../signalAlign/models/complement_median68pA_pop2.model");
    char *alignmentsFile = NULL;
    char *templateExpectationsFile = NULL;
    char *complementExpectationsFile = NULL;
    char *templateHdpOutfile = NULL;
    char *complementHdpOutfile = NULL;

    int64_t nbSamples, burnIn, thinning, samplingGridLength, j;
    bool verbose = FALSE;

    double baseGamma = NULL_HYPERPARAMETER;
    double middleGamma = NULL_HYPERPARAMETER;
    double leafGamma = NULL_HYPERPARAMETER;

    double baseGammaAlpha = NULL_HYPERPARAMETER;
    double baseGammaBeta = NULL_HYPERPARAMETER;

    double middleGammaAlpha = NULL_HYPERPARAMETER;
    double middleGammaBeta = NULL_HYPERPARAMETER;

    double leafGammaAlpha = NULL_HYPERPARAMETER;
    double leafGammaBeta = NULL_HYPERPARAMETER;

    double samplingGridStart, samplingGridEnd;

    int key;
    while (1) {
        static struct option long_options[] = {
                {"help",                        no_argument,        0,  'h'},
                {"verbose",                     no_argument,        0,  'o'},
                {"HdpType",                     required_argument,  0,  'p'},
                {"tempalteLookupTable",         required_argument,  0,  'T'},
                {"complementLookupTable",       required_argument,  0,  'C'},
                {"alignments",                  required_argument,  0,  'l'},
                {"templateExpectations",        required_argument,  0,  'E'},
                {"complementExpectations",      required_argument,  0,  'W'},
                {"templateHdp",                 required_argument,  0,  'v'},
                {"complementHdp",               required_argument,  0,  'w'},
                {"nbSamples",                   required_argument,  0,  'n'},
                {"burnIn",                      required_argument,  0,  'I'},
                {"thinning",                    required_argument,  0,  't'},
                {"baseGamma",                   required_argument,  0,  'B'},
                {"middleGamma",                 required_argument,  0,  'M'},
                {"leafGamma",                   required_argument,  0,  'L'},
                {"baseGammaAlpha",              required_argument,  0,  'g'},
                {"baseGammaBeta",               required_argument,  0,  'r'},
                {"middleGammaAlpha",            required_argument,  0,  'j'},
                {"middleGammaBeta",             required_argument,  0,  'y'},
                {"leafGammaAlpha",              required_argument,  0,  'i'},
                {"leafGammaBeta",               required_argument,  0,  'u'},
                {"samplingGridStart",           required_argument,  0,  's'},
                {"samplingGridEnd",             required_argument,  0,  'e'},
                {"samplingGridLength",          required_argument,  0,  'k'},
                {0, 0, 0, 0} };

        int option_index = 0;
        key = getopt_long(argc, argv, "h:o:p:T:C:l:E:W:v:w:n:I:t:B:M:L:g:r:j:y:i:u:s:e:k:",
                          long_options, &option_index);
        if (key == -1) {
            //usage();
            break;
        }
        switch (key) {
            case 'h':
                usage();
                return 1;
            case 'p':
                j = sscanf(optarg, "%" PRIi64 "", &hdpType);
                assert (j == 1);
                assert (hdpType >= 0);
                assert (hdpType <= 3);
                break;
            case 'T':
                templateLookupTable = stString_copy(optarg);
                break;
            case 'C':
                complementLookupTable = stString_copy(optarg);
                break;
            case 'l':
                alignmentsFile = stString_copy(optarg);
                break;
            case 'E':
                templateExpectationsFile = stString_copy(optarg);
                break;
            case 'W':
                complementExpectationsFile = stString_copy(optarg);
                break;
            case 'v':
                templateHdpOutfile = stString_copy(optarg);
                break;
            case 'w':
                complementHdpOutfile = stString_copy(optarg);
                break;
            case 'n':
                j = sscanf(optarg, "%" PRIi64 "", &nbSamples);
                assert (j == 1);
                assert (nbSamples >= 0);
                break;
            case 'I':
                j = sscanf(optarg, "%" PRIi64 "", &burnIn);
                assert (j == 1);
                assert (burnIn >= 0);
                break;
            case 't':
                j = sscanf(optarg, "%" PRIi64 "", &thinning);
                assert (j == 1);
                assert (thinning >= 0);
                break;
            case 'o':
                verbose = TRUE;
                break;
            case 'B':
                j = sscanf(optarg, "%lf", &baseGamma);
                assert (j == 1);
                assert (baseGamma >= 0);
                break;
            case 'M':
                j = sscanf(optarg, "%lf", &middleGamma);
                assert (j == 1);
                assert (middleGamma >= 0);
                break;
            case 'L':
                j = sscanf(optarg, "%lf", &leafGamma);
                assert (j == 1);
                assert (leafGamma >= 0);
                break;
            case 'g':
                j = sscanf(optarg, "%lf", &baseGammaAlpha);
                assert (j == 1);
                assert (baseGammaAlpha >= 0);
                break;
            case 'r':
                j = sscanf(optarg, "%lf", &baseGammaBeta);
                assert (j == 1);
                assert (baseGammaBeta >= 0);
                break;
            case 'j':
                j = sscanf(optarg, "%lf", &middleGammaAlpha);
                assert (j == 1);
                assert (middleGammaAlpha >= 0);
                break;
            case 'y':
                j = sscanf(optarg, "%lf", &middleGammaBeta);
                assert (j == 1);
                assert (middleGammaBeta >= 0);
                break;
            case 'i':
                j = sscanf(optarg, "%lf", &leafGammaAlpha);
                assert (j == 1);
                assert (leafGammaAlpha >= 0);
                break;
            case 'u':
                j = sscanf(optarg, "%lf", &leafGammaBeta);
                assert (j == 1);
                assert (leafGammaBeta >= 0);
                break;
            case 's':
                j = sscanf(optarg, "%lf", &samplingGridStart);
                assert (j == 1);
                break;
            case 'e':
                j = sscanf(optarg, "%lf", &samplingGridEnd);
                assert (j == 1);
                break;
            case 'k':
                j = sscanf(optarg, "%" PRIi64 "", &samplingGridLength);
                assert (j == 1);
                assert (samplingGridLength >= 0);
                break;
            default:
                usage();
                return 1;
        }
    }

    if ((templateHdpOutfile == NULL) || (complementHdpOutfile == NULL)) {
        st_errAbort("[buildHdpUtil] ERROR: Need to specify where to put the HDP files");
    }

    printStartMessage(hdpType, alignmentsFile, templateHdpOutfile, complementHdpOutfile);

    // option for building from alignment
    if (alignmentsFile != NULL) {
        if (!((hdpType >= 0) && (hdpType <= 8))) {
            st_errAbort("Invalid HDP type");
        }
        NanoporeHdpType type = (NanoporeHdpType) hdpType;
        nanoporeHdp_buildNanoporeHdpFromAlignment(type, templateLookupTable, complementLookupTable, alignmentsFile,
                                                  templateHdpOutfile, complementHdpOutfile,
                                                  nbSamples, burnIn, thinning, verbose,
                                                  baseGamma, middleGamma, leafGamma,
                                                  baseGammaAlpha, baseGammaBeta,
                                                  middleGammaAlpha, middleGammaBeta,
                                                  leafGammaAlpha, leafGammaBeta,
                                                  samplingGridStart, samplingGridEnd, samplingGridLength);
    } else {
#pragma omp parallel sections
        {
            {
                if (templateExpectationsFile != NULL) {
                    if (templateHdpOutfile == NULL) {
                        st_errAbort("Need to provide HDP file");
                    }
                    updateHdpFromAssignments(templateHdpOutfile, templateExpectationsFile, templateHdpOutfile,
                                             nbSamples, burnIn, thinning, verbose);
                }
            }
#pragma omp section
            {
                if (complementExpectationsFile != NULL) {
                    if (complementHdpOutfile == NULL) {
                        st_errAbort("Need to provide HDP file");
                    }
                    updateHdpFromAssignments(complementHdpOutfile, complementExpectationsFile, complementHdpOutfile,
                                             nbSamples, burnIn, thinning, verbose);
                }
            }
        }
        return 0;
    }
}
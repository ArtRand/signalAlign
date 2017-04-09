// Utility that builds a Nanopore HDP from an alignment, allows for experimenting with hyperparameters
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "pairwiseAligner.h"
#include "continuousHmm.h"

void usage() {
    fprintf(stderr, "\n\tbuildHdpUtil - make a new hierarchical Dirichlet process model from an alignment file\n\n");
    fprintf(stderr, "--help: display this message and exit\n");
    fprintf(stderr, "--verbose: enable verbose Gibbs sampling output\n");
    fprintf(stderr, "--oneD: flag for 1-D reads, only make a 'template' model\n");
    fprintf(stderr, "-a: kmer length to use in the model (usually 5 or 6)\n");
    fprintf(stderr, "-p: HDP type. Enum for the type of model to build (alphabet the model uses)\n");
    fprintf(stderr, "\tMost people will use 10, 11, 12, or 13\n");
    fprintf(stderr, "\t\t0: single level fixed (ACEGOT)\n");
    fprintf(stderr, "\t\t1: single level prior (ACEGOT)\n");
    fprintf(stderr, "\t\t2: multiset fixed (ACEGOT)\n");
    fprintf(stderr, "\t\t3: multiset prior (ACEGOT)\n");
    fprintf(stderr, "\t\t4: composiition fixed (ACEGOT)\n");
    fprintf(stderr, "\t\t5: composiition prior (ACEGOT)\n");
    fprintf(stderr, "\t\t6: middle nucleodides fixed (ACEGOT)\n");
    fprintf(stderr, "\t\t7: middle nucleodides prior (ACEGOT)\n");
    fprintf(stderr, "\t\t8: group multiset fixed (ACEGOT)\n");
    fprintf(stderr, "\t\t9: group multiset prior (ACEGOT)\n");
    fprintf(stderr, "\t\t10: single level prior (ACEGT)\n");
    fprintf(stderr, "\t\t11: multiset prior (ACEGT)\n");
    fprintf(stderr, "\t\t12: single level prior (ecoli) (ACEGIT)\n");
    fprintf(stderr, "\t\t13: multiset prior (ecoli) (ACEGIT)\n");
    fprintf(stderr, "-T: template lookup table (signalAlign HMM format)\n");
    fprintf(stderr, "-C: complement lookup table (signalAlign HMM format)\n");
    fprintf(stderr, "-l: alignments to take k-mer assignments from\n");
    fprintf(stderr, "-v: template HDP output file path\n");
    fprintf(stderr, "-w: complement HDP output file path\n");
    fprintf(stderr, "-n: number of samples\n");
    fprintf(stderr, "-I: burn in, samples to discard at the start\n");
    fprintf(stderr, "-t: thinning\n");
    fprintf(stderr, "-B: base gamma\n");
    fprintf(stderr, "-M: middle gamma\n");
    fprintf(stderr, "-L: leaf gamma\n");
    fprintf(stderr, "-g: base gamma alpha\n");
    fprintf(stderr, "-r: base gamma beta\n");
    fprintf(stderr, "-j: middle gamma alpha\n");
    fprintf(stderr, "-y: middle gamma beta\n");
    fprintf(stderr, "-i: leaf gamma alpha\n");
    fprintf(stderr, "-u: leaf gamma beta\n");
    fprintf(stderr, "-s: sampling grid start\n");
    fprintf(stderr, "-e: sampling grid end\n");
    fprintf(stderr, "-k: sampling grid length\n");
    exit(1);
}

void printStartMessage(int64_t hdpType, char *alignmentsFile, char *templateHdpOutfile, char *complementHdpOutfile, bool twoD) {
    fprintf(stderr, "Building Nanopore HDP\n");
    fprintf(stderr, "Making HDP type %"PRId64"\n", hdpType);
    if (alignmentsFile != NULL) {
        fprintf(stderr, "Using alignment from %s\n", alignmentsFile);
    }
    fprintf(stderr, "Putting template here: %s\n", templateHdpOutfile);
    if (twoD) {
        fprintf(stderr, "Putting complement here: %s\n", complementHdpOutfile);
    }
}

void updateHdpFromAssignments(const char *nHdpFile, const char *expectationsFile, const char *nHdpOutFile,
                              int64_t nbSamples, int64_t burnIn, int64_t thinning, bool verbose) {
    NanoporeHDP *nHdp = deserialize_nhdp(nHdpFile);
    Hmm *hdpHmm = hdpHmm_loadFromFile(expectationsFile, threeStateHdp, nHdp);
    hmmContinuous_destruct(hdpHmm, hdpHmm->type);

    fprintf(stderr, "signalAlign - Running Gibbs on HDP doing %"PRId64" samples %"PRId64"burn in %"PRId64"thinning\n",
            nbSamples, burnIn, thinning);
    execute_nhdp_gibbs_sampling(nHdp, nbSamples, burnIn, thinning, verbose);
    finalize_nhdp_distributions(nHdp);

    fprintf(stderr, "signalAlign - Serializing HDP to %s\n", nHdpOutFile);
    serialize_nhdp(nHdp, nHdpOutFile);
    destroy_nanopore_hdp(nHdp);
}

int main(int argc, char *argv[]) {
    int64_t hdpType = -1;
    char *templateLookupTable = NULL;
    char *complementLookupTable = NULL;
    char *alignmentsFile = NULL;
    char *templateHdpOutfile = NULL;
    char *complementHdpOutfile = NULL;

    int64_t nbSamples, burnIn, thinning, samplingGridLength, kmerLength, j;
    bool verbose = FALSE;
    bool twoD = TRUE;

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
                {"oneD",                        no_argument,        0,  'q'},
                {"kmerLength",                  required_argument,  0,  'a'},
                {"HdpType",                     required_argument,  0,  'p'},
                {"templateLookupTable",         required_argument,  0,  'T'},
                {"complementLookupTable",       required_argument,  0,  'C'},
                {"alignments",                  required_argument,  0,  'l'},
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
        key = getopt_long(argc, argv, "h:a:o:q:p:T:C:l:v:w:n:I:t:B:M:L:g:r:j:y:i:u:s:e:k:",
                          long_options, &option_index);
        if (key == -1) {
            //usage();
            break;
        }
        switch (key) {
            case 'h':
                usage();
                return 1;
            case 'a':
                j = sscanf(optarg, "%" PRIi64 "", &kmerLength);
                assert(j == 1);
                assert(kmerLength > 0);
                break;
            case 'p':
                j = sscanf(optarg, "%" PRIi64 "", &hdpType);
                assert (j == 1);
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
            case 'q':
                twoD = FALSE;
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
    
    (void) j;
    
    if ((templateHdpOutfile == NULL) || (complementHdpOutfile == NULL && twoD)) {
        st_errAbort("[buildHdpUtil] ERROR: Need to specify where to put the HDP files");
    }

    if ((templateLookupTable == NULL) || (complementLookupTable == NULL && twoD)) {
        st_errAbort("[buildHdpUtil] ERROR: Need lookup tables");
    }
    printStartMessage(hdpType, alignmentsFile, templateHdpOutfile, complementHdpOutfile, twoD);
    
    if (alignmentsFile == NULL) st_errAbort("[buildHdpUtil]Need to provide build alignment (assignments)");
    // option for building from alignment
    if (twoD) {
        if (!((hdpType >= 0) && (hdpType <= 13))) {
            st_errAbort("Invalid HDP type");
        }
        NanoporeHdpType type = (NanoporeHdpType) hdpType;
        nanoporeHdp_buildNanoporeHdpFromAlignment(type, kmerLength,
                                                  templateLookupTable, complementLookupTable, alignmentsFile,
                                                  templateHdpOutfile, complementHdpOutfile,
                                                  nbSamples, burnIn, thinning, verbose,
                                                  baseGamma, middleGamma, leafGamma,
                                                  baseGammaAlpha, baseGammaBeta,
                                                  middleGammaAlpha, middleGammaBeta,
                                                  leafGammaAlpha, leafGammaBeta,
                                                  samplingGridStart, samplingGridEnd, samplingGridLength);
            
    } else {
        if (hdpType != singleLevelPrior2 && hdpType != multisetPrior2) {
            st_errAbort("Invalid HDP type for 1D %i", hdpType);
        }
        NanoporeHdpType type = (NanoporeHdpType) hdpType;
        nanoporeHdp_buildOneDHdpFromAlignment(type, kmerLength,
                                              templateLookupTable,
                                              alignmentsFile,
                                              templateHdpOutfile,
                                              nbSamples, burnIn, thinning, verbose,
                                              baseGamma, middleGamma, leafGamma,
                                              baseGammaAlpha, baseGammaBeta,
                                              middleGammaAlpha, middleGammaBeta,
                                              leafGammaAlpha, leafGammaBeta,
                                              samplingGridStart, samplingGridEnd, samplingGridLength);
                                          
    }

    return 0;
}

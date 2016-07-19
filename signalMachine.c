#include <getopt.h>
#include <string.h>
#include "signalMachineUtils.h"
#include "pairwiseAligner.h"

#define STEP 6  // space between degenerate nucleotides in for error correction
#define ESTIMATE_PARAMS 1
#define ASSIGNMENT_THRESHOLD 0.1

void usage() {
    fprintf(stderr, "signalMachine binary, meant to be used through the signalAlign program.\n");
    fprintf(stderr, "See doc for runSignalAlign for help\n");
}

void printPairwiseAlignmentSummary(struct PairwiseAlignment *pA) {
    st_uglyf("contig 1: %s\n", pA->contig1);
    st_uglyf("strand 1: %lld\n", pA->strand1);
    st_uglyf("start  1: %lld\n", pA->start1);
    st_uglyf("end    1: %lld\n", pA->end1);
    st_uglyf("contig 2: %s\n", pA->contig2);
    st_uglyf("strand 2: %lld\n", pA->strand2);
    st_uglyf("start  2: %lld\n", pA->start2);
    st_uglyf("end    2: %lld\n", pA->end2);
}

void writePosteriorProbs(char *posteriorProbsFile, char *readFile, double *matchModel,
                         NanoporeReadAdjustmentParameters npp, double *events, char *target, bool forward,
                         char *contig, StateMachineType type, int64_t eventSequenceOffset,
                         int64_t referenceSequenceOffset, stList *alignedPairs, Strand strand) {
    // label for tsv output
    char *strandLabel;
    if (strand == template) {
        strandLabel = "t";
    }
    if (strand == complement) {
        strandLabel = "c";
    }

    // open the file for output
    FILE *fH = fopen(posteriorProbsFile, "a");
    for(int64_t i = 0; i < stList_length(alignedPairs); i++) {
        // grab the aligned pair
        stIntTuple *aPair = stList_get(alignedPairs, i);

        if (stIntTuple_length(aPair) != 4) {
            st_errAbort("Aligned pair tuples should have length 4, this one has length %lld\n",
                        stIntTuple_length(aPair));
        }

        int64_t x_i = stIntTuple_get(aPair, 1);
        int64_t x_adj = 0;  // x is the reference coordinate that we record in the aligned pairs
        if ((strand == template && forward) || (strand == complement && (!forward))) {
            x_adj = stIntTuple_get(aPair, 1) + referenceSequenceOffset;
        }
        if ((strand == complement && forward) || (strand == template && (!forward))) {
            int64_t refLength = (int64_t)strlen(target);
            int64_t refLengthInEvents = refLength - KMER_LENGTH;
            x_adj = refLengthInEvents - (x_i + (refLength - referenceSequenceOffset));
        }
        int64_t y = stIntTuple_get(aPair, 2) + eventSequenceOffset;             // event index
        double p = ((double)stIntTuple_get(aPair, 0)) / PAIR_ALIGNMENT_PROB_1;  // posterior prob
        char *pathKmer = (char *)stIntTuple_get(aPair, 3);

        // get the observations from the events
        double eventMean = events[(y * NB_EVENT_PARAMS)];
        double eventNoise = events[(y * NB_EVENT_PARAMS) + 1];
        double eventDuration = events[(y * NB_EVENT_PARAMS) + 2];

        // make the kmer string at the target index,
        char *k_i = st_malloc(KMER_LENGTH * sizeof(char));
        for (int64_t k = 0; k < KMER_LENGTH; k++) {
            k_i[k] = *(target + (x_i + k));
        }

        // get the kmer index // todo could make this the other function
        int64_t targetKmerIndex = emissions_discrete_getKmerIndexFromPtr(pathKmer);

        // get the expected event mean amplitude and noise
        double E_mean = matchModel[1 + (targetKmerIndex * MODEL_PARAMS)];
        double E_noise = matchModel[1 + (targetKmerIndex * MODEL_PARAMS + 2)];
        double scaled_Emean = E_mean * npp.scale + npp.shift;
        double scaled_Enoise = E_noise * npp.scale_sd;
        double descaledEventMean = emissions_signal_descaleEventMean_JordanStyle(eventMean, E_mean,
                                                                                 npp.scale, npp.shift, npp.var);

        // make reference kmer
        char *refKmer = NULL;
        if ((strand == template && forward) || (strand == complement && (!forward))) {
            refKmer = k_i;
        }
        if ((strand == complement && forward) || (strand == template && (!forward))) {
            refKmer = stString_reverseComplementString(k_i);
        }
        // write to disk
        fprintf(fH, "%s\t%"PRId64"\t%s\t%s\t%s\t%"PRId64"\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%s\n",
                contig, x_adj, refKmer, readFile, strandLabel, y, eventMean, eventNoise, eventDuration, k_i,
                scaled_Emean, scaled_Enoise, p, descaledEventMean, E_mean, pathKmer);

        // cleanup
        free(k_i);
    }
    fclose(fH);
}

void writePosteriorProbsSparse(char *posteriorProbsFile, char *readFile, char *target, bool forward, char *contig,
                               int64_t eventSequenceOffset, int64_t referenceSequenceOffset,
                               stList *alignedPairs, Strand strand) {
    /// / label for tsv output
    char *strandLabel;
    if (strand == template) {
        strandLabel = "t";
    }
    if (strand == complement) {
        strandLabel = "c";
    }

    // open the file for output
    FILE *fH = fopen(posteriorProbsFile, "a");
    for(int64_t i = 0; i < stList_length(alignedPairs); i++) {
        // grab the aligned pair
        stIntTuple *aPair = stList_get(alignedPairs, i);

        if (stIntTuple_length(aPair) != 4) {
            st_errAbort("Aligned pair tuples should have length 4, this one has length %lld\n",
                        stIntTuple_length(aPair));
        }

        int64_t x_i = stIntTuple_get(aPair, 1);
        int64_t x_adj = 0;  // x is the reference coordinate that we record in the aligned pairs
        if ((strand == template && forward) || (strand == complement && (!forward))) {
            x_adj = stIntTuple_get(aPair, 1) + referenceSequenceOffset;
        }
        if ((strand == complement && forward) || (strand == template && (!forward))) {
            int64_t refLength = (int64_t)strlen(target);
            int64_t refLengthInEvents = refLength - KMER_LENGTH;
            x_adj = refLengthInEvents - (x_i + (refLength - referenceSequenceOffset));
        }
        int64_t y = stIntTuple_get(aPair, 2) + eventSequenceOffset;             // event index
        double p = ((double)stIntTuple_get(aPair, 0)) / PAIR_ALIGNMENT_PROB_1;  // posterior prob
        char *pathKmer = (char *)stIntTuple_get(aPair, 3);

        // make the kmer string at the target index,
        char *k_i = st_malloc(KMER_LENGTH * sizeof(char));
        for (int64_t k = 0; k < KMER_LENGTH; k++) {
            k_i[k] = *(target + (x_i + k));
        }

        // get the expected event mean amplitude and noise
        // make reference kmer
        char *refKmer = NULL;
        if ((strand == template && forward) || (strand == complement && (!forward))) {
            refKmer = k_i;
        }
        if ((strand == complement && forward) || (strand == template && (!forward))) {
            refKmer = stString_reverseComplementString(k_i);
        }
        fprintf(fH, "%s\t%"PRId64"\t%s\t%s\t%s\t%"PRId64"\t%s\t%f\t%s\n",
                contig, x_adj, refKmer, readFile, strandLabel, y, k_i, p, pathKmer);
        // cleanup
        free(k_i);
    }
    fclose(fH);
}

StateMachine *buildStateMachine(const char *modelFile, NanoporeReadAdjustmentParameters npp, StateMachineType type,
                                NanoporeHDP *nHdp) {
    if ((type != threeState) && (type != threeStateHdp)) {
        st_errAbort("signalAlign - incompatible stateMachine type request");
    }

    if (type == threeState) {
        StateMachine *sM = getStateMachine3_descaled(modelFile, npp, !ESTIMATE_PARAMS);
        return sM;
    }
    if (type == threeStateHdp) {
        StateMachine *sM = getHdpStateMachine(nHdp, modelFile, npp);
        return sM;
    }
    else {
        st_errAbort("signalAlign - ERROR: buildStateMachine, didn't get correct input\n");
    }
    return 0;
}

inline void loadHmmRoutine(const char *hmmFile, StateMachine *sM, StateMachineType type, Hmm *expectations) {
    if ((type != threeState) && (type != threeStateHdp)) {
        st_errAbort("LoadSignalHmm : unupported stateMachineType");
    }
    hmmContinuous_loadSignalHmm(hmmFile, sM, type, expectations);
}

StateMachine *buildStateMachineAndLoadHmm(const char *modelFile, NanoporeReadAdjustmentParameters npp,
                                          StateMachineType type, NanoporeHDP *nHdp) {
    StateMachine *sM = buildStateMachine(modelFile, npp, type, nHdp);
    // commented out because now the model file has the transitions and the event model, so no longer need to
    // load the .hmm into the stateMachine
    //if (HmmFile != NULL) {
    //    loadHmmRoutine(HmmFile, sM, sM->type, hmmExpectations);
    //}
    return sM;
}

void updateHdpFromAssignments(const char *nHdpFile, const char *expectationsFile, const char *nHdpOutFile) {
    NanoporeHDP *nHdp = deserialize_nhdp(nHdpFile);
    Hmm *hdpHmm = hdpHmm_loadFromFile(expectationsFile, threeStateHdp, nHdp);
    hmmContinuous_destruct(hdpHmm, hdpHmm->type);

    fprintf(stderr, "signalAlign - Running Gibbs on HDP\n");
    execute_nhdp_gibbs_sampling(nHdp, 10000, 100000, 100, FALSE);
    finalize_nhdp_distributions(nHdp);

    fprintf(stderr, "signalAlign - Serializing HDP to %s\n", nHdpOutFile);
    serialize_nhdp(nHdp, nHdpOutFile);
    destroy_nanopore_hdp(nHdp);
}

static double totalScore(stList *alignedPairs) {
    double score = 0.0;
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        score += stIntTuple_get(aPair, 0);
    }
    return score;
}

double scoreByPosteriorProbabilityIgnoringGaps(stList *alignedPairs) {
    /*
     * Gives the average posterior match probability per base of the two sequences, ignoring indels.
     */
    return 100.0 * totalScore(alignedPairs) / ((double) stList_length(alignedPairs) * PAIR_ALIGNMENT_PROB_1);
}

stList *performSignalAlignment(StateMachine *sM, Sequence *eventSequence, int64_t *eventMap,
                               int64_t mapOffset, char *target, PairwiseAlignmentParameters *p,
                               stList *unmappedAnchors, DegenerateType degenerate) {
    if ((sM->type != threeState) && (sM->type != threeStateHdp)) {
        st_errAbort("signalAlign - You're trying to do the wrong king of alignment");
    }

    int64_t lX = sequence_correctSeqLength(strlen(target), kmer, sM->kmerLength);

    // remap anchor pairs
    stList *filteredRemappedAnchors = signalUtils_getRemappedAnchorPairs(unmappedAnchors, eventMap, mapOffset);

    // make sequences
    Sequence *sX = sequence_constructKmerSequence(
            lX, target, sequence_getKmer, sequence_sliceNucleotideSequence,
            (degenerate == canonicalVariants ? CANONICAL_NUCLEOTIDES :
             (degenerate == cytosineMethylation2 ? TWO_CYTOSINES : THREE_CYTOSINES)),
            (degenerate == canonicalVariants ? NB_CANONICAL_BASES :
             (degenerate == cytosineMethylation2 ? (NB_CYTOSINE_OPTIONS - 1) : NB_CYTOSINE_OPTIONS)),
            kmer);

    // do alignment
    stList *alignedPairs = getAlignedPairsUsingAnchors(sM, sX, eventSequence, filteredRemappedAnchors, p,
                                                       diagonalCalculationPosteriorMatchProbs, 1, 1);

    return alignedPairs;
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

void getSignalExpectations(StateMachine *sM, const char *model, const char *inputHmm, NanoporeHDP *nHdp,
                           Hmm *hmmExpectations, StateMachineType type,
                           NanoporeReadAdjustmentParameters npp, Sequence *eventSequence,
                           int64_t *eventMap, int64_t mapOffset, char *trainingTarget, PairwiseAlignmentParameters *p,
                           stList *unmappedAnchors, Strand strand, DegenerateType degenerate) {
    // correct sequence length
    int64_t lX = sequence_correctSeqLength(strlen(trainingTarget), event, sM->kmerLength);

    // remap the anchors
    stList *filteredRemappedAnchors = signalUtils_getRemappedAnchorPairs(unmappedAnchors, eventMap, mapOffset);

    Sequence *target = sequence_constructKmerSequence(
            lX, trainingTarget, sequence_getKmer, sequence_sliceNucleotideSequence,
            (degenerate == canonicalVariants ? CANONICAL_NUCLEOTIDES :
             (degenerate == cytosineMethylation2 ? TWO_CYTOSINES : THREE_CYTOSINES)),
            (degenerate == canonicalVariants ? NB_CANONICAL_BASES :
             (degenerate == cytosineMethylation2 ? (NB_CYTOSINE_OPTIONS - 1) : NB_CYTOSINE_OPTIONS)),
            kmer);

    getExpectationsUsingAnchors(sM, hmmExpectations, target, eventSequence, filteredRemappedAnchors, p,
                                diagonalCalculation_Expectations, 1, 1);
    //stateMachine_destruct(sM);
}

int main(int argc, char *argv[]) {
    StateMachineType sMtype = threeState;
    int64_t j = 0;
    int64_t diagExpansion = 50;
    double threshold = 0.01;
    int64_t constraintTrim = 14;
    int64_t degenerate;

    char *templateModelFile = stString_print("../models/testModel_template.model");
    char *complementModelFile = stString_print("../models/testModel_complement.model");
    char *readLabel = NULL;
    char *npReadFile = NULL;
    char *forwardReference = NULL;
    char *backwardReference = NULL;
    char *errorCorrectPath = NULL;
    char *posteriorProbsFile = NULL;
    char *templateHmmFile = NULL;
    char *complementHmmFile = NULL;
    char *templateExpectationsFile = NULL;
    char *complementExpectationsFile = NULL;
    char *templateHdp = NULL;
    char *complementHdp = NULL;
    bool sparseOutput = FALSE;

    int key;
    while (1) {
        static struct option long_options[] = {
                {"help",                    no_argument,        0,  'h'},
                {"sm3Hdp",                  no_argument,        0,  'd'},
                {"sparse_output",           no_argument,        0,  's'},
                {"degenerate",              required_argument,  0,  'o'},
                {"templateModel",           required_argument,  0,  'T'},
                {"complementModel",         required_argument,  0,  'C'},
                {"readLabel",               required_argument,  0,  'L'},
                {"npRead",                  required_argument,  0,  'q'},
                {"forward_reference",       required_argument,  0,  'f'},
                {"backward_reference",      required_argument,  0,  'b'},
                {"error_correct_path",      required_argument,  0,  'p'},
                {"posteriors",              required_argument,  0,  'u'},
                {"inTemplateHmm",           required_argument,  0,  'y'},
                {"inComplementHmm",         required_argument,  0,  'z'},
                {"templateHdp",             required_argument,  0,  'v'},
                {"complementHdp",           required_argument,  0,  'w'},
                {"templateExpectations",    required_argument,  0,  't'},
                {"complementExpectations",  required_argument,  0,  'c'},
                {"diagonalExpansion",       required_argument,  0,  'x'},
                {"threshold",               required_argument,  0,  'D'},
                {"constraintTrim",          required_argument,  0,  'm'},
                {0, 0, 0, 0} };

        int option_index = 0;

        key = getopt_long(argc, argv, "h:sd:o:p:a:T:C:L:q:f:b:p:u:y:z:v:w:t:c:x:D:m:",
                          long_options, &option_index);

        if (key == -1) {
            //usage();
            break;
        }
        switch (key) {
            case 'h':
                usage();
                return 1;
            case 's':
                sparseOutput = TRUE;
                break;
            case 'o':
                j = sscanf(optarg, "%" PRIi64 "", &degenerate);
                assert (j == 1);
                assert (constraintTrim >= 0);
                //constraintTrim = (int64_t)degenerate;
                break;
            case 'd':
                sMtype = threeStateHdp;
                break;
            case 'T':
                templateModelFile = stString_copy(optarg);
                break;
            case 'C':
                complementModelFile = stString_copy(optarg);
                break;
            case 'L':
                readLabel = stString_copy(optarg);
                break;
            case 'q':
                npReadFile = stString_copy(optarg);
                break;
            case 'f':
                forwardReference = stString_copy(optarg);
                break;
            case 'b':
                backwardReference= stString_copy(optarg);
                break;
            case 'p':
                errorCorrectPath = stString_copy(optarg);
                break;
            case 'u':
                posteriorProbsFile = stString_copy(optarg);
                break;
            case 't':
                templateExpectationsFile = stString_copy(optarg);
                break;
            case 'c':
                complementExpectationsFile = stString_copy(optarg);
                break;
            case 'y':
                templateHmmFile = stString_copy(optarg);
                break;
            case 'z':
                complementHmmFile = stString_copy(optarg);
                break;
            case 'v':
                templateHdp = stString_copy(optarg);
                break;
            case 'w':
                complementHdp = stString_copy(optarg);
                break;
            case 'x':
                j = sscanf(optarg, "%" PRIi64 "", &diagExpansion);
                assert (j == 1);
                assert (diagExpansion >= 0);
                diagExpansion = (int64_t)diagExpansion;
                break;
            case 'D':
                j = sscanf(optarg, "%lf", &threshold);
                assert (j == 1);
                assert (threshold >= 0);
                break;
            case 'm':
                j = sscanf(optarg, "%" PRIi64 "", &constraintTrim);
                assert (j == 1);
                assert (constraintTrim >= 0);
                constraintTrim = (int64_t)constraintTrim;
                break;
            default:
                usage();
                return 1;
        }
    }

    (void) j;  // silence unused variable warning.

    // Anchors //
    // get pairwise alignment from stdin, in exonerate CIGAR format
    FILE *fileHandleIn = stdin;
    // parse input CIGAR to get anchors
    struct PairwiseAlignment *pA;
    pA = cigarRead(fileHandleIn);

    // Alignment Parameters //
    // make the pairwise alignment parameters
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
    p->threshold = threshold;
    p->constraintDiagonalTrim = constraintTrim;
    p->diagonalExpansion = diagExpansion;
    st_uglyf("contstructed pA\n");
    // HDP routines //
    // load HDPs
    NanoporeHDP *nHdpT, *nHdpC;
    // check
    if ((templateHdp != NULL) || (complementHdp != NULL)) {
        if ((templateHdp == NULL) || (complementHdp == NULL)) {
            st_errAbort("Need to have template and complement HDPs");
        }
        if (sMtype != threeStateHdp) {
            fprintf(stderr, "[signalAlign] - Warning: this kind of stateMachine does not use the HDPs you gave\n");
        }
        fprintf(stderr, "[signalAlign] - using NanoporeHDPs\n");
    }

    #pragma omp parallel sections
    {
        {
            nHdpT = (templateHdp == NULL) ? NULL : deserialize_nhdp(templateHdp);
        }

        #pragma omp section
        {
            nHdpC = (complementHdp == NULL) ? NULL : deserialize_nhdp(complementHdp);
        }
    }

    ReferenceSequence *R;
    if (errorCorrectPath == NULL) { // not doing error correction
        if ((forwardReference == NULL) || (backwardReference == NULL)) {
            st_errAbort("[signalAlign] - ERROR: did not get reference files %s %s\n",
                        forwardReference, backwardReference);
        }
        R = signalUtils_ReferenceSequenceConstructFull(forwardReference, backwardReference, pA);
    } else {
        R = signalUtils_ReferenceSequenceConstructEmpty(pA);
    }
    st_uglyf("constructed reference sequences\n");
    // Nanopore Read //
    // load nanopore read
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);
    st_uglyf("Loaded nanoporeRead\n:");
    // constrain the event sequence to the positions given by the guide alignment
    Sequence *tEventSequence = makeEventSequenceFromPairwiseAlignment(npRead->templateEvents,
                                                                      pA->start2, pA->end2,
                                                                      npRead->templateEventMap);

    Sequence *cEventSequence = makeEventSequenceFromPairwiseAlignment(npRead->complementEvents,
                                                                      pA->start2, pA->end2,
                                                                      npRead->complementEventMap);
    st_uglyf("Constructed event sequences\n");
    // the aligned pairs start at (0,0) so we need to correct them based on the guide alignment later.
    // record the pre-zeroed alignment start and end coordinates here
    // for the events:
    int64_t tCoordinateShift = npRead->templateEventMap[pA->start2];
    int64_t cCoordinateShift = npRead->complementEventMap[pA->start2];

    // and for the reference:
    int64_t rCoordinateShift_t = pA->start1;
    int64_t rCoordinateShift_c = pA->end1;
    bool forward = pA->strand1;  // keep track of whether this is a forward mapped read or not

    stList *anchorPairs = signalUtils_guideAlignmentToRebasedAnchorPairs(pA, p);  // pA gets modified here, no turning back

    if (errorCorrectPath != NULL) {
        st_uglyf("Starting error correcting routine\n");

        //stList *aP = stList_construct3(0, &free);
        //char *fR, *bR, *perIterationOutputFile;
        //stList *templateAlignedPairs, *complementAlignedPairs;
//        #pragma omp parallel for
        for (int64_t i = 0; i < STEP; i++) {

            StateMachine *sMt = buildStateMachine(templateModelFile, npRead->templateParams, sMtype, nHdpT);
            StateMachine *sMc = buildStateMachine(complementModelFile, npRead->complementParams, sMtype, nHdpC);

            if (posteriorProbsFile == NULL) {
                st_errAbort("SignalAlign - didn't find output file path\n");
            }

            // load the first forward and backward reference sequences
            char *fR = stString_print("%sforward_sub%i.txt", errorCorrectPath, i);
            char *bR = stString_print("%sbackward_sub%i.txt", errorCorrectPath, i);

            if (!stFile_exists(fR) || !stFile_exists(bR)) {
                st_errAbort("Error finding error correct references %s %s\n", fR, bR);
            }
            // set referenceSequence to this iteration's sequence
            signalUtils_ReferenceSequenceSet(R, fR, bR);

            fprintf(stderr, "signalAlign - starting template alignment round %lld\n", i);

            // get aligned pairs
            stList *templateAlignedPairs = performSignalAlignment(sMt, tEventSequence, npRead->templateEventMap,
                                                          pA->start2, R->getTemplateTargetSequence(R),
                                                          p, anchorPairs,
                                                          degenerate);

            double templatePosteriorScore = scoreByPosteriorProbabilityIgnoringGaps(templateAlignedPairs);

            fprintf(stdout, "%s :: Iteration %lld, # alignedPairs (template): %lld, score: %f\n",
                    readLabel, i, stList_length(templateAlignedPairs), templatePosteriorScore);

            stList_sort(templateAlignedPairs, sortByXPlusYCoordinate2); //Ensure the coordinates are increasing

            char *perIterationOutputFile = stString_print("%s.%i", posteriorProbsFile, i);

            // write to file
            if (sparseOutput) {
                writePosteriorProbsSparse(perIterationOutputFile, readLabel, R->getTemplateTargetSequence(R),
                                          forward, pA->contig1, tCoordinateShift,
                                          rCoordinateShift_t, templateAlignedPairs, template);
            } else {
                writePosteriorProbs(perIterationOutputFile, readLabel, sMt->EMISSION_MATCH_MATRIX,
                                    npRead->templateParams, npRead->templateEvents, R->getTemplateTargetSequence(R),
                                    forward, pA->contig1, sMt->type, tCoordinateShift, rCoordinateShift_t,
                                    templateAlignedPairs, template);
            }

            fprintf(stderr, "signalAlign - starting complement alignment round %lld\n", i);

            stList *complementAlignedPairs = performSignalAlignment(sMc, cEventSequence,
                                                            npRead->complementEventMap, pA->start2,
                                                            R->getComplementTargetSequence(R),
                                                            p, anchorPairs, degenerate);

            double complementPosteriorScore = scoreByPosteriorProbabilityIgnoringGaps(complementAlignedPairs);

            // sort
            stList_sort(complementAlignedPairs, sortByXPlusYCoordinate2); //Ensure the coordinates are increasing

            fprintf(stdout, "%s :: Iteration %lld, # alignedPairs (complement): %lld, score: %f\n",
                    readLabel, i, stList_length(complementAlignedPairs), complementPosteriorScore);

            if (sparseOutput) {
                writePosteriorProbsSparse(perIterationOutputFile, readLabel, R->getComplementTargetSequence(R),
                                          forward, pA->contig1, cCoordinateShift,
                                          rCoordinateShift_c, complementAlignedPairs, complement);
            } else {
                writePosteriorProbs(perIterationOutputFile, readLabel, sMc->EMISSION_MATCH_MATRIX,
                                    npRead->complementParams, npRead->complementEvents,
                                    R->getComplementTargetSequence(R), forward, pA->contig1,
                                    sMc->type, cCoordinateShift, rCoordinateShift_c,
                                    complementAlignedPairs, complement);
            }
            stList_destruct(templateAlignedPairs);
            stList_destruct(complementAlignedPairs);
            stateMachine_destruct(sMt);
            stateMachine_destruct(sMc);
            free(perIterationOutputFile);
            free(fR);
            free(bR);
        }
//        #pragma omp critical
        {
            signalUtils_ReferenceSequenceDestruct(R);
            sequence_destruct(tEventSequence);
            sequence_destruct(cEventSequence);
            destructPairwiseAlignment(pA);
        }
        return 0;
    } else if ((templateExpectationsFile != NULL) && (complementExpectationsFile != NULL)) {
        st_uglyf("Starting expectations routine\n");
        // Expectation Routine //

        if (templateHmmFile == NULL) {
            st_errAbort("[signalMachine] ERROR: need to have input HMMs\n");
        }

        StateMachine *sMt = buildStateMachine(templateModelFile, npRead->templateParams, sMtype, nHdpT);

        // temporary way to 'turn off' estimates if I want to
        if (ESTIMATE_PARAMS) {                                                     //todo remove threshold, not used
            signalUtils_estimateNanoporeParams(sMt, npRead, &npRead->templateParams, ASSIGNMENT_THRESHOLD,
                                               signalUtils_templateOneDAssignmentsFromRead,
                                               nanopore_adjustTemplateEventsForDrift);
        }
        // make empty HMM to collect expectations
        Hmm *templateExpectations = hmmContinuous_getExpectationsHmm(sMt, p->threshold);

        // get expectations for template
        fprintf(stderr, "signalAlign - getting expectations for template\n");

        getSignalExpectations(sMt, templateModelFile, templateHmmFile, nHdpT,
                              templateExpectations, sMtype, npRead->templateParams,
                              tEventSequence, npRead->templateEventMap, pA->start2,
                              R->getTemplateTargetSequence(R), p,
                              anchorPairs, template, degenerate);

        if (sMtype == threeStateHdp) {
            fprintf(stderr, "signalAlign - got %" PRId64 "template HDP assignments\n",
                    hmmContinuous_howManyAssignments(templateExpectations));
        }

        // write to file
        fprintf(stderr, "signalAlign - writing expectations to file: %s\n", templateExpectationsFile);

        hmmContinuous_writeToFile(templateExpectationsFile, templateExpectations, sMtype);

        // get expectations for the complement
        fprintf(stderr, "signalAlign - getting expectations for complement\n");
        if (complementHmmFile == NULL) {
            st_errAbort("[signalMachine] ERROR: need to have input HMMs\n");
        }

        StateMachine *sMc = buildStateMachine(complementModelFile, npRead->complementParams, sMtype, nHdpC);

        if (ESTIMATE_PARAMS) {
            signalUtils_estimateNanoporeParams(sMc, npRead, &npRead->complementParams, ASSIGNMENT_THRESHOLD,
                                               signalUtils_complementOneDAssignmentsFromRead,
                                               nanopore_adjustComplementEventsForDrift);
        }

        Hmm *complementExpectations = hmmContinuous_getExpectationsHmm(sMc, p->threshold);

        getSignalExpectations(sMc, complementModelFile, complementHmmFile, nHdpC,
                              complementExpectations, sMtype,
                              npRead->complementParams, cEventSequence, npRead->complementEventMap, pA->start2,
                              R->getComplementTargetSequence(R), p, anchorPairs, complement, degenerate);

        if (sMtype == threeStateHdp) {
            fprintf(stderr, "signalAlign - got %"PRId64"complement HDP assignments\n",
                    hmmContinuous_howManyAssignments(complementExpectations));
        }

        // write to file
        fprintf(stderr, "signalAlign - writing expectations to file: %s\n", complementExpectationsFile);
        hmmContinuous_writeToFile(complementExpectationsFile, complementExpectations, sMtype);

        stateMachine_destruct(sMt);
        stateMachine_destruct(sMc);
        signalUtils_ReferenceSequenceDestruct(R);
        hmmContinuous_destruct(templateExpectations, sMtype);
        hmmContinuous_destruct(complementExpectations, sMtype);
        nanopore_nanoporeReadDestruct(npRead);
        sequence_destruct(tEventSequence);
        sequence_destruct(cEventSequence);
        pairwiseAlignmentBandingParameters_destruct(p);
        destructPairwiseAlignment(pA);
        stList_destruct(anchorPairs);
        return 0;
    } else {
        // Alignment Procedure //
        // Template alignment
        fprintf(stderr, "signalAlign - starting template alignment\n");

        // make template stateMachine
        StateMachine *sMt = buildStateMachine(templateModelFile, npRead->templateParams, sMtype, nHdpT);

        // re-estimate the nanoporeAdjustment parameters
        if (ESTIMATE_PARAMS) {
            signalUtils_estimateNanoporeParams(sMt, npRead, &npRead->templateParams, ASSIGNMENT_THRESHOLD,
                                               signalUtils_templateOneDAssignmentsFromRead,
                                               nanopore_adjustTemplateEventsForDrift);
        }

        stList *templateAlignedPairs = performSignalAlignment(sMt, tEventSequence, npRead->templateEventMap,
                                                              pA->start2, R->getTemplateTargetSequence(R),
                                                              p, anchorPairs,
                                                              degenerate);

        double templatePosteriorScore = scoreByPosteriorProbabilityIgnoringGaps(templateAlignedPairs);

        // sort
        stList_sort(templateAlignedPairs, sortByXPlusYCoordinate2); //Ensure the coordinates are increasing

        // write to file
        if (posteriorProbsFile != NULL) {
            if (sparseOutput) {
                writePosteriorProbsSparse(posteriorProbsFile, readLabel, R->getTemplateTargetSequence(R),
                                          forward, pA->contig1, tCoordinateShift,
                                          rCoordinateShift_t, templateAlignedPairs, template);
            } else {
                writePosteriorProbs(posteriorProbsFile, readLabel, sMt->EMISSION_MATCH_MATRIX,
                                    npRead->templateParams, npRead->templateEvents, R->getTemplateTargetSequence(R),
                                    forward, pA->contig1, sMt->type, tCoordinateShift, rCoordinateShift_t,
                                    templateAlignedPairs, template);
            }
        }

        // Complement alignment
        fprintf(stderr, "signalAlign - starting complement alignment\n");
        StateMachine *sMc = buildStateMachine(complementModelFile, npRead->complementParams, sMtype, nHdpC);

        if (ESTIMATE_PARAMS) {
            signalUtils_estimateNanoporeParams(sMc, npRead, &npRead->complementParams, ASSIGNMENT_THRESHOLD,
                                               signalUtils_complementOneDAssignmentsFromRead,
                                               nanopore_adjustComplementEventsForDrift);
        }

        stList *complementAlignedPairs = performSignalAlignment(sMc, cEventSequence,
                                                                npRead->complementEventMap, pA->start2,
                                                                R->getComplementTargetSequence(R),
                                                                p, anchorPairs, degenerate);

        double complementPosteriorScore = scoreByPosteriorProbabilityIgnoringGaps(complementAlignedPairs);

        // sort
        stList_sort(complementAlignedPairs, sortByXPlusYCoordinate2); //Ensure the coordinates are increasing

        // write to file
        if (posteriorProbsFile != NULL) {
            if (sparseOutput) {
                writePosteriorProbsSparse(posteriorProbsFile, readLabel, R->getComplementTargetSequence(R),
                                          forward, pA->contig1, cCoordinateShift,
                                          rCoordinateShift_c, complementAlignedPairs, complement);
            } else {
                writePosteriorProbs(posteriorProbsFile, readLabel, sMc->EMISSION_MATCH_MATRIX,
                                    npRead->complementParams, npRead->complementEvents,
                                    R->getComplementTargetSequence(R), forward, pA->contig1,
                                    sMc->type, cCoordinateShift, rCoordinateShift_c,
                                    complementAlignedPairs, complement);
            }
        }

        fprintf(stdout, "%s %"PRId64"\t%"PRId64"(%f)\t", readLabel, stList_length(anchorPairs),
                stList_length(templateAlignedPairs), templatePosteriorScore);
        fprintf(stdout, "%"PRId64"(%f)\n", stList_length(complementAlignedPairs), complementPosteriorScore);

        // final alignment clean up
        destructPairwiseAlignment(pA);
        nanopore_nanoporeReadDestruct(npRead);
        signalUtils_ReferenceSequenceDestruct(R);
        stateMachine_destruct(sMt);
        sequence_destruct(tEventSequence);
        stList_destruct(templateAlignedPairs);
        stateMachine_destruct(sMc);
        sequence_destruct(cEventSequence);
        stList_destruct(complementAlignedPairs);
        fprintf(stderr, "signalAlign - SUCCESS: finished alignment of query %s, exiting\n", readLabel);
    }
    return 0;
}

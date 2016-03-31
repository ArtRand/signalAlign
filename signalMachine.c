#include <getopt.h>
#include "pairwiseAlignment.h"
#include "pairwiseAligner.h"
#include "emissionMatrix.h"
#include "stateMachine.h"
#include "nanopore.h"
#include "continuousHmm.h"


void usage() {
    fprintf(stderr, "vanillaAlign binary, meant to be used through the signalAlign program.\n");
    fprintf(stderr, "See doc for signalAlign for help\n");
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
        int64_t x_adj;  // x is the reference coordinate that we record in the aligned pairs
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

        // get the kmer index
        int64_t targetKmerIndex = emissions_discrete_getKmerIndexFromKmer(k_i);

        // get the expected event mean amplitude and noise
        double E_mean = matchModel[1 + (targetKmerIndex * MODEL_PARAMS)];
        double E_noise = matchModel[1 + (targetKmerIndex * MODEL_PARAMS + 2)];
        double scaled_Emean = E_mean * npp.scale + npp.shift;
        double scaled_Enoise = E_noise * npp.scale_sd;
        double descaledEventMean = emissions_signal_descaleEventMean_JordanStyle(eventMean, E_mean,
                                                                                 npp.scale, npp.shift, npp.var);

        // make reference kmer
        char *refKmer;
        if ((strand == template && forward) || (strand == complement && (!forward))) {
            refKmer = k_i;
        }
        if ((strand == complement && forward) || (strand == template && (!forward))) {
            refKmer = stString_reverseComplementString(k_i);
        }
        // write to disk
        fprintf(fH, "%s\t%lld\t%s\t%s\t%s\t%lld\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%s\n",
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
        int64_t x_adj;  // x is the reference coordinate that we record in the aligned pairs
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
        char *refKmer;
        if ((strand == template && forward) || (strand == complement && (!forward))) {
            refKmer = k_i;
        }
        if ((strand == complement && forward) || (strand == template && (!forward))) {
            refKmer = stString_reverseComplementString(k_i);
        }
        fprintf(fH, "%s\t%lld\t%s\t%s\t%s\t%lld\t%s\t%f\t%s\n",
                contig, x_adj, refKmer, readFile, strandLabel, y, k_i, p, pathKmer);
        // cleanup
        free(k_i);
    }
    fclose(fH);
}

stList *getRemappedAnchorPairs(stList *unmappedAnchors, int64_t *eventMap, int64_t mapOffset) {
    stList *remapedAnchors = nanopore_remapAnchorPairsWithOffset(unmappedAnchors, eventMap, mapOffset);
    stList *filteredRemappedAnchors = filterToRemoveOverlap(remapedAnchors);
    return filteredRemappedAnchors;
}

StateMachine *buildStateMachine(const char *modelFile, NanoporeReadAdjustmentParameters npp, StateMachineType type,
                                Strand strand, NanoporeHDP *nHdp) {
    if ((type != threeState) && (type != threeStateHdp)) {
        st_errAbort("signalAlign - incompatable stateMachine type request");
    }

    if (type == threeState) {
        StateMachine *sM = getSM3_descaled(modelFile, npp);
        //StateMachine *sM = getStrawManStateMachine3(modelFile);
        //emissions_signal_scaleModel(sM, npp.scale, npp.shift, npp.var, npp.scale_sd, npp.var_sd);
        return sM;
    }
    if (type == threeStateHdp) {
        StateMachine *sM = getHdpMachine(nHdp, modelFile, npp);
        //StateMachine *sM = getHdpStateMachine3(nHdp, modelFile);
        return sM;
    }
    else {
        st_errAbort("signalAlign - ERROR: buildStateMachine, didn't get correct input\n");
    }
    return 0;
}

void updateHdpFromAssignments(const char *nHdpFile, const char *expectationsFile, const char *nHdpOutFile) {
    NanoporeHDP *nHdp = deserialize_nhdp(nHdpFile);
    Hmm *hdpHmm = hdpHmm_loadFromFile(expectationsFile, nHdp);
    hmmContinuous_destruct(hdpHmm, hdpHmm->type);

    fprintf(stderr, "signalAlign - Running Gibbs on HDP\n");
    execute_nhdp_gibbs_sampling(nHdp, 10000, 100000, 100, FALSE);
    finalize_nhdp_distributions(nHdp);

    fprintf(stderr, "signalAlign - Serializing HDP to %s\n", nHdpOutFile);
    serialize_nhdp(nHdp, nHdpOutFile);
    destroy_nanopore_hdp(nHdp);
}

void loadHmmRoutine(const char *hmmFile, StateMachine *sM, StateMachineType type, Hmm *expectations) {
    if ((type != threeState) && (type != threeStateHdp)) {
        st_errAbort("LoadSignalHmm : unupported stateMachineType");
    }
    hmmContinuous_loadSignalHmm(hmmFile, sM, type, expectations);
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

stList *performSignalAlignmentP(StateMachine *sM, Sequence *sY, int64_t *eventMap, int64_t mapOffset, char *target,
                                PairwiseAlignmentParameters *p, stList *unmappedAnchors,
                                void *(*targetGetFcn)(void *, int64_t),
                                void (*posteriorProbFcn)(StateMachine *sM, int64_t xay, DpMatrix *forwardDpMatrix,
                                                         DpMatrix *backwardDpMatrix, Sequence* sX, Sequence* sY,
                                                         double totalProbability, PairwiseAlignmentParameters *p,
                                                         void *extraArgs),
                                bool twoWay) {
    int64_t lX = sequence_correctSeqLength(strlen(target), kmer);

    // remap anchor pairs
    stList *filteredRemappedAnchors = getRemappedAnchorPairs(unmappedAnchors, eventMap, mapOffset);

    // make sequences
    Sequence *sX = sequence_constructReferenceSequence(lX, target, targetGetFcn, sequence_sliceNucleotideSequence2,
                                                       (twoWay ? TWO_CYTOSINES : THREE_CYTOSINES),
                                                       (twoWay ? TWO_WAY_CYTOSINE_OPTIONS : NB_CYTOSINE_OPTIONS),
                                                       kmer);

    // do alignment
    stList *alignedPairs = getAlignedPairsUsingAnchors(sM, sX, sY, filteredRemappedAnchors, p,
                                                       posteriorProbFcn, 1, 1);
    return alignedPairs;
}

stList *performSignalAlignment(StateMachine *sM, Sequence *eventSequence, int64_t *eventMap,
                               int64_t mapOffset, char *target, PairwiseAlignmentParameters *p,
                               stList *unmappedAncors, bool twoWay) {
    if ((sM->type != threeState) && (sM->type != threeStateHdp)) {
        st_errAbort("signalAlign - You're trying to do the wrong king of alignment");
    }

    stList *alignedPairs = performSignalAlignmentP(sM, eventSequence, eventMap, mapOffset, target, p,
                                                   unmappedAncors, sequence_getKmer3,
                                                   diagonalCalculationPosteriorMatchProbs, twoWay);
    return alignedPairs;
}

char *getSubSequence(char *seq, int64_t start, int64_t end, bool strand) {
    if (strand) {
        seq = stString_getSubString(seq, start, end - start);
        return seq;
    }
    seq = stString_getSubString(seq, end, start - end);
    return seq;
}

void rebasePairwiseAlignmentCoordinates(int64_t *start, int64_t *end, int64_t *strand,
                                        int64_t coordinateShift, bool flipStrand) {
    *start += coordinateShift;
    *end += coordinateShift;
    if (flipStrand) {
        *strand = *strand ? 0 : 1;
        int64_t i = *end;
        *end = *start;
        *start = i;
    }
}

stList *guideAlignmentToRebasedAnchorPairs(struct PairwiseAlignment *pA, PairwiseAlignmentParameters *p) {
    // check if we need to flip the reference
    bool flipStrand1 = !pA->strand1;
    int64_t refCoordShift = (pA->strand1 ? pA->start1 : pA->end1);

    // rebase the reference alignment to (0), but not the nanopore read, this is corrected when remapping the
    // anchorPairs
    rebasePairwiseAlignmentCoordinates(&(pA->start1), &(pA->end1), &(pA->strand1), -refCoordShift, flipStrand1);
    checkPairwiseAlignment(pA);

    //Convert input alignment into anchor pairs
    stList *unfilteredAnchorPairs = convertPairwiseForwardStrandAlignmentToAnchorPairs(
            pA, p->constraintDiagonalTrim);

    // sort
    stList_sort(unfilteredAnchorPairs, (int (*)(const void *, const void *)) stIntTuple_cmpFn);

    // filter
    stList *anchorPairs = filterToRemoveOverlap(unfilteredAnchorPairs);

    return anchorPairs;
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
    Sequence *eventS = sequence_construct2(endIdx - startIdx, elements, sequence_getEvent,
                                           sequence_sliceEventSequence2, event);

    return eventS;
}

void getSignalExpectations(const char *model, const char *inputHmm, NanoporeHDP *nHdp,
                           Hmm *hmmExpectations, StateMachineType type,
                           NanoporeReadAdjustmentParameters npp, Sequence *eventSequence,
                           int64_t *eventMap, int64_t mapOffset, char *trainingTarget, PairwiseAlignmentParameters *p,
                           stList *unmappedAnchors, Strand strand, bool twoWay) {
    // load match model, build stateMachine
    StateMachine *sM = buildStateMachine(model, npp, type, strand, nHdp);

    // load HMM
    if (inputHmm == NULL) {
        st_errAbort("[signalMachine] ERROR: need to have input HMMs\n");
    }

    fprintf(stderr, "signalAlign - loading HMM from file, %s\n", inputHmm);
    loadHmmRoutine(inputHmm, sM, type, hmmExpectations);

    // correct sequence length
    int64_t lX = sequence_correctSeqLength(strlen(trainingTarget), event);

    // remap the anchors
    stList *filteredRemappedAnchors = getRemappedAnchorPairs(unmappedAnchors, eventMap, mapOffset);

    Sequence *target = sequence_constructReferenceSequence(lX, trainingTarget,
                                                           sequence_getKmer3, sequence_sliceNucleotideSequence2,
                                                           (twoWay ? TWO_CYTOSINES : THREE_CYTOSINES),
                                                           (twoWay ? TWO_WAY_CYTOSINE_OPTIONS : NB_CYTOSINE_OPTIONS),
                                                           kmer);

    getExpectationsUsingAnchors(sM, hmmExpectations, target, eventSequence, filteredRemappedAnchors, p,
                                diagonalCalculation_Expectations, 1, 1);

    stateMachine_destruct(sM);
}

int main(int argc, char *argv[]) {
    StateMachineType sMtype = threeState;
    // HDP stuff
    int64_t j = 0;
    int64_t diagExpansion = 50;
    double threshold = 0.01;
    int64_t constraintTrim = 14;
    char *templateModelFile = stString_print("../../signalAlign/models/testModel_template.model");
    char *complementModelFile = stString_print("../../signalAlign/models/testModel_complement.model");
    char *readLabel = NULL;
    char *npReadFile = NULL;
    char *targetFile = NULL;
    char *posteriorProbsFile = NULL;
    char *templateHmmFile = NULL;
    char *complementHmmFile = NULL;
    char *templateExpectationsFile = NULL;
    char *complementExpectationsFile = NULL;
    char *templateHdp = NULL;
    char *complementHdp = NULL;
    char *cytosine_substitute = NULL;
    bool sparseOutput = FALSE;
    bool twoWay = FALSE;

    int key;
    while (1) {
        static struct option long_options[] = {
                {"help",                    no_argument,        0,  'h'},
                {"sm3Hdp",                  no_argument,        0,  'd'},
                {"twoWay",                  no_argument,        0,  'o'},
                {"sparse_output",           no_argument,        0,  's'},
                {"substitute",              required_argument,  0,  'M'},
                {"templateModel",           required_argument,  0,  'T'},
                {"complementModel",         required_argument,  0,  'C'},
                {"readLabel",               required_argument,  0,  'L'},
                {"npRead",                  required_argument,  0,  'q'},
                {"reference",               required_argument,  0,  'r'},
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

        key = getopt_long(argc, argv, "h:osd:p:M:a:T:C:L:q:r:u:y:z:v:w:t:c:x:D:m:",
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
                twoWay = TRUE;
                break;
            case 'd':
                sMtype = threeStateHdp;
                break;
            case 'M':
                cytosine_substitute = stString_copy(optarg);
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
            case 'r':
                targetFile = stString_copy(optarg);
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

    if (sMtype == threeState) {
        fprintf(stderr, "signalAlign - using three-state HMM model\n");
    }
    if (sMtype == threeStateHdp) {
        fprintf(stderr, "signalAlign - using three-state HMM-HDP model\n");
    }

    NanoporeHDP *nHdpT, *nHdpC;

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

    if ((templateHdp != NULL) || (complementHdp != NULL)) {
        if ((templateHdp == NULL) || (complementHdp == NULL)) {
            st_errAbort("Need to have template and complement HDPs");
        }
        if (sMtype != threeStateHdp) {
            fprintf(stderr, "signalAlign - Warning: this kind of stateMachine does not use the HDPs you gave\n");
        }
        fprintf(stderr, "signalAlign - using NanoporeHDPs\n");
    }

    // load reference sequence (reference sequence)
    FILE *reference = fopen(targetFile, "r");
    char *referenceSequence = stFile_getLineFromFile(reference);

    // load nanopore read
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);

    // make some params
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
    p->threshold = threshold;
    p->constraintDiagonalTrim = constraintTrim;
    p->diagonalExpansion = diagExpansion;
    // get pairwise alignment from stdin, in exonerate CIGAR format
    FILE *fileHandleIn = stdin;

    // parse input
    struct PairwiseAlignment *pA;
    pA = cigarRead(fileHandleIn);

    // put in to help with debugdebugging
    // printPairwiseAlignmentSummary(pA);

    // slice out the section of the reference we're aligning to
    char *trimmedRefSeq = getSubSequence(referenceSequence, pA->start1, pA->end1, pA->strand1);
    trimmedRefSeq = (pA->strand1 ? trimmedRefSeq : stString_reverseComplementString(trimmedRefSeq));

    // reverse complement for complement event sequence
    char *rc_trimmedRefSeq = stString_reverseComplementString(trimmedRefSeq);

    // change bases to methylated/hydroxymethylated, if asked to
    char *templateTargetSeq = cytosine_substitute == NULL ? trimmedRefSeq : stString_replace(trimmedRefSeq, "C",
                                                                                             cytosine_substitute);
    char *complementTargetSeq = cytosine_substitute == NULL ? rc_trimmedRefSeq : stString_replace(rc_trimmedRefSeq,
                                                                                                  "C",
                                                                                                  cytosine_substitute);

    // constrain the event sequence to the positions given by the guide alignment
    Sequence *tEventSequence = makeEventSequenceFromPairwiseAlignment(npRead->templateEvents,
                                                                      pA->start2, pA->end2,
                                                                      npRead->templateEventMap);

    Sequence *cEventSequence = makeEventSequenceFromPairwiseAlignment(npRead->complementEvents,
                                                                      pA->start2, pA->end2,
                                                                      npRead->complementEventMap);


    // the aligned pairs start at (0,0) so we need to correct them based on the guide alignment later.
    // record the pre-zeroed alignment start and end coordinates here

    // for the events:
    int64_t tCoordinateShift = npRead->templateEventMap[pA->start2];
    int64_t cCoordinateShift = npRead->complementEventMap[pA->start2];

    // and for the reference:
    int64_t rCoordinateShift_t = pA->start1;
    int64_t rCoordinateShift_c = pA->end1;
    bool forward = pA->strand1;  // keep track of whether this is a forward mapped read or not

    stList *anchorPairs = guideAlignmentToRebasedAnchorPairs(pA, p);

    if ((templateExpectationsFile != NULL) && (complementExpectationsFile != NULL)) {
        // Expectation Routine //
        // make empty HMM to collect expectations
        Hmm *templateExpectations = hmmContinuous_getEmptyHmm(sMtype, 0.0001, p->threshold,
                                                              npRead->templateParams.scale,
                                                              npRead->templateParams.shift,
                                                              npRead->templateParams.var);
        Hmm *complementExpectations = hmmContinuous_getEmptyHmm(sMtype, 0.0001, p->threshold,
                                                                npRead->complementParams.scale,
                                                                npRead->complementParams.shift,
                                                                npRead->complementParams.var);

        // get expectations for template
        fprintf(stderr, "signalAlign - getting expectations for template\n");
        getSignalExpectations(templateModelFile, templateHmmFile, nHdpT,
                              templateExpectations, sMtype, npRead->templateParams,
                              tEventSequence, npRead->templateEventMap, pA->start2, templateTargetSeq, p,
                              anchorPairs, template, twoWay);

        if (sMtype == threeStateHdp) {
            fprintf(stderr, "signalAlign - got %lld template HDP assignments\n",
                    hmmContinuous_howManyAssignments(templateExpectations));
        }

        // write to file
        fprintf(stderr, "signalAlign - writing expectations to file: %s\n", templateExpectationsFile);

        hmmContinuous_writeToFile(templateExpectationsFile, templateExpectations, sMtype);

        // get expectations for the complement
        fprintf(stderr, "signalAlign - getting expectations for complement\n");
        getSignalExpectations(complementModelFile, complementHmmFile, nHdpC,
                              complementExpectations, sMtype,
                              npRead->complementParams, cEventSequence, npRead->complementEventMap, pA->start2,
                              complementTargetSeq, p, anchorPairs, complement, twoWay);

        if (sMtype == threeStateHdp) {
            fprintf(stderr, "signalAlign - got %lld complement HDP assignments\n",
                    hmmContinuous_howManyAssignments(complementExpectations));
        }

        // write to file
        fprintf(stderr, "signalAlign - writing expectations to file: %s\n", complementExpectationsFile);
        hmmContinuous_writeToFile(complementExpectationsFile, complementExpectations, sMtype);

        hmmContinuous_destruct(templateExpectations, sMtype);
        hmmContinuous_destruct(complementExpectations, sMtype);
        nanopore_nanoporeReadDestruct(npRead);
        sequence_sequenceDestroy(tEventSequence);
        sequence_sequenceDestroy(cEventSequence);
        pairwiseAlignmentBandingParameters_destruct(p);
        destructPairwiseAlignment(pA);
        stList_destruct(anchorPairs);
        if (nHdpT != NULL) {
            destroy_nanopore_hdp(nHdpT);
        }
        if (nHdpC != NULL) {
            destroy_nanopore_hdp(nHdpC);
        }

        return 0;
    } else {
        // Alignment Procedure //
        StateMachine *sMt;
        stList *templateAlignedPairs, *complementAlignedPairs;
        double templatePosteriorScore, complementPosteriorScore;

        // Template alignment
        fprintf(stderr, "signalAlign - starting template alignment\n");

        // make template stateMachine
        sMt = buildStateMachine(templateModelFile, npRead->templateParams, sMtype, template, nHdpT);
        // load HMM, if given
        if (templateHmmFile != NULL) {
            fprintf(stderr, "loading HMM from file, %s\n", templateHmmFile);
            loadHmmRoutine(templateHmmFile, sMt, sMt->type, NULL);
        }
        // get aligned pairs
        templateAlignedPairs = performSignalAlignment(sMt, tEventSequence, npRead->templateEventMap,
                                                      pA->start2, templateTargetSeq, p, anchorPairs, twoWay);

        templatePosteriorScore = scoreByPosteriorProbabilityIgnoringGaps(templateAlignedPairs);

        // sort
        stList_sort(templateAlignedPairs, sortByXPlusYCoordinate2); //Ensure the coordinates are increasing

        // write to file
        if (posteriorProbsFile != NULL) {
            if (sparseOutput) {
                writePosteriorProbsSparse(posteriorProbsFile, readLabel, trimmedRefSeq, forward, pA->contig1,
                                          tCoordinateShift, rCoordinateShift_t, templateAlignedPairs, template);
            } else {
                writePosteriorProbs(posteriorProbsFile, readLabel, sMt->EMISSION_MATCH_MATRIX,
                                    npRead->templateParams, npRead->templateEvents, trimmedRefSeq, forward,
                                    pA->contig1, sMt->type, tCoordinateShift, rCoordinateShift_t,
                                    templateAlignedPairs, template);
            }
        }

        // Complement alignment
        fprintf(stderr, "signalAlign - starting complement alignment\n");
        StateMachine *sMc = buildStateMachine(complementModelFile, npRead->complementParams, sMtype, complement, nHdpC);
        if (complementHmmFile != NULL) {
            fprintf(stderr, "loading HMM from file, %s\n", complementHmmFile);
            loadHmmRoutine(complementHmmFile, sMc, sMc->type, NULL);
        }
        // get aligned pairs
        complementAlignedPairs = performSignalAlignment(sMc, cEventSequence,
                                                        npRead->complementEventMap, pA->start2,
                                                        complementTargetSeq, p, anchorPairs, twoWay);

        complementPosteriorScore = scoreByPosteriorProbabilityIgnoringGaps(complementAlignedPairs);

        // sort
        stList_sort(complementAlignedPairs, sortByXPlusYCoordinate2); //Ensure the coordinates are increasing

        // write to file
        if (posteriorProbsFile != NULL) {
            if (sparseOutput) {
                writePosteriorProbsSparse(posteriorProbsFile, readLabel, rc_trimmedRefSeq, forward, pA->contig1,
                                          cCoordinateShift, rCoordinateShift_c, complementAlignedPairs, complement);
            } else {
                writePosteriorProbs(posteriorProbsFile, readLabel, sMc->EMISSION_MATCH_MATRIX,
                                    npRead->complementParams, npRead->complementEvents, rc_trimmedRefSeq,
                                    forward, pA->contig1, sMc->type,cCoordinateShift, rCoordinateShift_c,
                                    complementAlignedPairs, complement);
            }
        }

        fprintf(stdout, "%s %lld\t%lld(%f)\t", readLabel, stList_length(anchorPairs),
                stList_length(templateAlignedPairs), templatePosteriorScore);
        fprintf(stdout, "%lld(%f)\n", stList_length(complementAlignedPairs), complementPosteriorScore);

        // final alignment clean up
        stateMachine_destruct(sMt);
        sequence_sequenceDestroy(tEventSequence);
        stList_destruct(templateAlignedPairs);
        stateMachine_destruct(sMc);
        sequence_sequenceDestroy(cEventSequence);
        stList_destruct(complementAlignedPairs);
        fprintf(stderr, "signalAlign - SUCCESS: finished alignment of query %s, exiting\n", readLabel);
    }
    return 0;
}

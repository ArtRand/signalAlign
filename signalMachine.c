#include <getopt.h>
#include <string.h>
#include "signalMachineUtils.h"
#include "pairwiseAligner.h"

#define ESTIMATE_PARAMS 1
#define ASSIGNMENT_THRESHOLD 0.1


typedef enum {
    full = 0,
    variantCaller = 1,
    assignments = 2
} OutputFormat;

void usage() {
    fprintf(stderr, "\n\tsignalMachine - Align ONT ionic current to a reference sequence\n\n");
    fprintf(stderr, "--help: Display this super useful message and exit\n");
    fprintf(stderr, "--sm3Hdp, -d: Flag, enable HMM-HDP model\n");
    fprintf(stderr, "--twoD, -e: Flag, use 2D workflow (enables complement alignment)\n");
    fprintf(stderr, "-s: Output format, 0=full, 1=variantCaller, 2=assignments\n");
    fprintf(stderr, "-o: Degernate, 0=C/E, 1=C/E/O, 2=A/I, 3=A/C/G/T\n");
    fprintf(stderr, "-T: Template HMM model\n");
    fprintf(stderr, "-C: Complement HMM model\n");
    fprintf(stderr, "-L: Read (output) label\n");
    fprintf(stderr, "-q: NanoporeRead (in npRead format)\n");
    fprintf(stderr, "-f: Forward reference to align to as a flat file\n");
    fprintf(stderr, "-b: Backward reference to align to as a flat fiel\n");
    fprintf(stderr, "-p: Guide alignment file, containing CIGARs in EXONERATE format\n");
    fprintf(stderr, "-u: Posteriors (output) file path, place to put the output\n");
    fprintf(stderr, "-v: TemplateHDP file\n");
    fprintf(stderr, "-w: Complement HDP file\n");
    fprintf(stderr, "-t: Template expectations (HMM transitions) output location\n");
    fprintf(stderr, "-c: Complement expectations (HMM transitions) output location\n");
    fprintf(stderr, "-x: Diagonal expansion, how much to expand the dynamic programming envelope\n");
    fprintf(stderr, "-D: Posterior probability threshold, keep aligned pairs with posterior prob >= this\n");
    fprintf(stderr, "-m: Constranint trim, how much to trim the guide alignment anchors by\n\n");
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

inline int64_t adjustReferenceCoordinate(int64_t x_i, int64_t referenceSeqOffset,
                                         int64_t referenceLengthInKmers, int64_t referenceLength,
                                         Strand strand, bool forward) {
    if ((strand == template && forward) || (strand == complement && !forward)) {
        return x_i + referenceSeqOffset;
    } else {
        return referenceLengthInKmers - (x_i + (referenceLength - referenceSeqOffset));
    }
}

inline char *makeReferenceKmer(const char *k_i, Strand strand, bool forward) {
    if ((strand == template && forward) || (strand == complement && !forward)) {
        return stString_copy(k_i);
    } else {
        return stString_reverseComplementString(k_i);
    }
}

inline char *kmerFromString(const char *string, int64_t start, int64_t kmerLength) {
    char *k_i = st_malloc(kmerLength * sizeof(char));
    for (int64_t i = 0; i < kmerLength; i++) {
        k_i[i] = *(string + (start + i));
    }
    k_i[kmerLength] = '\0';
    return k_i;
}

inline int64_t adjustQueryPosition(int64_t unadjustedQueryPosition, int64_t kmerLength, Strand strand, bool forward) {
    if ((strand == template && forward) || (strand == complement && !forward)) {
        return unadjustedQueryPosition;
    } else {
        return (kmerLength - 1) - unadjustedQueryPosition;
    }
}

void writePosteriorProbsFull(char *posteriorProbsFile, char *readLabel, StateMachine *sM,
                             NanoporeReadAdjustmentParameters npp, double *events, char *target, bool forward,
                             char *contig, int64_t eventSequenceOffset, int64_t referenceSequenceOffset,
                             stList *alignedPairs, Strand strand) {
    // label for tsv output
    char *strandLabel = strand == template ? "t" : "c";

    // open the file for output
    FILE *fH = fopen(posteriorProbsFile, "a");

    // get some lengths outside the loop
    int64_t refLength = (int64_t )strlen(target);
    int64_t refLengthInKmers = refLength - sM->kmerLength;

    for(int64_t i = 0; i < stList_length(alignedPairs); i++) {
        // grab the aligned pair
        stIntTuple *aPair = stList_get(alignedPairs, i);

        if (stIntTuple_length(aPair) != 4) {
            st_errAbort("Aligned pair tuples should have length 4, this one has length %lld\n",
                        stIntTuple_length(aPair));
        }

        // nucleotide sequence coordinate
        int64_t x_i = stIntTuple_get(aPair, 1);
        // adjust back to reference coordinates
        int64_t x_adj = adjustReferenceCoordinate(x_i, referenceSequenceOffset, refLengthInKmers, refLength,
                                                  strand, forward);
        // event index, adjust to to entire event sequence coordinates (event sequence is trimmed during alignment)
        int64_t y = stIntTuple_get(aPair, 2) + eventSequenceOffset;
        // posterior probability
        double p = ((double)stIntTuple_get(aPair, 0)) / PAIR_ALIGNMENT_PROB_1;
        // path (variant-called) kmer
        char *pathKmer = (char *)stIntTuple_get(aPair, 3);

        double eventMean = sequence_getEventMean(events, y);
        double eventNoise = sequence_getEventNoise(events, y);
        double eventDuration = sequence_getEventDuration(events, y);

        // make the kmer string at the target index,
        char *k_i = kmerFromString(target, x_i, sM->kmerLength);

        int64_t targetKmerIndex = kmer_id(pathKmer, sM->alphabet, sM->alphabetSize, sM->kmerLength);

        // get the expected event mean amplitude and noise
        double E_mean = sM->EMISSION_MATCH_MATRIX[(targetKmerIndex * MODEL_PARAMS)];
        double E_noise = sM->EMISSION_MATCH_MATRIX[(targetKmerIndex * MODEL_PARAMS + 2)];
        double scaled_Emean = E_mean * npp.scale + npp.shift;
        double scaled_Enoise = E_noise * npp.scale_sd;
        double descaledEventMean = emissions_signal_descaleEventMean_JordanStyle(eventMean, E_mean,
                                                                                 npp.scale, npp.shift, npp.var);

        // make reference kmer
        char *refKmer = makeReferenceKmer(k_i, strand, forward);

        // write to file
        fprintf(fH, "%s\t%"PRId64"\t%s\t%s\t%s\t%"PRId64"\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%s\n",
                contig, x_adj, refKmer, readLabel, strandLabel, y, eventMean, eventNoise, eventDuration, k_i,
                scaled_Emean, scaled_Enoise, p, descaledEventMean, E_mean, pathKmer);

        // cleanup
        free(k_i);
        free(refKmer);
    }
    fclose(fH);
}

void writePosteriorProbsVC(char *posteriorProbsFile, char *readLabel, StateMachine *sM, char *target, bool forward,
                           int64_t eventSequenceOffset, int64_t referenceSequenceOffset, stList *alignedPairs,
                           Strand strand, double posteriorScore) {
    // label for tsv output
    char *strandLabel = strand == template ? "t" : "c";
    char *forwardLabel = forward ? "forward" : "backward";

    // open the file for output
    FILE *fH = fopen(posteriorProbsFile, "a");

    // get some lengths outside the loop
    int64_t refLength = (int64_t )strlen(target);
    int64_t refLengthInKmers = refLength - sM->kmerLength;

    for(int64_t i = 0; i < stList_length(alignedPairs); i++) {
        // grab the aligned pair
        stIntTuple *aPair = stList_get(alignedPairs, i);

        if (stIntTuple_length(aPair) != 4) {
            st_errAbort("Aligned pair tuples should have length 4, this one has length %lld\n",
                        stIntTuple_length(aPair));
        }

        // trimmed nucleotide sequence coordinate
        int64_t x_i = stIntTuple_get(aPair, 1);
        // make the kmer string at the target index,
        char *k_i = kmerFromString(target, x_i, sM->kmerLength);
        char *refKmer = makeReferenceKmer(k_i, strand, forward);
        stList *queryPositions = path_findDegeneratePositions(refKmer, sM->kmerLength);
        // check if this aligned pair reports on a query position
        if (stList_length(queryPositions) == 0) {
            free(k_i);
            free(refKmer);
            stList_destruct(queryPositions);
            continue;
        }
        // adjust back to reference coordinates
        int64_t x_adj = adjustReferenceCoordinate(x_i, referenceSequenceOffset, refLengthInKmers, refLength,
                                                  strand, forward);
        // event index, adjust to to entire event sequence coordinates (event sequence is trimmed during alignment)
        int64_t y = stIntTuple_get(aPair, 2) + eventSequenceOffset;
        // posterior probability
        double p = ((double)stIntTuple_get(aPair, 0)) / PAIR_ALIGNMENT_PROB_1;
        // path (variant-called) kmer
        char *pathKmer = (char *)stIntTuple_get(aPair, 3);
        // get the base that was called in this aligned pair
        int64_t nQueryPositions = stList_length(queryPositions);
        for (int64_t q = 0; q < nQueryPositions; q++) {
            // position in the reference kmer eg. AGXGG -> 2
            int64_t unadjustedQueryPosition = *(int64_t  *)stList_get(queryPositions, q);
            // position in the pathKmer
            int64_t queryPosition = adjustQueryPosition(unadjustedQueryPosition, sM->kmerLength,
                                                        strand, forward);
            // called base
            char base = pathKmer[queryPosition];
            // position in the reference we're reporting on
            int64_t reportPosition = x_adj + unadjustedQueryPosition;
            fprintf(fH, "%"PRId64"\t%"PRId64"\t%c\t%f\t%s\t%s\t%s\t%f\n", y, reportPosition, base, p,
                    strandLabel, forwardLabel, readLabel, posteriorScore);
        }
        free(k_i);
        free(refKmer);
        stList_destruct(queryPositions);
    }
    fclose(fH);
}

void writeAssignments(char *posteriorProbsFile, StateMachine *sM, double *events, int64_t eventSequenceOffset,
                      NanoporeReadAdjustmentParameters npp, stList *alignedPairs, Strand strand) {
    // label for tsv output
    char *strandLabel = strand == template ? "t" : "c";

    // open the file for output
    FILE *fH = fopen(posteriorProbsFile, "a");

    for(int64_t i = 0; i < stList_length(alignedPairs); i++) {
        // grab the aligned pair
        stIntTuple *aPair = stList_get(alignedPairs, i);

        if (stIntTuple_length(aPair) != 4) {
            st_errAbort("Aligned pair tuples should have length 4, this one has length %lld\n",
                        stIntTuple_length(aPair));
        }

        // event index, adjust to to entire event sequence coordinates (event sequence is trimmed during alignment)
        int64_t y = stIntTuple_get(aPair, 2) + eventSequenceOffset;
        // posterior probability
        double p = ((double)stIntTuple_get(aPair, 0)) / PAIR_ALIGNMENT_PROB_1;
        // path (variant-called) kmer
        char *pathKmer = (char *)stIntTuple_get(aPair, 3);

        // get the observed event mean
        double eventMean = sequence_getEventMean(events, y);
        // get the kmer index
        int64_t targetKmerIndex = kmer_id(pathKmer, sM->alphabet, sM->alphabetSize, sM->kmerLength);
        // get the expected mean from the model
        double E_mean = sM->EMISSION_MATCH_MATRIX[(targetKmerIndex * MODEL_PARAMS)];
        // descale the observed mean
        double descaledEventMean = emissions_signal_descaleEventMean_JordanStyle(eventMean, E_mean,
                                                                                 npp.scale, npp.shift, npp.var);
        fprintf(fH, "%s\t%s\t%lf\t%lf\n", pathKmer, strandLabel, descaledEventMean, p);
    }
    fclose(fH);
}

void outputAlignment(
        OutputFormat fmt,
        char *posteriorProbsFile,
        char *readLabel,
        StateMachine *sM,
        NanoporeReadAdjustmentParameters npp,
        double *events,
        char *target,
        bool forward,
        char *contig,
        int64_t eventSequenceOffset,
        int64_t referenceSequenceOffset,
        stList *alignedPairs, 
        double posteriorScore,
        Strand strand) {
    switch (fmt) {
        case full:
            writePosteriorProbsFull(posteriorProbsFile, readLabel, sM, npp, events, target, forward, contig,
                                    eventSequenceOffset, referenceSequenceOffset, alignedPairs, strand);
            break;
        case variantCaller:
            writePosteriorProbsVC(posteriorProbsFile, readLabel, sM, target, forward, eventSequenceOffset,
                                  referenceSequenceOffset, alignedPairs, strand, posteriorScore);
            break;
        case assignments:
            writeAssignments(posteriorProbsFile, sM, events, eventSequenceOffset, npp, alignedPairs, strand);
            break;
        default:
            fprintf(stderr, "signalAlign - No valid output format provided\n");
            return;
    }
}

StateMachine *buildStateMachine(const char *modelFile, NanoporeReadAdjustmentParameters npp, StateMachineType type,
                                NanoporeHDP *nHdp) {
    if ((type != threeState) && (type != threeStateHdp)) {
        st_errAbort("signalAlign - incompatible stateMachine type request");
    }
    if (!stFile_exists(modelFile)) {
        st_errAbort("signalAlign - ERROR: couldn't find model file here: %s\n", modelFile);
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
    Sequence *sX = sequence_constructReferenceKmerSequence(lX, target, sequence_getKmer,
                                                           sequence_sliceNucleotideSequence, degenerate, kmer);

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

void getSignalExpectations(StateMachine *sM, Hmm *hmmExpectations, Sequence *eventSequence,
                           int64_t *eventMap, int64_t mapOffset, char *trainingTarget, PairwiseAlignmentParameters *p,
                           stList *unmappedAnchors, DegenerateType degenerate) {
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

}

int main(int argc, char *argv[]) {
    StateMachineType sMtype = threeState;
    int64_t j = 0;
    int64_t diagExpansion = 50;
    double threshold = 0.01;
    int64_t constraintTrim = 14;
    int64_t degenerate;
    int64_t outFmt;
    bool twoD = FALSE;
    char *templateModelFile = NULL;
    char *complementModelFile = NULL;
    char *readLabel = NULL;
    char *npReadFile = NULL;
    char *forwardReference = NULL;
    char *backwardReference = NULL;
    char *exonerateCigarFile= NULL;
    char *posteriorProbsFile = NULL;
    char *templateExpectationsFile = NULL;
    char *complementExpectationsFile = NULL;
    char *templateHdp = NULL;
    char *complementHdp = NULL;

    int key;
    while (1) {
        static struct option long_options[] = {
                {"help",                    no_argument,        0,  'h'},
                {"sm3Hdp",                  no_argument,        0,  'd'},
                {"sparse_output",           no_argument,        0,  's'},
                {"twoD",                    no_argument,        0,  'e'},
                {"degenerate",              required_argument,  0,  'o'},
                {"templateModel",           required_argument,  0,  'T'},
                {"complementModel",         required_argument,  0,  'C'},
                {"readLabel",               required_argument,  0,  'L'},
                {"npRead",                  required_argument,  0,  'q'},
                {"forward_reference",       required_argument,  0,  'f'},
                {"backward_reference",      required_argument,  0,  'b'},
                {"exonerate_cigar_file",    required_argument,  0,  'p'},
                {"posteriors",              required_argument,  0,  'u'},
                {"templateHdp",             required_argument,  0,  'v'},
                {"complementHdp",           required_argument,  0,  'w'},
                {"templateExpectations",    required_argument,  0,  't'},
                {"complementExpectations",  required_argument,  0,  'c'},
                {"diagonalExpansion",       required_argument,  0,  'x'},
                {"threshold",               required_argument,  0,  'D'},
                {"constraintTrim",          required_argument,  0,  'm'},
                {0, 0, 0, 0} };

        int option_index = 0;

        key = getopt_long(argc, argv, "h:d:e:s:o:a:T:C:L:q:f:b:p:u:v:w:t:c:x:D:m:",
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
                j = sscanf(optarg, "%" PRIi64 "", &outFmt);
                assert (j == 1);
                break;
            case 'e':
                twoD = TRUE;
                break;
            case 'o':
                j = sscanf(optarg, "%" PRIi64 "", &degenerate);
                assert (j == 1);
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
                exonerateCigarFile = stString_copy(optarg);
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
    
    // check for models
    if ((templateModelFile == NULL) || (complementModelFile == NULL && twoD)) {
        st_errAbort("Missing model files, exiting\n");
        return 1;
    }

    if (exonerateCigarFile == NULL) {
        st_errAbort("[signalMachine]ERROR: Need to provide input guide alignments, exiting\n");
        return 1;
    }

    // Anchors //
    // get pairwise alignment from stdin, in exonerate CIGAR format
    //FILE *fileHandleIn = stdin;
    if (!stFile_exists(exonerateCigarFile)) {
        st_errAbort("[signalMachine]ERROR: Didn't find input alignment file, looked %s\n", exonerateCigarFile);
    } else {
        st_uglyf("[signalMachine]NOTICE: Using guide alignments from %s\n", exonerateCigarFile);
    }
    
    FILE *fileHandleIn = fopen(exonerateCigarFile, "r");
    // parse input CIGAR to get anchors
    struct PairwiseAlignment *pA;
    pA = cigarRead(fileHandleIn);

    // Alignment Parameters //
    // make the pairwise alignment parameters
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
    p->threshold = threshold;
    p->constraintDiagonalTrim = constraintTrim;
    p->diagonalExpansion = diagExpansion;

    // HDP routines //
    // load HDPs
    NanoporeHDP *nHdpT, *nHdpC;
    // check
    if ((templateHdp != NULL) || (complementHdp != NULL)) {
        if ((templateHdp == NULL) || (complementHdp == NULL && twoD)) {
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
    if ((forwardReference == NULL) || (backwardReference == NULL)) {
        st_errAbort("[signalAlign] - ERROR: did not get reference files %s %s\n",
                    forwardReference, backwardReference);
    }
    R = signalUtils_ReferenceSequenceConstructFull(forwardReference, backwardReference, pA);
    
    // Nanopore Read //
    // load nanopore read
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);

    // constrain the event sequence to the positions given by the guide alignment
    Sequence *tEventSequence = makeEventSequenceFromPairwiseAlignment(npRead->templateEvents,
                                                                      pA->start2, pA->end2,
                                                                      (twoD ? npRead->templateEventMap : 
                                                                       npRead->templateStrandEventMap));

    Sequence *cEventSequence;
    if (twoD) {
        cEventSequence = makeEventSequenceFromPairwiseAlignment(npRead->complementEvents,
                                                                pA->start2, pA->end2,
                                                                npRead->complementEventMap);
    } else {
        cEventSequence = NULL;
    }

    // the aligned pairs start at (0,0) so we need to correct them based on the guide alignment later.
    // record the pre-zeroed alignment start and end coordinates here
    // for the events:
    int64_t tCoordinateShift = twoD ? npRead->templateEventMap[pA->start2] : npRead->templateStrandEventMap[pA->start2];
    int64_t cCoordinateShift = twoD ? npRead->complementEventMap[pA->start2] : 0;

    // and for the reference:
    int64_t rCoordinateShift_t = pA->start1;
    int64_t rCoordinateShift_c = twoD ? pA->end1 : 0;
    bool forward = pA->strand1;  // keep track of whether this is a forward mapped read or not

    stList *anchorPairs = signalUtils_guideAlignmentToRebasedAnchorPairs(pA, p);  // pA gets modified here, no turning back

    if ((templateExpectationsFile != NULL) || (complementExpectationsFile != NULL)) {
        st_uglyf("Starting expectations routine\n");
        // Expectation Routine //
        StateMachine *sMt = buildStateMachine(templateModelFile, npRead->templateParams, sMtype, nHdpT);

        // temporary way to 'turn off' estimates if I want to
        if (ESTIMATE_PARAMS) {                                                     //todo remove threshold, not used
            signalUtils_estimateNanoporeParams(sMt, npRead, &npRead->templateParams, ASSIGNMENT_THRESHOLD,
                                               signalUtils_templateOneDAssignmentsFromRead,
                                               nanopore_adjustTemplateEventsForDrift);
        }
        // make empty HMM to collect expectations
        Hmm *templateExpectations = hmmContinuous_getExpectationsHmm(sMt, p->threshold, 0.001, 0.001);

        // get expectations for template
        fprintf(stderr, "signalAlign - getting expectations for template\n");

        getSignalExpectations(sMt, templateExpectations, tEventSequence,
                              (twoD ? npRead->templateEventMap : npRead->templateStrandEventMap),
                              pA->start2,
                              R->getTemplateTargetSequence(R),
                              p, anchorPairs, degenerate);

        if (sMtype == threeStateHdp) {
            fprintf(stderr, "signalAlign - got %" PRId64 "template HDP assignments\n",
                    hmmContinuous_howManyAssignments(templateExpectations));
        }

        // write to file
        fprintf(stderr, "signalAlign - writing expectations to file: %s\n", templateExpectationsFile);

        hmmContinuous_writeToFile(templateExpectationsFile, templateExpectations, sMtype);

        // get expectations for the complement
        StateMachine *sMc;
        Hmm *complementExpectations = NULL;
        if (twoD) {
            fprintf(stderr, "signalAlign - getting expectations for complement\n");

            sMc = buildStateMachine(complementModelFile, npRead->complementParams, sMtype, nHdpC);

            if (ESTIMATE_PARAMS) {
                signalUtils_estimateNanoporeParams(sMc, npRead, &npRead->complementParams, ASSIGNMENT_THRESHOLD,
                                                   signalUtils_complementOneDAssignmentsFromRead,
                                                   nanopore_adjustComplementEventsForDrift);
            }

            complementExpectations = hmmContinuous_getExpectationsHmm(sMc, p->threshold, 0.001, 0.001);
            
            getSignalExpectations(sMc, complementExpectations, cEventSequence, npRead->complementEventMap,
                                  pA->start2,
                                  R->getComplementTargetSequence(R),
                                  p, anchorPairs, degenerate);

            if (sMtype == threeStateHdp) {
                fprintf(stderr, "signalAlign - got %"PRId64"complement HDP assignments\n",
                        hmmContinuous_howManyAssignments(complementExpectations));
            }
            // write to file
            fprintf(stderr, "signalAlign - writing expectations to file: %s\n", complementExpectationsFile);
            hmmContinuous_writeToFile(complementExpectationsFile, complementExpectations, sMtype);
        }


        stateMachine_destruct(sMt);
        signalUtils_ReferenceSequenceDestruct(R);
        hmmContinuous_destruct(templateExpectations, sMtype);
        nanopore_nanoporeReadDestruct(npRead);
        sequence_destruct(tEventSequence);
        pairwiseAlignmentBandingParameters_destruct(p);
        destructPairwiseAlignment(pA);
        stList_destruct(anchorPairs);
        if (twoD) {
            stateMachine_destruct(sMc);
            sequence_destruct(cEventSequence);
            hmmContinuous_destruct(complementExpectations, sMtype);
        }
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
        if (sMtype == threeStateHdp) {
            stateMachine3_setModelToHdpExpectedValues(sMt, nHdpT);
        }

        stList *templateAlignedPairs = performSignalAlignment(sMt, tEventSequence, 
                                                              (twoD ? npRead->templateEventMap : npRead->templateStrandEventMap),
                                                              pA->start2, R->getTemplateTargetSequence(R),
                                                              p, anchorPairs,
                                                              degenerate);
        double templatePosteriorScore = scoreByPosteriorProbabilityIgnoringGaps(templateAlignedPairs);

        // sort
        stList_sort(templateAlignedPairs, sortByXPlusYCoordinate2); //Ensure the coordinates are increasing

        // write to file
        if (posteriorProbsFile != NULL) {
            outputAlignment(outFmt, posteriorProbsFile, readLabel, sMt, npRead->templateParams, npRead->templateEvents,
                            R->getTemplateTargetSequence(R), forward, pA->contig1, tCoordinateShift, rCoordinateShift_t,
                            templateAlignedPairs, templatePosteriorScore,template);
        }

        stList *complementAlignedPairs;
        double complementPosteriorScore = 0.0;
        StateMachine *sMc;
        if (twoD) {
            // Complement alignment
            fprintf(stderr, "signalAlign - starting complement alignment\n");
            sMc = buildStateMachine(complementModelFile, npRead->complementParams, sMtype, nHdpC);

            if (ESTIMATE_PARAMS) {
                signalUtils_estimateNanoporeParams(sMc, npRead, &npRead->complementParams, ASSIGNMENT_THRESHOLD,
                                                   signalUtils_complementOneDAssignmentsFromRead,
                                                   nanopore_adjustComplementEventsForDrift);
            }

            if (sMtype == threeStateHdp) {
                stateMachine3_setModelToHdpExpectedValues(sMc, nHdpC);
            }

            complementAlignedPairs = performSignalAlignment(sMc, cEventSequence,
                                                            npRead->complementEventMap, pA->start2,
                                                            R->getComplementTargetSequence(R),
                                                            p, anchorPairs, degenerate);

            complementPosteriorScore = scoreByPosteriorProbabilityIgnoringGaps(complementAlignedPairs);

            // sort
            stList_sort(complementAlignedPairs, sortByXPlusYCoordinate2); //Ensure the coordinates are increasing

            // write to file
            if (posteriorProbsFile != NULL) {
                outputAlignment(outFmt, posteriorProbsFile, readLabel, sMc, npRead->complementParams,
                                npRead->complementEvents, R->getComplementTargetSequence(R), forward, pA->contig1,
                                cCoordinateShift, rCoordinateShift_c, complementAlignedPairs, complementPosteriorScore,
                                complement);
            }

        }

        fprintf(stdout, "%s %"PRId64"\t%"PRId64"(%f)\t", readLabel, stList_length(anchorPairs),
                stList_length(templateAlignedPairs), templatePosteriorScore);
        if (twoD) {
            fprintf(stdout, "%"PRId64"(%f)\n", stList_length(complementAlignedPairs), complementPosteriorScore);
        } else {
            fprintf(stdout, "\n");
        }

        // final alignment clean up
        destructPairwiseAlignment(pA);
        nanopore_nanoporeReadDestruct(npRead);
        signalUtils_ReferenceSequenceDestruct(R);
        stateMachine_destruct(sMt);
        sequence_destruct(tEventSequence);
        stList_destruct(templateAlignedPairs);
        if (twoD) {
            stateMachine_destruct(sMc);
            sequence_destruct(cEventSequence);
            stList_destruct(complementAlignedPairs);
        }
        fprintf(stderr, "signalAlign - SUCCESS: finished alignment of query %s, exiting\n", readLabel);
    }
    return 0;
}

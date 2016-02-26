#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <inttypes.h>
#include "nanopore.h"
#include "pairwiseAligner.h"


static NanoporeRead *nanopore_nanoporeReadConstruct(int64_t readLength,
                                                    int64_t nbTemplateEvents,
                                                    int64_t nbComplementEvents) {
    NanoporeRead *npRead = st_malloc(sizeof(NanoporeRead));

    npRead->readLength = readLength;
    npRead->nbTemplateEvents = nbTemplateEvents;
    npRead->nbComplementEvents = nbComplementEvents;

    npRead->twoDread = st_malloc(npRead->readLength * sizeof(char));

    // the map contains the index of the event corresponding to each kmer in the read sequence so
    // the length of the map has to be the same as the read sequence, not the number of events
    // there can be events that aren't mapped to a kmer
    npRead->templateEventMap = st_malloc(npRead->readLength * sizeof(int64_t));
    npRead->templateEvents = st_malloc(npRead->nbTemplateEvents * NB_EVENT_PARAMS * sizeof(double));

    npRead->complementEventMap = st_malloc(npRead->readLength * sizeof(int64_t));
    npRead->complementEvents = st_malloc(npRead->nbComplementEvents * NB_EVENT_PARAMS * sizeof(double));

    npRead->scaled = TRUE;
    // return
    return npRead;
}

static void nanopore_descaleEvents(int64_t nb_events, double *events, double scale, double shift) {
    for (int64_t i = 0; i < nb_events; i += NB_EVENT_PARAMS) {
        events[i] = (events[i] - shift) / scale;
    }
}

NanoporeRead *nanopore_loadNanoporeReadFromFile(const char *nanoporeReadFile) {
    FILE *fH = fopen(nanoporeReadFile, "r");

    // line 1 [2D read length] [# of template events] [# of complement events]
    //        [template scale] [template shift] [template var] [template scale_sd] [template var_sd]
    //        [complement scale] [complement shift] [complement var] [complement scale_sd] [complement var_sd] \n
    char *string = stFile_getLineFromFile(fH);
    stList *tokens = stString_split(string);
    int64_t npRead_readLength, npRead_nbTemplateEvents, npRead_nbComplementEvents;
    assert(stList_length(tokens) == 12);

    int64_t j = sscanf(stList_get(tokens, 0), "%lld", &npRead_readLength);
    if (j != 1) {
        st_errAbort("error parsing nanopore read length\n");
    }
    j = sscanf(stList_get(tokens, 1), "%lld", &npRead_nbTemplateEvents);
    if (j != 1) {
        st_errAbort("error parsing nanopore template event lengths\n");
    }
    j = sscanf(stList_get(tokens, 2), "%lld", &npRead_nbComplementEvents);
    if (j != 1) {
        st_errAbort("error parsing nanopore complement event lengths\n");
    }
    NanoporeRead *npRead = nanopore_nanoporeReadConstruct(npRead_readLength,
                                                          npRead_nbTemplateEvents,
                                                          npRead_nbComplementEvents);
    // template params
    j = sscanf(stList_get(tokens, 3), "%lf", &(npRead->templateParams.scale));
    if (j != 1) {
        st_errAbort("error parsing nanopore template scale\n");
    }
    j = sscanf(stList_get(tokens, 4), "%lf", &(npRead->templateParams.shift));
    if (j != 1) {
        st_errAbort("error parsing nanopore template shift\n");
    }
    j = sscanf(stList_get(tokens, 5), "%lf", &(npRead->templateParams.var));
    if (j != 1) {
        st_errAbort("error parsing nanopore template var\n");
    }
    j = sscanf(stList_get(tokens, 6), "%lf", &(npRead->templateParams.scale_sd));
    if (j != 1) {
        st_errAbort("error parsing nanopore template scale_sd\n");
    }
    j = sscanf(stList_get(tokens, 7), "%lf", &(npRead->templateParams.var_sd));
    if (j != 1) {
        st_errAbort("error parsing nanopore template var_sd\n");
    }

    // complement params
    j = sscanf(stList_get(tokens, 8), "%lf", &(npRead->complementParams.scale));
    if (j != 1) {
        st_errAbort("error parsing nanopore complement scale\n");
    }
    j = sscanf(stList_get(tokens, 9), "%lf", &(npRead->complementParams.shift));
    if (j != 1) {
        st_errAbort("error parsing nanopore complement shift\n");
    }
    j = sscanf(stList_get(tokens, 10), "%lf", &(npRead->complementParams.var));
    if (j != 1) {
        st_errAbort("error parsing nanopore complement var\n");
    }
    j = sscanf(stList_get(tokens, 11), "%lf", &(npRead->complementParams.scale_sd));
    if (j != 1) {
        st_errAbort("error parsing nanopore complement scale_sd\n");
    }
    j = sscanf(stList_get(tokens, 12), "%lf", &(npRead->complementParams.var_sd));
    if (j != 1) {
        st_errAbort("error parsing nanopore template var_sd\n");
    }

    // cleanup line 1
    free(string);
    stList_destruct(tokens);

    // line 2 [2D read sequence] \n
    string = stFile_getLineFromFile(fH);
    j = sscanf(string, "%s", npRead->twoDread);
    if (j != 1) {
        st_errAbort("error parsing read from npRead file\n");
    }
    free(string);

    // line 3 [template event map] \n
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    // check for correctness
    if (stList_length(tokens) != npRead->readLength) {
        st_errAbort(
                "template event map is not the correct length, should be %lld, got %lld",
                npRead->readLength,
                stList_length(tokens));
    }
    for (int64_t i = 0; i < npRead->readLength; i++) {
        j = sscanf(stList_get(tokens, i), "%lld", &(npRead->templateEventMap[i]));
        if (j != 1) {
            st_errAbort("error loading in template eventMap\n");
        }
    }
    free(string);
    stList_destruct(tokens);

    // line 4 [template events (mean, stddev, length)] \n
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    // check
    if (stList_length(tokens) != (npRead->nbTemplateEvents * NB_EVENT_PARAMS)) {
        st_errAbort(
                "incorrect number of template events, should be %lld, got %lld",
                npRead->nbTemplateEvents,
                stList_length(tokens));
    }
    for (int64_t i = 0; i < (npRead->nbTemplateEvents * NB_EVENT_PARAMS); i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(npRead->templateEvents[i]));
        if (j != 1) {
            st_errAbort("error loading in template events\n");
        }
    }
    free(string);
    stList_destruct(tokens);

    // line 5 [complement event map] \n
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    // check for correctness
    if (stList_length(tokens) != npRead->readLength) {
        st_errAbort(
                "complament event map is not the correct length, should be %lld, got %lld",
                npRead->readLength,
                stList_length(tokens));
    }
    for (int64_t i = 0; i < npRead->readLength; i++) {
        j = sscanf(stList_get(tokens, i), "%lld", &(npRead->complementEventMap[i]));
        if (j != 1) {
            st_errAbort("error loading in complement eventMap\n");
        }
    }
    free(string);
    stList_destruct(tokens);

    // line 6 [complement events (mean, stddev, length)] \n
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    // check
    if (stList_length(tokens) != (npRead->nbComplementEvents * NB_EVENT_PARAMS)) {
        st_errAbort(
                "incorrect number of complement events, should be %lld, got %lld",
                npRead->nbComplementEvents,
                stList_length(tokens));
    }
    for (int64_t i = 0; i < (npRead->nbComplementEvents * NB_EVENT_PARAMS); i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(npRead->complementEvents[i]));
        if (j != 1) {
            st_errAbort("error loading in complement events\n");
        }
    }
    free(string);
    stList_destruct(tokens);

    fclose(fH);
    return npRead;
}

stList *nanopore_remapAnchorPairs(stList *anchorPairs, int64_t *eventMap) {
    stList *mappedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    for (int64_t i = 0; i < stList_length(anchorPairs); i++) {
        stIntTuple *pair = stList_get(anchorPairs, i);
        stList_append(mappedPairs,
                      stIntTuple_construct2(stIntTuple_get(pair, 0), eventMap[stIntTuple_get(pair, 1)]));
    }

    return mappedPairs;
}

stList *nanopore_remapAnchorPairsWithOffset(stList *unmappedPairs, int64_t *eventMap, int64_t mapOffset) {
    stList *mappedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    for (int64_t i = 0; i < stList_length(unmappedPairs); i++) {
        stIntTuple *pair = stList_get(unmappedPairs, i);

        stList_append(mappedPairs,
                      stIntTuple_construct2(stIntTuple_get(pair, 0), eventMap[stIntTuple_get(pair, 1)] -
                              eventMap[mapOffset]));
    }

    return mappedPairs;
}

void nanopore_descaleNanoporeRead(NanoporeRead *npRead) {
    nanopore_descaleEvents(npRead->nbTemplateEvents, npRead->templateEvents,
                           npRead->templateParams.scale,
                           npRead->templateParams.shift);
    nanopore_descaleEvents(npRead->nbComplementEvents, npRead->complementEvents,
                           npRead->complementParams.scale,
                           npRead->complementParams.shift);
    npRead->scaled = FALSE;
}

void nanopore_nanoporeReadDestruct(NanoporeRead *npRead) {
    free(npRead->twoDread);
    free(npRead->templateEventMap);
    free(npRead->templateEvents);
    free(npRead->complementEventMap);
    free(npRead->complementEvents);
    free(npRead);
}

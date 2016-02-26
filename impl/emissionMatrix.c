#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include "emissionMatrix.h"

/*
 * Code-generated emissions matrix for 2-mers
 */

void emissions_kmer_setMatchProbsToDefaults(double *emissionMatchProbs) {
    /*
     * This sets the match probabilities to default values for matching kmers, not really used anymore...
     */

    const double M=-2.1149196655034745; //log(0.12064298095701059);
    const double V=-4.5691014376830479; //log(0.010367271172731285);
    const double S=-3.9833860032220842; //log(0.01862247669752685);
    const double N=-2.772588722;        //log(0.25**2)

    //Symmetric matrix of emission probabilities.

    const double i[MATRIX_SIZE] = {
        M+M, M+V, M+S, M+V, M+N, V+M, V+V, V+S, V+V, V+N, S+M, S+V, S+S, S+V, S+N, V+M, V+V, V+S, V+V, V+N, N+M, N+V, N+S, N+V, N+N,
        M+V, M+M, M+V, M+S, M+N, V+V, V+M, V+V, V+S, V+N, S+V, S+M, S+V, S+S, S+N, V+V, V+M, V+V, V+S, V+N, N+V, N+M, N+V, N+S, N+N,
        M+S, M+V, M+M, M+V, M+N, V+S, V+V, V+M, V+V, V+N, S+S, S+V, S+M, S+V, S+N, V+S, V+V, V+M, V+V, V+N, N+S, N+V, N+M, N+V, N+N,
        M+V, M+S, M+V, M+M, M+N, V+V, V+S, V+V, V+M, V+N, S+V, S+S, S+V, S+M, S+N, V+V, V+S, V+V, V+M, V+N, N+V, N+S, N+V, N+M, N+N,
        M+N, M+N, M+N, M+N, M+N, V+N, V+N, V+N, V+N, V+N, S+N, S+N, S+N, S+N, S+N, V+N, V+N, V+N, V+N, V+N, N+N, N+N, N+N, N+N, N+N,
        V+M, V+V, V+S, V+V, V+N, M+M, M+V, M+S, M+V, M+N, V+M, V+V, V+S, V+V, V+N, S+M, S+V, S+S, S+V, S+N, N+M, N+V, N+S, N+V, N+N,
        V+V, V+M, V+V, V+S, V+N, M+V, M+M, M+V, M+S, M+N, V+V, V+M, V+V, V+S, V+N, S+V, S+M, S+V, S+S, S+N, N+V, N+M, N+V, N+S, N+N,
        V+S, V+V, V+M, V+V, V+N, M+S, M+V, M+M, M+V, M+N, V+S, V+V, V+M, V+V, V+N, S+S, S+V, S+M, S+V, S+N, N+S, N+V, N+M, N+V, N+N,
        V+V, V+S, V+V, V+M, V+N, M+V, M+S, M+V, M+M, M+N, V+V, V+S, V+V, V+M, V+N, S+V, S+S, S+V, S+M, S+N, N+V, N+S, N+V, N+M, N+N,
        V+N, V+N, V+N, V+N, V+N, M+N, M+N, M+N, M+N, M+N, V+N, V+N, V+N, V+N, V+N, S+N, S+N, S+N, S+N, S+N, N+N, N+N, N+N, N+N, N+N,
        S+M, S+V, S+S, S+V, S+N, V+M, V+V, V+S, V+V, V+N, M+M, M+V, M+S, M+V, M+N, V+M, V+V, V+S, V+V, V+N, N+M, N+V, N+S, N+V, N+N,
        S+V, S+M, S+V, S+S, S+N, V+V, V+M, V+V, V+S, V+N, M+V, M+M, M+V, M+S, M+N, V+V, V+M, V+V, V+S, V+N, N+V, N+M, N+V, N+S, N+N,
        S+S, S+V, S+M, S+V, S+N, V+S, V+V, V+M, V+V, V+N, M+S, M+V, M+M, M+V, M+N, V+S, V+V, V+M, V+V, V+N, N+S, N+V, N+M, N+V, N+N,
        S+V, S+S, S+V, S+M, S+N, V+V, V+S, V+V, V+M, V+N, M+V, M+S, M+V, M+M, M+N, V+V, V+S, V+V, V+M, V+N, N+V, N+S, N+V, N+M, N+N,
        S+N, S+N, S+N, S+N, S+N, V+N, V+N, V+N, V+N, V+N, M+N, M+N, M+N, M+N, M+N, V+N, V+N, V+N, V+N, V+N, N+N, N+N, N+N, N+N, N+N,
        V+M, V+V, V+S, V+V, V+N, S+M, S+V, S+S, S+V, S+N, V+M, V+V, V+S, V+V, V+N, M+M, M+V, M+S, M+V, M+N, N+M, N+V, N+S, N+V, N+N,
        V+V, V+M, V+V, V+S, V+N, S+V, S+M, S+V, S+S, S+N, V+V, V+M, V+V, V+S, V+N, M+V, M+M, M+V, M+S, M+N, N+V, N+M, N+V, N+S, N+N,
        V+S, V+V, V+M, V+V, V+N, S+S, S+V, S+M, S+V, S+N, V+S, V+V, V+M, V+V, V+N, M+S, M+V, M+M, M+V, M+N, N+S, N+V, N+M, N+V, N+N,
        V+V, V+S, V+V, V+M, V+N, S+V, S+S, S+V, S+M, S+N, V+V, V+S, V+V, V+M, V+N, M+V, M+S, M+V, M+M, M+N, N+V, N+S, N+V, N+M, N+N,
        V+N, V+N, V+N, V+N, V+N, S+N, S+N, S+N, S+N, S+N, V+N, V+N, V+N, V+N, V+N, M+N, M+N, M+N, M+N, M+N, N+N, N+N, N+N, N+N, N+N,
        N+M, N+V, N+S, N+V, N+N, N+M, N+V, N+S, N+V, N+N, N+M, N+V, N+S, N+V, N+N, N+M, N+V, N+S, N+V, N+N, N+M, N+V, N+S, N+V, N+N,
        N+V, N+M, N+V, N+S, N+N, N+V, N+M, N+V, N+S, N+N, N+V, N+M, N+V, N+S, N+N, N+V, N+M, N+V, N+S, N+N, N+V, N+M, N+V, N+S, N+N,
        N+S, N+V, N+M, N+V, N+N, N+S, N+V, N+M, N+V, N+N, N+S, N+V, N+M, N+V, N+N, N+S, N+V, N+M, N+V, N+N, N+S, N+V, N+M, N+V, N+N,
        N+V, N+S, N+V, N+M, N+N, N+V, N+S, N+V, N+M, N+N, N+V, N+S, N+V, N+M, N+N, N+V, N+S, N+V, N+M, N+N, N+V, N+S, N+V, N+M, N+N,
        N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N, N+N,
        };

    memcpy(emissionMatchProbs, i, sizeof(double)*MATRIX_SIZE);
}


void emissions_kmer_setGapProbsToDefaults(double *emissionGapProbs) {
    /*
     * This is used to set the emissions to reasonable values.
     */
    //const double G = -1.6094379124341003; //log(0.2)
    const double G = -3.2188758248682006; //log(0.2)+log(0.2)
    const double i[25] = {
        G, G, G, G, G,
        G, G, G, G, G,
        G, G, G, G, G,
        G, G, G, G, G,
        G, G, G, G, G,
        };

    memcpy(emissionGapProbs, i, sizeof(double)*25);
}

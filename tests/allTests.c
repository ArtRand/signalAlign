/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <CuTest.h>
#include "sonLib.h"
#include "CuTest.h"
#include "sonLibCommon.h"

CuSuite *NanoporeHdpTestSuite(void);
CuSuite *HdpTestSuite(void);
CuSuite *variableOrderPairwiseAlignerTestSuite(void);
CuSuite *signalPairwiseAlignerTestSuite(void);
CuSuite *stateMachineAlignmentTestSuite(void);
// multiple alignment not implemented, yet
//CuSuite* multipleAlignerTestSuite(void);
//CuSuite* pairwiseAlignmentLongTestSuite(void);  // legacy to remind myself


CuSuite *stBaseAlignerRunAllTests(void) {
    CuString *output = CuStringNew();
    CuSuite *suite = CuSuiteNew();
    CuSuiteAddSuite(suite, NanoporeHdpTestSuite());
    CuSuiteAddSuite(suite, HdpTestSuite());
    CuSuiteAddSuite(suite, signalPairwiseAlignerTestSuite());
    CuSuiteAddSuite(suite, variableOrderPairwiseAlignerTestSuite());
    CuSuiteAddSuite(suite, stateMachineAlignmentTestSuite());

    // coming soon..
    //CuSuiteAddSuite(suite, multipleAlignerTestSuite());
    //CuSuiteAddSuite(suite, pairwiseAlignmentLongTestSuite());

    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s\n", output->buffer);
    CuStringDelete(output);

    return suite;
}

int main(int argc, char *argv[]) {
    if(argc == 2) {
        st_setLogLevelFromString(argv[1]);
    }
    CuSuite *allTests = stBaseAlignerRunAllTests();
    int good = allTests->failCount > 0;
    CuSuiteDelete(allTests);
    return good;
}

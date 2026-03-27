/* linbox/tests/test-frobenius-suite.C
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 * Written by Omesh Dhar Dwivedi <odd23@drexel.edu>
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 * ========LICENCE========
 */

#include "linbox/linbox-config.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>

#include "linbox/ring/modular.h"
#include "linbox/util/commentator.h"
#include "linbox/ring/ntl.h"

#include "linbox/algorithms/frobenius-large.h"
#include "linbox/algorithms/frobenius-large-bf.h"
#include "linbox/algorithms/frobenius-large-dense.h"

#include "test-frobenius-suite.h"

using namespace LinBox;

// Parse algorithm mask from string like "0", "01", "012", "2" etc.
// bit0 = Toeplitz, bit1 = Butterfly, bit2 = Dense
// Default (empty string or "012") = all three
int parseAlgoMask(const char *s) {
    if (s == nullptr || strlen(s) == 0) return 7; // all
    int mask = 0;
    for (size_t i = 0; i < strlen(s); ++i) {
        if (s[i] == '0') mask |= 1;
        if (s[i] == '1') mask |= 2;
        if (s[i] == '2') mask |= 4;
    }
    return mask ? mask : 7;
}

int main(int argc, char **argv) {
    // Defaults
    uint64_t p    = 10000019;
    uint64_t e    = 1;
    size_t   k    = 0;
    int      seed = time(NULL);
    size_t   s    = 0;   // 0 = use default test cases
    size_t   w    = 0;
    size_t   h    = 0;
    int      poly = 0;   // 0=x, 1=x-1, 2=x+1
    char     algoStr[16] = "";  // empty = all

    static Argument args[] = {
        { 'k', "-k K", "Number of invariant factors to compute (0 = all)",
            TYPE_INT, &k },
        { 'p', "-p P", "Characteristic of field GF(p^e)",
            TYPE_INT, &p },
        { 'e', "-e E", "Extension degree of field GF(p^e)",
            TYPE_INT, &e },
        { 'r', "-r R", "Random seed",
            TYPE_INT, &seed },
        { 's', "-s S", "Custom: number of distinct block sizes (0 = use defaults)",
            TYPE_INT, &s },
        { 'w', "-w W", "Custom: width (repetitions per block size)",
            TYPE_INT, &w },
        { 'H', "-H H", "Custom: height (step between block sizes)",
            TYPE_INT, &h },
        { 'q', "-q Q", "Custom: polynomial (0=x, 1=x-1, 2=x+1)",
            TYPE_INT, &poly },
        { 'a', "-a A", "Algorithms to run: string of digits 0=Toeplitz 1=Butterfly 2=Dense (default: all)",
            TYPE_STR, algoStr },
        END_OF_ARGUMENTS
    };

    parseArguments(argc, argv, args);
    srand(seed);

    typedef NTL_zz_p   Field;
    typedef NTL_zz_pX  PolyRing;

    Field   F(p, e);
    PolyRing R(F);

    FrobeniusLarge<PolyRing>         FT(R);
    FrobeniusLargeButterfly<PolyRing> FB(R);
    FrobeniusLargeDense<PolyRing>    FD(R);

    int algoMask = parseAlgoMask(algoStr[0] ? algoStr : nullptr);

    std::cout << "=== Frobenius Test Suite ===" << std::endl;
    std::cout << "Field: GF(" << p << "^" << e << ")"
              << "  seed=" << seed
              << "  k=" << k << std::endl;
    std::cout << "Algorithms: "
              << ((algoMask & 1) ? "Toeplitz " : "")
              << ((algoMask & 2) ? "Butterfly " : "")
              << ((algoMask & 4) ? "Dense " : "")
              << std::endl;

    bool pass = runSuite(
        FT, FB, FD,
        F, R,
        k,
        algoMask,
        seed,
        s, w, h, poly);

        
    return pass ? 0 : -1;
}
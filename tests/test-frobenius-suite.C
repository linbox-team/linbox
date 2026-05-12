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
#include "linbox/algorithms/frobenius-large-search.h"
#include "linbox/algorithms/invariant-factors.h"

#include "test-frobenius-suite.h"

using namespace LinBox;

// bit0=Toeplitz, bit1=Butterfly, bit2=Dense, bit3=Search, bit4=LIFs
int parseAlgoMask(const char *s) {
    if (s == nullptr || strlen(s) == 0) return 31; // all
    int mask = 0;
    for (size_t i = 0; i < strlen(s); ++i) {
        if (s[i] == '0') mask |= 1;
        if (s[i] == '1') mask |= 2;
        if (s[i] == '2') mask |= 4;
        if (s[i] == '3') mask |= 8;
        if (s[i] == '4') mask |= 16;
    }
    return mask ? mask : 31;
}

int main(int argc, char **argv) {
    int    p    = 10000019;
    int    e    = 1;
    size_t k    = 0;
    int    seed = time(NULL);
    size_t s    = 0;
    size_t w    = 0;
    size_t h    = 0;
    int    poly = 0;
    // TYPE_STR expects std::string* not char[] — using char[] caused
    // parseArguments to call std::string::assign() on a zero-filled char
    // array, crashing in memmove (EXC_BAD_ACCESS address=0x0).
    std::string algoStr = "";

    static Argument args[] = {
        { 'k', "-k K", "Number of invariant factors to compute (0 = all)", TYPE_INT, &k },
        { 'p', "-p P", "Characteristic of field GF(p^e)",                  TYPE_INT, &p },
        { 'e', "-e E", "Extension degree of field GF(p^e)",                 TYPE_INT, &e },
        { 'r', "-r R", "Random seed",                                       TYPE_INT, &seed },
        { 's', "-s S", "Custom: number of distinct block sizes",            TYPE_INT, &s },
        { 'w', "-w W", "Custom: width (repetitions per block size)",        TYPE_INT, &w },
        { 'H', "-H H", "Custom: height (step between block sizes)",         TYPE_INT, &h },
        { 'q', "-q Q", "Custom: polynomial (0=x, 1=x-1, 2=x+1)",          TYPE_INT, &poly },
        { 'a', "-a A", "Algorithms: 0=Toeplitz 1=Butterfly 2=Dense 3=Search 4=LIFs", TYPE_STR, &algoStr },
        END_OF_ARGUMENTS
    };

    parseArguments(argc, argv, args);
    srand(seed);

    typedef NTL_zz_p  Field;
    typedef NTL_zz_pX PolyRing;

    Field    F(p, e);
    PolyRing R(F);

    FrobeniusLarge<PolyRing>          FT(R);
    FrobeniusLargeButterfly<PolyRing> FB(R);
    FrobeniusLargeDense<PolyRing>     FD(R);
    FrobeniusLargeSearch<PolyRing>    FS(R);
    InvariantFactors<Field, PolyRing> IFD(F, R);

    int algoMask = parseAlgoMask(algoStr.empty() ? nullptr : algoStr.c_str());

    std::cout << "=== Frobenius Test Suite ===" << std::endl;
    std::cout << "Field: GF(" << p << "^" << e << ")"
              << "  seed=" << seed << "  k=" << k << std::endl;
    std::cout << "Algorithms: "
              << ((algoMask &  1) ? "Toeplitz "  : "")
              << ((algoMask &  2) ? "Butterfly " : "")
              << ((algoMask &  4) ? "Dense "     : "")
              << ((algoMask &  8) ? "Search "    : "")
              << ((algoMask & 16) ? "LIFs "      : "")
              << std::endl;

    bool pass = runSuite(FT, FB, FD, FS, IFD, F, R, k, algoMask, seed, s, w, h, poly);
    return pass ? 0 : -1;
}
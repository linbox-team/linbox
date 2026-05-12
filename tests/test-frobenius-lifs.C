/* linbox/tests/test-frobenius-lifs.C
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
 * 
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

#include "linbox/algorithms/invariant-factors.h"

#include "test-frobenius-suite.h"

using namespace LinBox;

int main(int argc, char **argv) {
    int    p    = 3;
    size_t k    = 0;
    int    seed = time(NULL);
    size_t s    = 0;
    size_t w    = 0;
    size_t h    = 0;
    int    poly = 0;

    static Argument args[] = {
        { 'p', "-p P", "Characteristic of field GF(p)",
            TYPE_INT, &p },
        { 'k', "-k K", "Number of invariant factors to compute (0 = all)",
            TYPE_INT, &k },
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
        END_OF_ARGUMENTS
    };

    parseArguments(argc, argv, args);
    srand(seed);

    typedef Givaro::Modular<double> Field;
    typedef NTL_zz_pX PolyRing;
    typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;
    typedef typename PolyRing::Element Polynomial;

    Field    F(p);
    PolyRing R(p);

    InvariantFactors<Field, PolyRing> IFD(F, R);

    std::cout << "=== Frobenius LIFs Test Suite ===" << std::endl;
    std::cout << "Field: GF(" << p << ")"
              << "  seed=" << seed
              << "  k=" << k << std::endl;

    std::vector<TestParams> cases = defaultTestCases();
    if (s > 0 && w > 0 && h > 0) {
        std::string name = "custom_s" + std::to_string(s)
                         + "_w" + std::to_string(w)
                         + "_h" + std::to_string(h);
        cases.insert(cases.begin(), {s, w, h, poly, name});
    }

    std::vector<TestResult> results;
    bool allPass = true;

    for (auto &tc : cases) {
        size_t n = tc.w * tc.h * tc.s * (tc.s + 1) / 2;
        SparseMat M(F, n, n);
        std::vector<Polynomial> expected;
        buildJordanMatrix(M, expected, F, R, tc.s, tc.w, tc.h, tc.poly);
        sprayMatrix(M, F, seed);

        auto res = runOne(IFD, "LIFs", F, R, M, expected, tc.name, k);
        allPass &= res.pass;
        results.push_back(res);
    }

    printTable(results);
    return allPass ? 0 : -1;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
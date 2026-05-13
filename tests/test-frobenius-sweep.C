/* linbox/tests/test-frobenius-sweep.C
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 * Written by Omesh Dhar Dwivedi <odd23@drexel.edu>
 *
 * Sweeps matrix sizes by varying s from s_min to s_max (with fixed w, h, poly)
 * to find where algorithm crossovers occur. Focuses on ButterflySearch vs
 * DenseSearch by default but any algorithm subset can be selected.
 *
 * Matrix dimension: n = w * h * s * (s+1) / 2
 *
 * Usage:
 *   ./test-frobenius-sweep [-S smin] [-E smax] [-w W] [-H H] [-q Q]
 *                          [-k K] [-p P] [-r R] [-a ALGOS]
 *
 * -a accepts digits:
 *   0=Toeplitz  1=ToeplitzSearch  2=Butterfly  3=ButterflySearch
 *   4=Dense     5=DenseSearch     6=LIFs
 *   default: ButterflySearch + DenseSearch (35)
 *
 * ========LICENCE========
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License.
 * ========LICENCE========
 */

#include "linbox/linbox-config.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>

#include "linbox/ring/modular.h"
#include "linbox/util/commentator.h"
#include "linbox/ring/ntl.h"

#include "linbox/algorithms/frobenius-large-generic.h"

#include "linbox/algorithms/invariant-factors.h"

#include "test-frobenius-suite.h"

using namespace LinBox;

// bit0=Toeplitz, bit1=ToeplitzSearch, bit2=Butterfly, bit3=ButterflySearch,
// bit4=Dense,    bit5=DenseSearch,    bit6=LIFs
int parseAlgoMask(const char *s) {
    if (s == nullptr || strlen(s) == 0) return 40; // default: ButterflySearch + DenseSearch
    int mask = 0;
    for (size_t i = 0; i < strlen(s); ++i) {
        if (s[i] == '0') mask |= 1;
        if (s[i] == '1') mask |= 2;
        if (s[i] == '2') mask |= 4;
        if (s[i] == '3') mask |= 8;
        if (s[i] == '4') mask |= 16;
        if (s[i] == '5') mask |= 32;
        if (s[i] == '6') mask |= 64;
    }
    return mask ? mask : 40;
}

struct SweepResult {
    size_t s;
    size_t n;
    size_t nnz;
    double t_toeplitz;
    double t_toeplitz_search;
    double t_butterfly;
    double t_butterfly_search;
    double t_dense;
    double t_dense_search;
    double t_lifs;
};

template<class FrobeniusObject, class SparseMat, class PolyRing>
double timeOne(FrobeniusObject &FO, SparseMat &M, size_t k, const PolyRing &R) {
    typedef typename PolyRing::Element Polynomial;
    std::vector<Polynomial> fs;
    Givaro::Timer T; T.clear(); T.start();
    FO.frobeniusInvariants(fs, M, k);
    T.stop();
    return T.usertime();
}

int main(int argc, char **argv) {
    int    p      = 10000019;
    int    e      = 1;
    size_t k      = 0;       // 0 = use log2(n)
    int    seed   = time(NULL);
    size_t s_min  = 4;
    size_t s_max  = 30;
    size_t w      = 4;
    size_t h      = 4;
    int    poly   = 0;
    char   algoStr[16] = "";

    static Argument args[] = {
        { 'S', "-S S", "Minimum s (controls min matrix size)",               TYPE_INT, &s_min },
        { 'E', "-E E", "Maximum s (controls max matrix size)",               TYPE_INT, &s_max },
        { 'w', "-w W", "Width: repetitions per block size",                  TYPE_INT, &w },
        { 'H', "-H H", "Height: step between block sizes",                   TYPE_INT, &h },
        { 'q', "-q Q", "Polynomial (0=x, 1=x-1, 2=x+1)",                   TYPE_INT, &poly },
        { 'k', "-k K", "Number of invariant factors (0 = log2(n))",         TYPE_INT, &k },
        { 'p', "-p P", "Field characteristic",                               TYPE_INT, &p },
        { 'e', "-e E", "Extension degree",                                   TYPE_INT, &e },
        { 'r', "-r R", "Random seed",                                        TYPE_INT, &seed },
        { 'a', "-a A", "Algorithms: 0=Toeplitz 1=ToeplitzSearch 2=Butterfly 3=ButterflySearch 4=Dense 5=DenseSearch 6=LIFs (default: 35)",
            TYPE_STR, algoStr },
        END_OF_ARGUMENTS
    };

    parseArguments(argc, argv, args);
    srand(seed);

    typedef NTL_zz_p  Field;
    typedef NTL_zz_pX PolyRing;
    typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;
    typedef PolyRing::Element Polynomial;

    Field    F(p, e);
    PolyRing R(F);

    int algoMask = parseAlgoMask(algoStr[0] ? algoStr : nullptr);

    FrobeniusLarge<PolyRing>                FT(R);
    FrobeniusLargeSearch<PolyRing>          FTS(R);
    FrobeniusLargeButterfly<PolyRing>       FB(R);
    FrobeniusLargeButterflySearch<PolyRing> FBS(R);
    FrobeniusLargeDense<PolyRing>           FD(R);
    FrobeniusLargeDenseSearch<PolyRing>     FDS(R);
    InvariantFactors<Field, PolyRing>       IFD(F, R);

    std::cout << "=== Frobenius Algorithm Sweep ===" << std::endl;
    std::cout << "Field: GF(" << p << "^" << e << ")"
              << "  seed=" << seed
              << "  w=" << w << "  h=" << h
              << "  poly=" << poly << std::endl;
    std::cout << "Algorithms: "
              << ((algoMask &  1) ? "Toeplitz "        : "")
              << ((algoMask &  2) ? "ToeplitzSearch "  : "")
              << ((algoMask &  4) ? "Butterfly "       : "")
              << ((algoMask &  8) ? "ButterflySearch " : "")
              << ((algoMask & 16) ? "Dense "           : "")
              << ((algoMask & 32) ? "DenseSearch "     : "")
              << ((algoMask & 64) ? "LIFs "            : "")
              << std::endl;
    std::cout << std::string(100, '-') << std::endl;

    // Print header
    std::cout << std::left
              << std::setw(6)  << "s"
              << std::setw(8)  << "n"
              << std::setw(6)  << "k";
    if (algoMask &  1) std::cout << std::setw(12) << "Toeplitz";
    if (algoMask &  2) std::cout << std::setw(14) << "ToepSearch";
    if (algoMask &  4) std::cout << std::setw(12) << "Butterfly";
    if (algoMask &  8) std::cout << std::setw(14) << "BflySearch";
    if (algoMask & 16) std::cout << std::setw(12) << "Dense";
    if (algoMask & 32) std::cout << std::setw(14) << "DenseSearch";
    if (algoMask & 64) std::cout << std::setw(12) << "LIFs";
    std::cout << std::setw(14) << "Winner" << "\n";
    std::cout << std::string(100, '-') << std::endl;

    for (size_t s = s_min; s <= s_max; ++s) {
        size_t n = w * h * s * (s + 1) / 2;
        size_t keff = (k == 0) ? (size_t)std::max(1.0, std::log2((double)n)) : k;

        SparseMat M(F, n, n);
        std::vector<Polynomial> expected;
        buildJordanMatrix(M, expected, F, R, s, w, h, poly);
        sprayMatrix(M, F, seed);

        // time each selected algorithm
        double times[7] = {-1,-1,-1,-1,-1,-1,-1};
        if (algoMask &  1) times[0] = timeOne(FT,  M, keff, R);
        if (algoMask &  2) times[1] = timeOne(FTS, M, keff, R);
        if (algoMask &  4) times[2] = timeOne(FB,  M, keff, R);
        if (algoMask &  8) times[3] = timeOne(FBS, M, keff, R);
        if (algoMask & 16) times[4] = timeOne(FD,  M, keff, R);
        if (algoMask & 32) times[5] = timeOne(FDS, M, keff, R);
        if (algoMask & 64) times[6] = timeOne(IFD, M, keff, R);

        // find winner among selected
        const char* labels[7] = {"Toeplitz","ToepSearch","Butterfly","BflySearch","Dense","DenseSearch","LIFs"};
        double best = 1e18;
        int winner = -1;
        for (int i = 0; i < 7; ++i) {
            if (times[i] >= 0 && times[i] < best) { best = times[i]; winner = i; }
        }

        std::cout << std::left
                  << std::setw(6)  << s
                  << std::setw(8)  << n
                  << std::setw(6)  << keff;
        for (int i = 0; i < 7; ++i) {
            if (times[i] < 0) continue;
            int w_col = (i == 1 || i == 3 || i == 5) ? 14 : 12;
            std::cout << std::setw(w_col) << std::fixed << std::setprecision(4) << times[i];
        }
        std::cout << std::setw(14) << (winner >= 0 ? labels[winner] : "?") << "\n";
        std::cout.flush();
    }

    return 0;
}
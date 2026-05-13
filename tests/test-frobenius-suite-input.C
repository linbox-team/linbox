/* linbox/tests/test-frobenius-suite-input.C
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 * Written by Omesh Dhar Dwivedi <odd23@drexel.edu>
 *
 * Speed benchmark: read a matrix from a given path in SMS format
 * and compare timing of all preconditioner variants + LIFs.
 * No correctness check — pure timing.
 *
 * Usage:
 *   ./test-frobenius-suite-input -f /path/to/matrix.mat [-p 10007] [-k 0] [-r 42] [-a 0123456] [-n 5]
 *
 * -a accepts digits:
 *   0=Toeplitz  1=ToeplitzSearch  2=Butterfly  3=ButterflySearch
 *   4=Dense     5=DenseSearch     6=LIFs
 *   default: all seven
 *
 * ========LICENCE========
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License.
 * ========LICENCE========
 */

#ifndef DISABLE_COMMENTATOR
#define DISABLE_COMMENTATOR
#endif

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <iomanip>
#include <cmath>

#include "linbox/ring/modular.h"
#include "linbox/util/commentator.h"
#include "linbox/ring/ntl.h"
#include "linbox/matrix/sparse-matrix.h"
#include "givaro/givtimer.h"

#include "linbox/algorithms/frobenius-large.h"
#include "linbox/algorithms/frobenius-large-search.h"
#include "linbox/algorithms/frobenius-large-bf.h"
#include "linbox/algorithms/frobenius-large-bf-search.h"
#include "linbox/algorithms/frobenius-large-dense.h"
#include "linbox/algorithms/frobenius-large-dense-search.h"
#include "linbox/algorithms/invariant-factors.h"

using namespace LinBox;

template<typename Solver, typename SparseMat, typename PolyRing>
void bench(
    Solver &solver,
    SparseMat &M,
    size_t keff,
    int nruns,
    const PolyRing &R,
    const char *label)
{
    typedef typename PolyRing::Element Polynomial;

    // warmup (not timed)
    std::vector<Polynomial> fs;
    solver.frobeniusInvariants(fs, M, keff);

    double total = 0.0;
    size_t nfactors = 0;
    for (int run = 0; run < nruns; ++run) {
        fs.clear();
        Givaro::Timer T; T.clear(); T.start();
        solver.frobeniusInvariants(fs, M, keff);
        T.stop();
        total += T.usertime();
        nfactors = fs.size();
    }
    std::cout << std::left << std::setw(18) << label
              << "  factors=" << std::setw(6) << nfactors
              << "  avg_time=" << std::fixed << std::setprecision(6) << total/nruns << "s"
              << "  total=" << total << "s" << std::endl;
}

// bit0=Toeplitz, bit1=ToeplitzSearch, bit2=Butterfly, bit3=ButterflySearch,
// bit4=Dense,    bit5=DenseSearch,    bit6=LIFs
int parseAlgoMask(const char *s) {
    if (s == nullptr || strlen(s) == 0) return 127; // all
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
    return mask ? mask : 127;
}

int main(int argc, char **argv) {

    commentator().setMaxDetailLevel(-1);
    commentator().setMaxDepth(-1);
    commentator().setReportStream(std::clog);

    int    p       = 10007;
    int    e       = 1;
    size_t k       = 0;
    int    seed    = time(NULL);
    int    nruns   = 5;
    char   filepath[512] = "";
    char   algoStr[16]   = "";

    static Argument args[] = {
        { 'f', "-f F", "Input matrix file path (SMS format)",
            TYPE_STR, filepath },
        { 'p', "-p P", "Field characteristic",
            TYPE_INT, &p },
        { 'e', "-e E", "Field extension degree",
            TYPE_INT, &e },
        { 'k', "-k K", "Number of invariant factors (0 = log2(n))",
            TYPE_INT, &k },
        { 'r', "-r R", "Random seed",
            TYPE_INT, &seed },
        { 'n', "-n N", "Number of timing runs (first is warmup, not counted)",
            TYPE_INT, &nruns },
        { 'a', "-a A", "Algorithms: 0=Toeplitz 1=ToeplitzSearch 2=Butterfly 3=ButterflySearch 4=Dense 5=DenseSearch 6=LIFs (default: all)",
            TYPE_STR, algoStr },
        END_OF_ARGUMENTS
    };

    parseArguments(argc, argv, args);
    srand(seed);

    if (strlen(filepath) == 0) {
        std::cerr << "Error: no input file specified. Use -f <path/to/matrix.mat>" << std::endl;
        return -1;
    }

    std::ifstream input(filepath);
    if (!input) {
        std::cerr << "Error: could not open " << filepath << std::endl;
        return -1;
    }

    typedef NTL_zz_p  Field;
    typedef NTL_zz_pX PolyRing;
    typedef SparseMatrix<Field> SparseMat;

    Field    F(p, e);
    PolyRing R(F);

    SparseMat M(F);
    M.read(input);
    input.close();

    size_t n   = M.rowdim();
    size_t nnz = M.size();

    size_t keff = (k == 0) ? (size_t)std::max(1.0, std::log2((double)n)) : k;

    int algoMask = parseAlgoMask(algoStr[0] ? algoStr : nullptr);

    std::cout << "=== Frobenius Speed Benchmark ===" << std::endl;
    std::cout << "File:      " << filepath << std::endl;
    std::cout << "n:         " << n        << std::endl;
    std::cout << "nnz:       " << nnz      << std::endl;
    std::cout << "Field:     GF(" << p << "^" << e << ")" << std::endl;
    std::cout << "k:         " << keff << (k == 0 ? " (log2(n) default)" : "") << std::endl;
    std::cout << "seed:      " << seed     << std::endl;
    std::cout << "runs:      " << nruns << " (+ 1 warmup)" << std::endl;
    std::cout << "Algorithms:"
              << ((algoMask &  1) ? " Toeplitz"        : "")
              << ((algoMask &  2) ? " ToeplitzSearch"  : "")
              << ((algoMask &  4) ? " Butterfly"       : "")
              << ((algoMask &  8) ? " ButterflySearch" : "")
              << ((algoMask & 16) ? " Dense"           : "")
              << ((algoMask & 32) ? " DenseSearch"     : "")
              << ((algoMask & 64) ? " LIFs"            : "")
              << std::endl;
    std::cout << std::string(70, '-') << std::endl;

    FrobeniusLarge<PolyRing>                FT(R);
    FrobeniusLargeSearch<PolyRing>          FTS(R);
    FrobeniusLargeButterfly<PolyRing>       FB(R);
    FrobeniusLargeButterflySearch<PolyRing> FBS(R);
    FrobeniusLargeDense<PolyRing>           FD(R);
    FrobeniusLargeDenseSearch<PolyRing>     FDS(R);
    InvariantFactors<Field, PolyRing>       IFD(F, R);

    if (algoMask &  1) bench(FT,  M, keff, nruns, R, "Toeplitz");
    if (algoMask &  2) bench(FTS, M, keff, nruns, R, "ToeplitzSearch");
    if (algoMask &  4) bench(FB,  M, keff, nruns, R, "Butterfly");
    if (algoMask &  8) bench(FBS, M, keff, nruns, R, "ButterflySearch");
    if (algoMask & 16) bench(FD,  M, keff, nruns, R, "Dense");
    if (algoMask & 32) bench(FDS, M, keff, nruns, R, "DenseSearch");
    if (algoMask & 64) bench(IFD, M, keff, nruns, R, "LIFs");

    return 0;
}
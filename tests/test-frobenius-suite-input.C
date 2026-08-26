/* linbox/tests/test-frobenius-suite-input.C
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 * Written by Omesh Dhar Dwivedi <odd23@drexel.edu>
 *
 * Speed benchmark: read a matrix from a given path in SMS format
 * and compare timing of all preconditioner variants + LIFs.
 * No correctness check — pure timing.
 *
 * Usage:
 *   ./test-frobenius-suite-input -f /path/to/matrix.sms [-p 10007] [-k 0] [-r 42] [-a 0123456] [-n 5]
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

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <ctime>

#include "linbox/ring/modular.h"
#include "linbox/util/commentator.h"
#include "linbox/ring/ntl.h"
#include "linbox/matrix/sparse-matrix.h"
#include "givaro/givtimer.h"

#include "linbox/algorithms/frobenius-large-generic.h"
#include "linbox/algorithms/invariant-factors.h"

using namespace LinBox;

template<typename Solver, typename SparseMat, typename PolyRing>
std::vector<typename PolyRing::Element> bench(
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

    // The final timed result is returned only after the timer has stopped.
    return fs;
}

// Compare the requested leading factors. LIFs returns k+2 entries, while
// Dense may omit trailing unit factors after an early exit.
template<typename PolyRing>
bool sameFirstK(
    const std::vector<typename PolyRing::Element> &dense,
    const std::vector<typename PolyRing::Element> &lifs,
    size_t k,
    const PolyRing &R,
    size_t &badIndex)
{
    if (lifs.size() < k) {
        badIndex = lifs.size();
        return false;
    }

    for (size_t i = 0; i < k; ++i) {
        if (i < dense.size()) {
            // LinBox's Frobenius tests compare polynomial elements directly.
            if (dense[i] != lifs[i]) {
                badIndex = i;
                return false;
            }
        } else if (!R.isOne(lifs[i])) {
            // A missing Dense entry denotes a trailing unit invariant factor.
            badIndex = i;
            return false;
        }
    }

    badIndex = k;
    return true;
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

    int         p       = 10007;
    int         e       = 1;
    size_t      k       = 0;
    int         seed    = time(NULL);
    int         nruns   = 5;
    std::string filepath;
    std::string algoStr;

    static Argument args[] = {
        { 'f', "-f F", "Input matrix file path (SMS format)",
            TYPE_STR, &filepath },
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
            TYPE_STR, &algoStr },
        END_OF_ARGUMENTS
    };

    parseArguments(argc, argv, args);
    srand(seed);

    if (filepath.empty()) {
        std::cerr << "Error: no input file specified. Use -f <path/to/matrix.sms>" << std::endl;
        return -1;
    }

    std::ifstream input(filepath.c_str());
    if (!input) {
        std::cerr << "Error: could not open " << filepath << std::endl;
        return -1;
    }

    typedef NTL_zz_p  Field;
    typedef NTL_zz_pX PolyRing;
    typedef SparseMatrix<Field> SparseMat;
    typedef PolyRing::Element Polynomial;

    Field    F(p, e);
    PolyRing R(F);

    SparseMat M(F);
    M.read(input);
    input.close();

    size_t n   = M.rowdim();
    size_t nnz = M.size();

    size_t keff = (k == 0) ? (size_t)std::max(1.0, std::log2((double)n)) : k;

    int algoMask = parseAlgoMask(algoStr.empty() ? nullptr : algoStr.c_str());

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

    std::vector<Polynomial> denseResult;
    std::vector<Polynomial> lifsResult;

    if (algoMask &  1) bench(FT,  M, keff, nruns, R, "Toeplitz");
    if (algoMask &  2) bench(FTS, M, keff, nruns, R, "ToeplitzSearch");
    if (algoMask &  4) bench(FB,  M, keff, nruns, R, "Butterfly");
    if (algoMask &  8) bench(FBS, M, keff, nruns, R, "ButterflySearch");
    if (algoMask & 16) denseResult = bench(FD,  M, keff, nruns, R, "Dense");
    if (algoMask & 32) bench(FDS, M, keff, nruns, R, "DenseSearch");
    if (algoMask & 64) lifsResult = bench(IFD, M, keff, nruns, R, "LIFs");

    // Optional untimed cross-check. It runs only when both Dense and LIFs
    // were requested together (for example, with -a 46).
    if ((algoMask & 16) && (algoMask & 64)) {
        size_t badIndex = 0;
        if (!sameFirstK(denseResult, lifsResult, keff, R, badIndex)) {
            std::cerr << "Correctness check: MISMATCH at invariant factor "
                      << (badIndex + 1)
                      << " (Dense returned " << denseResult.size()
                      << " entries; LIFs returned " << lifsResult.size()
                      << ")" << std::endl;

            if (badIndex < denseResult.size())
                R.write(std::cerr << "Dense: ", denseResult[badIndex]) << std::endl;
            else
                std::cerr << "Dense: 1 (implicit trailing unit)" << std::endl;

            if (badIndex < lifsResult.size())
                R.write(std::cerr << "LIFs:  ", lifsResult[badIndex]) << std::endl;
            else
                std::cerr << "LIFs:  <missing>" << std::endl;

            std::cerr << "This compares one randomized trial from each method; "
                      << "rerun before diagnosing an implementation bug."
                      << std::endl;
            return 2;
        }

        std::cout << "Correctness check: PASS (first " << keff
                  << " invariant factors; one randomized trial)"
                  << std::endl;
    }

    return 0;
}
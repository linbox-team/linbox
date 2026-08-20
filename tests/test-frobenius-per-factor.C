/* linbox/tests/test-frobenius-per-factor.C
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 * Written by Omesh Dhar Dwivedi <odd23@drexel.edu>
 *
 * For a given matrix, times kthInvariantFactor directly (no search overhead)
 * for each k=1,...,kmax and prints a per-factor timing table.
 *
 * Usage:
 *   ./test-frobenius-per-factor -f /path/to/matrix.mat [-p 10007] [-k 0] [-n 3] [-a 12]
 *
 * -k: max number of factors to test (0 = log2(n))
 * -n: number of timing runs per (algorithm, k) pair (averaged)
 * -a: 0=Toeplitz 1=Butterfly 2=Dense (default: 12 = Butterfly + Dense)
 */

#ifndef DISABLE_COMMENTATOR
#define DISABLE_COMMENTATOR
#endif

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstdlib>

#include "linbox/ring/modular.h"
#include "linbox/util/commentator.h"
#include "linbox/ring/ntl.h"
#include "linbox/matrix/sparse-matrix.h"
#include "givaro/givtimer.h"

#include "linbox/algorithms/frobenius-large-generic.h"

using namespace LinBox;

typedef NTL_zz_p  Field;
typedef NTL_zz_pX PolyRing;
typedef SparseMatrix<Field> SparseMat;
typedef PolyRing::Element Polynomial;

template<class Solver>
double timeKth(Solver &solver, const SparseMat &M, const Polynomial &f1, size_t k, int nruns) {
    double total = 0.0;
    for (int run = 0; run < nruns; ++run) {
        Polynomial fk;
        size_t Ck = 0;
        Givaro::Timer T; T.clear(); T.start();
        solver.kthInvariantFactor(fk, Ck, M, f1, k);
        T.stop();
        total += T.usertime();
    }
    return total / nruns;
}

int parseAlgoMask(const char *s) {
    if (s == nullptr || strlen(s) == 0) return 6; // Butterfly + Dense
    int mask = 0;
    for (size_t i = 0; i < strlen(s); ++i) {
        if (s[i] == '0') mask |= 1;
        if (s[i] == '1') mask |= 2;
        if (s[i] == '2') mask |= 4;
    }
    return mask ? mask : 6;
}

int main(int argc, char **argv) {

    commentator().setMaxDetailLevel(-1);
    commentator().setMaxDepth(-1);
    commentator().setReportStream(std::clog);

    int    p        = 10007;
    int    e        = 1;
    size_t k        = 0;
    int    seed     = time(NULL);
    int    nruns    = 3;
    std::string filepath;
    std::string algoStr;

    static Argument args[] = {
        { 'f', "-f F", "Input matrix file path (SMS format)",        TYPE_STR, &filepath },
        { 'p', "-p P", "Field characteristic",                       TYPE_INT, &p },
        { 'e', "-e E", "Field extension degree",                     TYPE_INT, &e },
        { 'k', "-k K", "Max k to test (0 = log2(n))",               TYPE_INT, &k },
        { 'r', "-r R", "Random seed",                                TYPE_INT, &seed },
        { 'n', "-n N", "Runs per (algo, k) pair",                   TYPE_INT, &nruns },
        { 'a', "-a A", "Algorithms: 0=Toeplitz 1=Butterfly 2=Dense (default: 12)",
            TYPE_STR, &algoStr },
        END_OF_ARGUMENTS
    };

    parseArguments(argc, argv, args);
    srand(seed);

    if (filepath.empty()) {
        std::cerr << "Error: no input file. Use -f <path>" << std::endl;
        return -1;
    }
    std::ifstream input(filepath.c_str());
    if (!input) {
        std::cerr << "Error: could not open " << filepath << std::endl;
        return -1;
    }

    Field    F(p, e);
    PolyRing R(F);

    SparseMat M(F);
    M.read(input);
    input.close();

    size_t n   = M.rowdim();
    size_t nnz = M.size();
    size_t kmax = (k == 0) ? (size_t)std::max(1.0, std::log2((double)n)) : k;

    int algoMask = parseAlgoMask(algoStr.empty() ? nullptr : algoStr.c_str());

    std::cout << "=== Per-Factor Timing ===" << std::endl;
    std::cout << "File:      " << filepath << std::endl;
    std::cout << "n:         " << n << std::endl;
    std::cout << "nnz:       " << nnz << std::endl;
    std::cout << "Field:     GF(" << p << "^" << e << ")" << std::endl;
    std::cout << "kmax:      " << kmax << (k == 0 ? " (log2(n) default)" : "") << std::endl;
    std::cout << "runs/cell: " << nruns << std::endl;
    std::cout << "Algorithms:"
              << ((algoMask & 1) ? " Toeplitz"  : "")
              << ((algoMask & 2) ? " Butterfly" : "")
              << ((algoMask & 4) ? " Dense"     : "")
              << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    FrobeniusLarge<PolyRing>          FT(R);
    FrobeniusLargeButterfly<PolyRing> FB(R);
    FrobeniusLargeDense<PolyRing>     FD(R);

    // compute f1 = minpoly(A) once, reused as modulus for all k
    Polynomial f1;
    FT.minpoly(f1, M);

    // header
    std::cout << std::left << std::setw(6) << "k";
    if (algoMask & 1) std::cout << std::setw(14) << "Toeplitz";
    if (algoMask & 2) std::cout << std::setw(14) << "Butterfly";
    if (algoMask & 4) std::cout << std::setw(14) << "Dense";
    std::cout << "Winner\n";
    std::cout << std::string(60, '-') << "\n";

    for (size_t ki = 1; ki <= kmax; ++ki) {
        double t[3] = {-1.0, -1.0, -1.0};

        if (algoMask & 1) t[0] = timeKth(FT, M, f1, ki, nruns);
        if (algoMask & 2) t[1] = timeKth(FB, M, f1, ki, nruns);
        if (algoMask & 4) t[2] = timeKth(FD, M, f1, ki, nruns);

        const char* labels[3] = {"Toeplitz", "Butterfly", "Dense"};
        double best = 1e18;
        int winner = -1;
        for (int i = 0; i < 3; ++i)
            if (t[i] >= 0 && t[i] < best) { best = t[i]; winner = i; }

        std::cout << std::left << std::setw(6) << ki;
        for (int i = 0; i < 3; ++i)
            if (t[i] >= 0)
                std::cout << std::setw(14) << std::fixed << std::setprecision(6) << t[i];
        std::cout << (winner >= 0 ? labels[winner] : "?") << "\n";
        std::cout.flush();
    }

    return 0;
}
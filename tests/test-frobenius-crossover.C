/* linbox/tests/test-frobenius-crossover.C
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 * Written by Omesh Dhar Dwivedi <odd23@drexel.edu>
 *
 * Constructs synthetic Jordan matrices (like test-frobenius-suite) and finds
 * the crossover k where Butterfly.kthInvariantFactor starts beating Dense,
 * across a sweep of matrix sizes and shapes.
 *
 * Usage:
 *   ./test-frobenius-crossover [-L smin] [-U smax] [-p P] [-n N] [-t T]
 *
 * -L/-U: range of s values (matrix size n = w*h*s*(s+1)/2)
 * -n N:  timing runs per (algo, k) pair
 * -t T:  max k to probe (0 = 3*log2(n))
 * -p P:  field characteristic
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

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>

#include "linbox/ring/modular.h"
#include "linbox/util/commentator.h"
#include "linbox/ring/ntl.h"
#include "givaro/givtimer.h"

#include "linbox/algorithms/frobenius-large-generic.h"
#include "test-frobenius-suite.h"

using namespace LinBox;

typedef NTL_zz_p  Field;
typedef NTL_zz_pX PolyRing;
typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;
typedef PolyRing::Element Polynomial;

template<class Solver>
double timeKth(Solver &solver, const SparseMat &M, const Polynomial &f1,
               size_t k, int nruns)
{
    double total = 0.0;
    for (int r = 0; r < nruns; ++r) {
        Polynomial fk;
        size_t rawDegree = 0;
        Givaro::Timer T; T.clear(); T.start();
        solver.kthInvariantFactor(fk, rawDegree, M, f1, k);
        T.stop();
        total += T.usertime();
    }
    return total / nruns;
}

struct CrossoverResult {
    size_t k;           // first k where butterfly < dense (0 = not found)
    double t_butterfly;
    double t_dense;
    double log2n;
};

CrossoverResult findCrossover(
    FrobeniusLargeButterfly<PolyRing> &FB,
    FrobeniusLargeDense<PolyRing>     &FD,
    const SparseMat &M,
    const Polynomial &f1,
    size_t kmax,
    int nruns)
{
    CrossoverResult res = {0, 0.0, 0.0, std::log2((double)M.rowdim())};
    for (size_t k = 1; k <= kmax; ++k) {
        double tb = timeKth(FB, M, f1, k, nruns);
        double td = timeKth(FD, M, f1, k, nruns);
        if (tb < td) {
            res.k = k; res.t_butterfly = tb; res.t_dense = td;
            return res;
        }
    }
    return res;
}

int main(int argc, char **argv) {

    commentator().setMaxDetailLevel(-1);
    commentator().setMaxDepth(-1);
    commentator().setReportStream(std::clog);

    int    p        = 10000019;
    int    e        = 1;
    int    seed     = time(NULL);
    size_t s_min    = 3;
    size_t s_max    = 20;
    int    nruns    = 3;
    size_t kmax_arg = 0;

    static Argument args[] = {
        { 'p', "-p P", "Field characteristic",                          TYPE_INT, &p },
        { 'e', "-e E", "Field extension degree",                        TYPE_INT, &e },
        { 'r', "-r R", "Random seed",                                   TYPE_INT, &seed },
        { 'L', "-L L", "Min s (matrix size n = w*h*s*(s+1)/2)",        TYPE_INT, &s_min },
        { 'U', "-U U", "Max s",                                         TYPE_INT, &s_max },
        { 'n', "-n N", "Timing runs per (algo, k) pair",                TYPE_INT, &nruns },
        { 't', "-t T", "Max k to probe (0 = 3*log2(n))",               TYPE_INT, &kmax_arg },
        END_OF_ARGUMENTS
    };

    parseArguments(argc, argv, args);
    srand(seed);

    Field    F(p, e);
    PolyRing R(F);

    FrobeniusLargeButterfly<PolyRing> FB(R);
    FrobeniusLargeDense<PolyRing>     FD(R);

    // Same shapes as test-frobenius-suite defaults, swept over s_min..s_max
    struct Shape { size_t w; size_t h; int poly; std::string name; };
    std::vector<Shape> shapes = {
        {1, 10, 0, "tall_x"},
        {5,  1, 1, "flat_xm1"},
        {4,  4, 0, "tri_x"},
        {1, 10, 2, "tall_xp1"},
        {5,  1, 0, "flat_x"},
    };

    std::cout << "=== Butterfly vs Dense: Crossover k by Matrix Type and Size ===" << std::endl;
    std::cout << "Field: GF(" << p << "^" << e << ")  seed=" << seed
              << "  s=" << s_min << ".." << s_max
              << "  runs/cell=" << nruns << std::endl;
    std::cout << std::string(98, '-') << std::endl;
    std::cout << std::left
              << std::setw(14) << "Type"
              << std::setw(5)  << "s"
              << std::setw(8)  << "n"
              << std::setw(8)  << "nnz"
              << std::setw(9)  << "log2(n)"
              << std::setw(12) << "Crossover k"
              << std::setw(14) << "Butterfly(s)"
              << std::setw(12) << "Dense(s)"
              << std::setw(10) << "Dense/Bfly"
              << "\n";
    std::cout << std::string(98, '-') << "\n";
    std::cout.flush();

    for (auto &shape : shapes) {
        bool first = true;
        for (size_t s = s_min; s <= s_max; ++s) {
            size_t n = shape.w * shape.h * s * (s + 1) / 2;
            if (n < 2) continue;

            size_t kmax = (kmax_arg > 0)
                ? kmax_arg
                : (size_t)std::max(3.0, 3.0 * std::log2((double)n));

            SparseMat M(F, n, n);
            std::vector<Polynomial> expected;
            buildJordanMatrix(M, expected, F, R, s, shape.w, shape.h, shape.poly);
            sprayMatrix(M, F, seed);

            Polynomial f1;
            FD.minpoly(f1, M);

            CrossoverResult cr = findCrossover(FB, FD, M, f1, kmax, nruns);

            std::cout << std::left
                      << std::setw(14) << (first ? shape.name : "")
                      << std::setw(5)  << s
                      << std::setw(8)  << n
                      << std::setw(8)  << M.size()
                      << std::setw(9)  << std::fixed << std::setprecision(1) << cr.log2n;
            first = false;

            if (cr.k == 0) {
                std::cout << std::setw(12) << ">"+std::to_string(kmax)
                          << "  (Dense wins throughout)\n";
            } else {
                std::cout << std::setw(12) << cr.k
                          << std::setw(14) << std::fixed << std::setprecision(6) << cr.t_butterfly
                          << std::setw(12) << std::fixed << std::setprecision(6) << cr.t_dense
                          << std::setw(10) << std::fixed << std::setprecision(2) << cr.t_dense/cr.t_butterfly
                          << "x\n";
            }
            std::cout.flush();
        }
        std::cout << std::string(98, '-') << "\n";
    }

    return 0;
}
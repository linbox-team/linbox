/* linbox/tests/test-frobenius-suite.h
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

#ifndef __LINBOX_test_frobenius_suite_H
#define __LINBOX_test_frobenius_suite_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "linbox/matrix/sparse-matrix.h"
#include "givaro/givtimer.h"

using namespace LinBox;

// ============================================================
// Matrix shape descriptor
// ============================================================

struct TestParams {
    size_t s;
    size_t w;
    size_t h;
    int    poly;
    std::string name;
};

// Compute n = w*h*s*(s+1)/2 without overflowing size_t.
inline bool checkedCaseOrder(size_t s, size_t w, size_t h, size_t &n)
{
    if (s == 0 || w == 0 || h == 0) return false;

    size_t a = s;
    size_t b = s + 1;
    if (b == 0) return false;
    if ((a & 1U) == 0) a /= 2;
    else                b /= 2;

    const size_t maxSize = std::numeric_limits<size_t>::max();
    if (a != 0 && b > maxSize / a) return false;
    size_t value = a * b;
    if (w > maxSize / value) return false;
    value *= w;
    if (h > maxSize / value) return false;
    n = value * h;
    return true;
}

// ============================================================
// Polynomial helpers
// ============================================================

template<class Field>
typename Field::Element evalPoly(const Field &F, int poly, typename Field::Element val) {
    typedef typename Field::Element Element;
    Element res;
    switch (poly) {
        case 0: return val;
        case 1: { Element one; F.init(one, 1); F.sub(res, val, one); return res; }
        case 2: { Element one; F.init(one, 1); F.add(res, val, one); return res; }
        default: return val;
    }
}

// ============================================================
// Jordan block construction
// ============================================================
template<class Field>
void writeJordanBlock(
    SparseMatrix<Field, SparseMatrixFormat::CSR> &M,
    const Field &F,
    int poly,
    size_t sz,
    size_t row0,
    size_t col0)
{
    typedef typename Field::Element Element;
    if (poly != 0) {
        Element lambda;
        if (poly == 1) F.init(lambda, 1);
        else           F.init(lambda, -1);
        for (size_t i = 0; i < sz; ++i)
            M.setEntry(row0 + i, col0 + i, lambda);
    }
    Element one; F.init(one, 1);
    for (size_t i = 0; i + 1 < sz; ++i)
        M.setEntry(row0 + i + 1, col0 + i, one);
}

// ============================================================
// Build Jordan form matrix
// ============================================================
template<class Field, class PolyRing>
bool buildJordanMatrix(
    SparseMatrix<Field, SparseMatrixFormat::CSR> &M,
    std::vector<typename PolyRing::Element> &expectedFactors,
    const Field &F,
    const PolyRing &R,
    size_t s,
    size_t w,
    size_t h,
    int poly)
{
    typedef typename PolyRing::Element Polynomial;
    typedef typename PolyRing::Coeff   Coeff;

    size_t n = 0;
    if (!checkedCaseOrder(s, w, h, n)) return false;

    M.resize(n, n);
    expectedFactors.clear();

    size_t row0 = 0;
    for (size_t lev = s; lev >= 1; --lev) {
        size_t bsz = lev * h;
        for (size_t rep = 0; rep < w; ++rep) {
            writeJordanBlock(M, F, poly, bsz, row0, row0);
            row0 += bsz;

            Polynomial base, factor;
            R.assign(base, R.zero);
            switch (poly) {
                case 0: R.setCoeff(base, 0, (Coeff)0);  R.setCoeff(base, 1, (Coeff)1);  break;
                case 1: R.setCoeff(base, 0, (Coeff)-1); R.setCoeff(base, 1, (Coeff)1);  break;
                case 2: R.setCoeff(base, 0, (Coeff)1);  R.setCoeff(base, 1, (Coeff)1);  break;
            }
            R.assign(factor, R.one);
            for (size_t i = 0; i < bsz; ++i)
                R.mulin(factor, base);

            expectedFactors.push_back(factor);
        }
    }
    M.finalize();
    return true;
}

// ============================================================
// Spray: random elementary similarity transforms
// ============================================================
template<class Field>
void sprayMatrix(
    SparseMatrix<Field, SparseMatrixFormat::CSR> &M,
    const Field &F,
    int seed)
{
    typedef typename Field::Element Element;
    typedef typename Field::RandIter RandIter;

    size_t n       = M.rowdim();
    // There is no distinct j for a 1x1 matrix.
    if (n < 2) return;

    size_t nSprays = (size_t)std::ceil(std::sqrt((double)n));
    RandIter RI(F, 0, seed);

    for (size_t spray = 0; spray < nSprays; ++spray) {
        size_t i = (size_t)(std::rand() % n);
        size_t j;
        do { j = (size_t)(std::rand() % n); } while (j == i);

        Element alpha;
        do { RI.random(alpha); } while (F.isZero(alpha));

        std::vector<std::pair<size_t,Element>> rowJ;
        for (size_t col = 0; col < n; ++col) {
            Element val; F.init(val, 0);
            M.getEntry(val, j, col);
            if (!F.isZero(val)) rowJ.push_back({col, val});
        }
        for (typename std::vector<std::pair<size_t, Element> >::const_iterator
                 it = rowJ.begin(); it != rowJ.end(); ++it) {
            const size_t col = it->first;
            const Element &val = it->second;
            Element cur; F.init(cur, 0);
            M.getEntry(cur, i, col);
            Element contrib; F.mul(contrib, alpha, val);
            F.addin(cur, contrib);
            M.setEntry(i, col, cur);
        }

        std::vector<std::pair<size_t,Element>> colI;
        for (size_t row = 0; row < n; ++row) {
            Element val; F.init(val, 0);
            M.getEntry(val, row, i);
            if (!F.isZero(val)) colI.push_back({row, val});
        }
        Element negAlpha; F.neg(negAlpha, alpha);
        for (typename std::vector<std::pair<size_t, Element> >::const_iterator
                 it = colI.begin(); it != colI.end(); ++it) {
            const size_t row = it->first;
            const Element &val = it->second;
            Element cur; F.init(cur, 0);
            M.getEntry(cur, row, j);
            Element contrib; F.mul(contrib, negAlpha, val);
            F.addin(cur, contrib);
            M.setEntry(row, j, cur);
        }
    }
    M.finalize();
}

// ============================================================
// Timing and result records
// ============================================================
struct RunTiming {
    double wall;
    double user;
    double system;
    size_t returnedEntries;
};

struct TimingTotals {
    double wall;
    double user;
    double system;
};

struct TestResult {
    std::string algoName;
    std::string caseName;
    size_t      n;
    size_t      nnz;
    std::vector<RunTiming> timings;
    double      walltime;
    double      usertime;  // Legacy field: average user CPU time.
    double      systime;
    bool        pass;
};

inline TimingTotals timingTotals(const std::vector<RunTiming> &timings)
{
    TimingTotals totals = { 0.0, 0.0, 0.0 };
    for (size_t i = 0; i < timings.size(); ++i) {
        totals.wall += timings[i].wall;
        totals.user += timings[i].user;
        totals.system += timings[i].system;
    }
    return totals;
}

inline void printRunTiming(
    const std::string &algoName,
    const std::string &caseName,
    int run,
    int nruns,
    const RunTiming &timing)
{
    std::cout << "Timing"
              << "  algorithm=" << algoName
              << "  case=" << caseName
              << "  run=" << (run + 1) << "/" << nruns
              << "  wall_time=" << std::fixed << std::setprecision(6)
              << timing.wall << "s"
              << "  user_time=" << timing.user << "s"
              << "  system_time=" << timing.system << "s"
              << "  returned_entries=" << timing.returnedEntries
              << std::endl;
}

inline void printTimingSummary(const TestResult &result)
{
    const TimingTotals totals = timingTotals(result.timings);
    const double runs = static_cast<double>(result.timings.size());

    std::cout << "Timing"
              << "  algorithm=" << result.algoName
              << "  case=" << result.caseName
              << "  summary"
              << "  runs=" << result.timings.size()
              << "  avg_wall_time=" << std::fixed << std::setprecision(6)
              << totals.wall / runs << "s"
              << "  total_wall_time=" << totals.wall << "s"
              << "  avg_user_time=" << totals.user / runs << "s"
              << "  total_user_time=" << totals.user << "s"
              << "  avg_system_time=" << totals.system / runs << "s"
              << "  total_system_time=" << totals.system << "s"
              << std::endl;
}

template<class PolyRing>
bool checkComputedFactors(
    const std::vector<typename PolyRing::Element> &computed,
    const std::vector<typename PolyRing::Element> &expected,
    const PolyRing &R,
    const std::string &algoName,
    const std::string &caseName,
    size_t k,
    size_t matrixOrder,
    int run,
    int nruns)
{
    // expected stores the non-unit prefix. Both a missing solver entry and a
    // position beyond that prefix denote an implicit trailing unit factor.
    const size_t nCheck = (k > 0) ? k : matrixOrder;
    bool pass = true;
    size_t badIndex = nCheck;
    for (size_t i = 0; i < nCheck; ++i) {
        const bool expectedExplicit = i < expected.size();
        const bool computedExplicit = i < computed.size();

        if (expectedExplicit) {
            if (!computedExplicit || expected[i] != computed[i]) {
                pass = false;
                badIndex = i;
                break;
            }
        } else if (computedExplicit && !R.isOne(computed[i])) {
            pass = false;
            badIndex = i;
            break;
        }
    }

    if (!pass) {
        std::cout << "\n[FAIL] " << algoName << " on " << caseName
                  << " (run " << (run + 1) << "/" << nruns << ")\n";
        std::cout << "  requested positions: " << nCheck
                  << ", computed entries: " << computed.size()
                  << ", expected non-unit entries: " << expected.size()
                  << "\n";
        std::cout << "  first mismatch [" << badIndex << "]:\n";
        std::cout << "    computed: ";
        if (badIndex < computed.size()) R.write(std::cout, computed[badIndex]);
        else                            std::cout << "1 (implicit)";
        std::cout << "\n    expected: ";
        if (badIndex < expected.size()) R.write(std::cout, expected[badIndex]);
        else                            std::cout << "1 (implicit)";
        std::cout << "\n";
    }
    return pass;
}

// ============================================================
// Run one algorithm on one test case (no warmup calls)
// ============================================================
template<class FrobeniusObject, class Field, class PolyRing>
TestResult runOne(
    FrobeniusObject &FO,
    const std::string &algoName,
    const Field &F,
    const PolyRing &R,
    const SparseMatrix<Field, SparseMatrixFormat::CSR> &M,
    const std::vector<typename PolyRing::Element> &expected,
    const std::string &caseName,
    size_t k,
    int nruns = 1)
{
    typedef typename PolyRing::Element Polynomial;
    (void)F;

    TestResult result;
    result.algoName = algoName;
    result.caseName = caseName;
    result.n = M.rowdim();
    result.nnz = M.size();
    result.pass = true;
    result.timings.reserve(static_cast<size_t>(nruns));

    std::vector<Polynomial> computed;
    for (int run = 0; run < nruns; ++run) {
        computed.clear();
        Givaro::Timer timer;
        timer.clear();
        timer.start();
        FO.frobeniusInvariants(computed, M, k);
        timer.stop();

        RunTiming timing;
        timing.wall = timer.realtime();
        timing.user = timer.usertime();
        timing.system = timer.systime();
        timing.returnedEntries = computed.size();
        result.timings.push_back(timing);
        printRunTiming(algoName, caseName, run, nruns, timing);

        if (!checkComputedFactors(computed, expected, R, algoName,
                                  caseName, k, M.rowdim(), run, nruns))
            result.pass = false;
    }

    const TimingTotals totals = timingTotals(result.timings);
    result.walltime = totals.wall / nruns;
    result.usertime = totals.user / nruns;
    result.systime = totals.system / nruns;
    printTimingSummary(result);
    return result;
}

// LIFs has a separate path because frobeniusInvariants() hardcodes b=k+2.
// The explicit-b API returns the opposite order, so the reverse is timed too.
template<class FrobeniusObject, class Field, class PolyRing>
TestResult runOneLIFs(
    FrobeniusObject &FO,
    const std::string &algoName,
    const Field &F,
    const PolyRing &R,
    const SparseMatrix<Field, SparseMatrixFormat::CSR> &M,
    const std::vector<typename PolyRing::Element> &expected,
    const std::string &caseName,
    size_t k,
    size_t blockSize,
    int nruns)
{
    typedef typename PolyRing::Element Polynomial;
    (void)F;

    TestResult result;
    result.algoName = algoName;
    result.caseName = caseName;
    result.n = M.rowdim();
    result.nnz = M.size();
    result.pass = true;
    result.timings.reserve(static_cast<size_t>(nruns));

    std::cout << "LIFs"
              << "  case=" << caseName
              << "  requested_k=";
    if (k == 0) std::cout << "all";
    else        std::cout << k;
    std::cout << "  block_b=" << blockSize << std::endl;

    std::vector<Polynomial> computed;
    for (int run = 0; run < nruns; ++run) {
        computed.clear();
        Givaro::Timer timer;
        timer.clear();
        timer.start();
        FO.largestInvariantFactors(computed, M, blockSize);
        std::reverse(computed.begin(), computed.end());
        timer.stop();

        RunTiming timing;
        timing.wall = timer.realtime();
        timing.user = timer.usertime();
        timing.system = timer.systime();
        timing.returnedEntries = computed.size();
        result.timings.push_back(timing);
        printRunTiming(algoName, caseName, run, nruns, timing);

        if (!checkComputedFactors(computed, expected, R, algoName,
                                  caseName, k, M.rowdim(), run, nruns))
            result.pass = false;
    }

    const TimingTotals totals = timingTotals(result.timings);
    result.walltime = totals.wall / nruns;
    result.usertime = totals.user / nruns;
    result.systime = totals.system / nruns;
    printTimingSummary(result);
    return result;
}

// ============================================================
// Print aggregate results table (wall time is the primary runtime)
// ============================================================
inline void printTable(const std::vector<TestResult> &results)
{
    std::cout << "\n";
    std::cout << std::left
              << std::setw(29) << "Algorithm"
              << std::setw(26) << "Test Case"
              << std::setw(9)  << "n"
              << std::setw(11) << "nnz"
              << std::setw(7)  << "Runs"
              << std::setw(14) << "Avg wall (s)"
              << std::setw(14) << "Avg user (s)"
              << std::setw(14) << "Avg sys (s)"
              << std::setw(8)  << "Pass"
              << "\n";
    std::cout << std::string(132, '-') << "\n";

    for (size_t i = 0; i < results.size(); ++i) {
        const TestResult &result = results[i];
        const TimingTotals totals = timingTotals(result.timings);
        const double runs = static_cast<double>(result.timings.size());

        std::cout << std::left
                  << std::setw(29) << result.algoName
                  << std::setw(26) << result.caseName
                  << std::setw(9)  << result.n
                  << std::setw(11) << result.nnz
                  << std::setw(7)  << result.timings.size()
                  << std::setw(14) << std::fixed << std::setprecision(6)
                  << totals.wall / runs
                  << std::setw(14) << totals.user / runs
                  << std::setw(14) << totals.system / runs
                  << std::setw(8)  << (result.pass ? "PASS" : "FAIL")
                  << "\n";
    }
    std::cout << "\n";
}

// ============================================================
// Default test cases
// ============================================================
inline std::vector<TestParams> defaultTestCases() {
    return {
        {3, 1, 10, 0, "tall_x"},
        {3, 5,  1, 1, "flat_xm1"},
        {4, 4,  4, 0, "tri_x"},
        {3, 1, 10, 2, "tall_xp1"},
        {3, 5,  1, 0, "flat_x"},
    };
}

// ============================================================
// Main test suite runner
// Algorithms bitmask:
//   bit0 = Toeplitz
//   bit1 = ToeplitzSearch
//   bit2 = Butterfly
//   bit3 = ButterflySearch
//   bit4 = Dense
//   bit5 = DenseSearch
//   bit6 = LIFs
//   bit7 = ToeplitzFactorAware
//   bit8 = ButterflyFactorAware
//   bit9 = DenseFactorAware
// ============================================================
template<
    class FrobeniusToeplitz,
    class FrobeniusToeplitzSearch,
    class FrobeniusButterfly,
    class FrobeniusButterflySearch,
    class FrobeniusDense,
    class FrobeniusDenseSearch,
    class FrobeniusLifs,
    class FrobeniusToeplitzFactorAware,
    class FrobeniusButterflyFactorAware,
    class FrobeniusDenseFactorAware,
    class Field,
    class PolyRing>
bool runSuite(
    FrobeniusToeplitz        &FT,
    FrobeniusToeplitzSearch  &FTS,
    FrobeniusButterfly       &FB,
    FrobeniusButterflySearch &FBS,
    FrobeniusDense           &FD,
    FrobeniusDenseSearch     &FDS,
    FrobeniusLifs            &IFD,
    FrobeniusToeplitzFactorAware  &FTF,
    FrobeniusButterflyFactorAware &FBF,
    FrobeniusDenseFactorAware     &FDF,
    const Field              &F,
    const PolyRing           &R,
    size_t k,
    int    algoMask,
    int    seed,
    int    nruns,
    size_t lifsBlock,
    size_t custom_s    = 0,
    size_t custom_w    = 0,
    size_t custom_h    = 0,
    int    custom_poly = -1)
{
    typedef typename PolyRing::Element Polynomial;
    typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;

    std::vector<TestParams> cases = defaultTestCases();
    if (custom_s > 0 && custom_w > 0 && custom_h > 0) {
        int p = (custom_poly >= 0) ? custom_poly : 0;
        std::string name = "custom_s" + std::to_string(custom_s)
                         + "_w" + std::to_string(custom_w)
                         + "_h" + std::to_string(custom_h);
        cases.insert(cases.begin(), {custom_s, custom_w, custom_h, p, name});
    }

    if (nruns <= 0) {
        std::cerr << "Error: number of timed runs must be positive."
                  << std::endl;
        return false;
    }

    // Preflight every case before constructing a matrix or running a solver.
    std::vector<size_t> caseOrders;
    caseOrders.reserve(cases.size());
    for (size_t i = 0; i < cases.size(); ++i) {
        size_t n = 0;
        if (!checkedCaseOrder(cases[i].s, cases[i].w, cases[i].h, n)) {
            std::cerr << "Error: invalid or overflowing dimensions for case "
                      << cases[i].name << "." << std::endl;
            return false;
        }
        if (k > n) {
            std::cerr << "Error: requested k=" << k << " exceeds n=" << n
                      << " for case " << cases[i].name << "." << std::endl;
            return false;
        }
        if ((algoMask & 64) && k > 0 &&
            (lifsBlock < k || lifsBlock > n)) {
            std::cerr << "Error: LIFs block size b=" << lifsBlock
                      << " must satisfy k <= b <= n for case "
                      << cases[i].name << " (k=" << k << ", n=" << n
                      << ")." << std::endl;
            return false;
        }
        caseOrders.push_back(n);
    }

    std::vector<TestResult> results;
    bool allPass = true;

    for (size_t caseIndex = 0; caseIndex < cases.size(); ++caseIndex) {
        const TestParams &tc = cases[caseIndex];
        const size_t n = caseOrders[caseIndex];
        SparseMat M(F, n, n);
        std::vector<Polynomial> expected;
        if (!buildJordanMatrix(M, expected, F, R, tc.s, tc.w, tc.h,
                               tc.poly)) {
            std::cerr << "Error: could not construct case " << tc.name
                      << "." << std::endl;
            return false;
        }
        sprayMatrix(M, F, seed);

        if (algoMask & 1) {
            TestResult r = runOne(FT, "Toeplitz", F, R, M, expected,
                                  tc.name, k, nruns);
            allPass = r.pass && allPass; results.push_back(r);
        }
        if (algoMask & 2) {
            TestResult r = runOne(FTS, "ToeplitzSearch", F, R, M,
                                  expected, tc.name, k, nruns);
            allPass = r.pass && allPass; results.push_back(r);
        }
        if (algoMask & 4) {
            TestResult r = runOne(FB, "Butterfly", F, R, M, expected,
                                  tc.name, k, nruns);
            allPass = r.pass && allPass; results.push_back(r);
        }
        if (algoMask & 8) {
            TestResult r = runOne(FBS, "ButterflySearch", F, R, M,
                                  expected, tc.name, k, nruns);
            allPass = r.pass && allPass; results.push_back(r);
        }
        if (algoMask & 16) {
            TestResult r = runOne(FD, "Dense", F, R, M, expected,
                                  tc.name, k, nruns);
            allPass = r.pass && allPass; results.push_back(r);
        }
        if (algoMask & 32) {
            TestResult r = runOne(FDS, "DenseSearch", F, R, M,
                                  expected, tc.name, k, nruns);
            allPass = r.pass && allPass; results.push_back(r);
        }
        if (algoMask & 64) {
            const size_t caseBlock = (k == 0) ? n : lifsBlock;
            TestResult r = runOneLIFs(IFD, "LIFs", F, R, M, expected,
                                      tc.name, k, caseBlock, nruns);
            allPass = r.pass && allPass; results.push_back(r);
        }
        if (algoMask & 128) {
            TestResult r = runOne(FTF, "ToeplitzFactorAware", F, R, M,
                                  expected, tc.name, k, nruns);
            allPass = r.pass && allPass; results.push_back(r);
        }
        if (algoMask & 256) {
            TestResult r = runOne(FBF, "ButterflyFactorAware", F, R, M,
                                  expected, tc.name, k, nruns);
            allPass = r.pass && allPass; results.push_back(r);
        }
        if (algoMask & 512) {
            TestResult r = runOne(FDF, "DenseFactorAware", F, R, M,
                                  expected, tc.name, k, nruns);
            allPass = r.pass && allPass; results.push_back(r);
        }
    }

    printTable(results);
    return allPass;
}

#endif // __LINBOX_test_frobenius_suite_H

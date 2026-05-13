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

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cassert>

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
void buildJordanMatrix(
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

    size_t d = 1;
    size_t n = d * w * h * s * (s + 1) / 2;

    M.resize(n, n);
    expectedFactors.clear();

    size_t row0 = 0;
    for (size_t lev = s; lev >= 1; --lev) {
        size_t bsz = d * lev * h;
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
    size_t nSprays = (size_t)std::ceil(std::sqrt((double)n));
    RandIter RI(F, 0, seed);

    for (size_t spray = 0; spray < nSprays; ++spray) {
        size_t i = (size_t)(rand() % n);
        size_t j;
        do { j = (size_t)(rand() % n); } while (j == i);

        Element alpha;
        do { RI.random(alpha); } while (F.isZero(alpha));

        std::vector<std::pair<size_t,Element>> rowJ;
        for (size_t col = 0; col < n; ++col) {
            Element val; F.init(val, 0);
            M.getEntry(val, j, col);
            if (!F.isZero(val)) rowJ.push_back({col, val});
        }
        for (auto &[col, val] : rowJ) {
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
        for (auto &[row, val] : colI) {
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
// Result struct
// ============================================================
struct TestResult {
    std::string algoName;
    std::string caseName;
    size_t      n;
    size_t      nnz;
    double      usertime;
    bool        pass;
};

// ============================================================
// Run one algorithm on one test case
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
    size_t k)
{
    typedef typename PolyRing::Element Polynomial;

    std::vector<Polynomial> computed;
    Givaro::Timer T;

    T.clear(); T.start();
    FO.frobeniusInvariants(computed, M, k);
    T.stop();

    // When k > 0, only verify the first k factors.
    // When k = 0, verify all expected factors.
    size_t n_check = (k > 0) ? std::min(k, expected.size()) : expected.size();

    bool pass = true;
    if (computed.size() < n_check) {
        pass = false;
    } else {
        for (size_t i = 0; i < n_check; ++i)
            if (expected[i] != computed[i]) { pass = false; break; }
    }

    if (!pass) {
        std::cout << "\n[FAIL] " << algoName << " on " << caseName << "\n";
        std::cout << "  computed (" << computed.size() << "):\n";
        for (size_t i = 0; i < computed.size(); ++i) {
            std::cout << "    [" << i << "] "; R.write(std::cout, computed[i]); std::cout << "\n";
        }
        std::cout << "  expected first " << n_check << " of " << expected.size() << ":\n";
        for (size_t i = 0; i < n_check; ++i) {
            std::cout << "    [" << i << "] "; R.write(std::cout, expected[i]); std::cout << "\n";
        }
    }

    TestResult res;
    res.algoName  = algoName;
    res.caseName  = caseName;
    res.n         = M.rowdim();
    res.nnz       = M.size();
    res.usertime  = T.usertime();
    res.pass      = pass;
    return res;
}

// ============================================================
// Print results table
// ============================================================
void printTable(const std::vector<TestResult> &results) {
    std::cout << "\n";
    std::cout << std::left
              << std::setw(20) << "Algorithm"
              << std::setw(28) << "Test Case"
              << std::setw(8)  << "n"
              << std::setw(10) << "nnz"
              << std::setw(12) << "Time (s)"
              << std::setw(8)  << "Pass"
              << "\n";
    std::cout << std::string(86, '-') << "\n";
    for (auto &r : results) {
        std::cout << std::left
                  << std::setw(20) << r.algoName
                  << std::setw(28) << r.caseName
                  << std::setw(8)  << r.n
                  << std::setw(10) << r.nnz
                  << std::setw(12) << std::fixed << std::setprecision(4) << r.usertime
                  << std::setw(8)  << (r.pass ? "PASS" : "FAIL")
                  << "\n";
    }
    std::cout << "\n";
}

// ============================================================
// Default test cases
// ============================================================
std::vector<TestParams> defaultTestCases() {
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
// ============================================================
template<
    class FrobeniusToeplitz,
    class FrobeniusToeplitzSearch,
    class FrobeniusButterfly,
    class FrobeniusButterflySearch,
    class FrobeniusDense,
    class FrobeniusDenseSearch,
    class FrobeniusLifs,
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
    const Field              &F,
    const PolyRing           &R,
    size_t k,
    int    algoMask,
    int    seed,
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

    std::vector<TestResult> results;
    bool allPass = true;

    for (auto &tc : cases) {
        size_t n = tc.w * tc.h * tc.s * (tc.s + 1) / 2;
        SparseMat M(F, n, n);
        std::vector<Polynomial> expected;
        buildJordanMatrix(M, expected, F, R, tc.s, tc.w, tc.h, tc.poly);
        sprayMatrix(M, F, seed);

        if (algoMask &  1) { auto r = runOne(FT,  "Toeplitz",        F, R, M, expected, tc.name, k); allPass &= r.pass; results.push_back(r); }
        if (algoMask &  2) { auto r = runOne(FTS, "ToeplitzSearch",  F, R, M, expected, tc.name, k); allPass &= r.pass; results.push_back(r); }
        if (algoMask &  4) { auto r = runOne(FB,  "Butterfly",       F, R, M, expected, tc.name, k); allPass &= r.pass; results.push_back(r); }
        if (algoMask &  8) { auto r = runOne(FBS, "ButterflySearch", F, R, M, expected, tc.name, k); allPass &= r.pass; results.push_back(r); }
        if (algoMask & 16) { auto r = runOne(FD,  "Dense",           F, R, M, expected, tc.name, k); allPass &= r.pass; results.push_back(r); }
        if (algoMask & 32) { auto r = runOne(FDS, "DenseSearch",     F, R, M, expected, tc.name, k); allPass &= r.pass; results.push_back(r); }
        if (algoMask & 64) { auto r = runOne(IFD, "LIFs",            F, R, M, expected, tc.name, k); allPass &= r.pass; results.push_back(r); }
    }

    printTable(results);
    return allPass;
}

#endif // __LINBOX_test_frobenius_suite_H
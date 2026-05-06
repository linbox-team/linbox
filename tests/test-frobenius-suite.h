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
    size_t s;       // number of distinct block sizes
    size_t w;       // width: repetitions per block size
    size_t h;       // height: step between successive block sizes
    int    poly;    // 0 = x, 1 = x-1, 2 = x+1
    std::string name;
};

// ============================================================
// Polynomial helpers
// ============================================================

// Evaluate p_poly at a field element val.
// poly=0: x        -> val
// poly=1: x-1      -> val - 1
// poly=2: x+1      -> val + 1
template<class Field>
typename Field::Element evalPoly(const Field &F, int poly, typename Field::Element val) {
    typedef typename Field::Element Element;
    Element res;
    switch (poly) {
        case 0: // x
            return val;
        case 1: { // x - 1
            Element one; F.init(one, 1);
            F.sub(res, val, one);
            return res;
        }
        case 2: { // x + 1
            Element one; F.init(one, 1);
            F.add(res, val, one);
            return res;
        }
        default:
            return val;
    }
}

// ============================================================
// Jordan block construction (degree 1 polynomials only for now)
// p=0: nilpotent Jordan block (eigenvalue 0), size sz
// p=1: Jordan block eigenvalue 1, size sz
// p=2: Jordan block eigenvalue -1, size sz
// Writes into M at row/col offset (row0, col0)
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

    // Diagonal entry (eigenvalue)
    if (poly != 0) {
        Element lambda;
        if (poly == 1) F.init(lambda, 1);
        else           F.init(lambda, -1); // x+1 -> eigenvalue -1
        for (size_t i = 0; i < sz; ++i)
            M.setEntry(row0 + i, col0 + i, lambda);
    }
    // Superdiagonal (subdiagonal in column-first convention — we use row-first)
    // Standard Jordan block: 1s on the subdiagonal (below diagonal)
    Element one; F.init(one, 1);
    for (size_t i = 0; i + 1 < sz; ++i)
        M.setEntry(row0 + i + 1, col0 + i, one);
}

// ============================================================
// Build Jordan form matrix from parameters (s, w, h, poly)
// Block size list: d*s*h, d*(s-1)*h, ..., d*h  each repeated w times
// For degree-1 polys d=1.
// Returns matrix and the expected invariant factor list as
// coefficient vectors (lowest degree first).
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

    // degree is 1 for p0, p1, p2
    size_t d = 1;

    // Compute total dimension
    // n = d * w * h * s*(s+1)/2
    size_t n = d * w * h * s * (s + 1) / 2;

    M.resize(n, n);

    expectedFactors.clear();

    // Block sizes: for level lev = s, s-1, ..., 1
    // block size = d * lev * h, repeated w times
    size_t row0 = 0;
    for (size_t lev = s; lev >= 1; --lev) {
        size_t bsz = d * lev * h;
        for (size_t rep = 0; rep < w; ++rep) {
            writeJordanBlock(M, F, poly, bsz, row0, row0);
            row0 += bsz;

            // Build expected invariant factor: p^bsz
            // For poly=0: x^bsz  -> coeffs [0,0,...,0,1]
            // For poly=1: (x-1)^bsz
            // For poly=2: (x+1)^bsz
            // We build the polynomial in R
            Polynomial base, factor;
            R.assign(base, R.zero);
            // set base = p_poly
            switch (poly) {
                case 0: // x
                    R.setCoeff(base, 0, (Coeff)0);
                    R.setCoeff(base, 1, (Coeff)1);
                    break;
                case 1: // x - 1
                    R.setCoeff(base, 0, (Coeff)-1);
                    R.setCoeff(base, 1, (Coeff)1);
                    break;
                case 2: // x + 1
                    R.setCoeff(base, 0, (Coeff)1);
                    R.setCoeff(base, 1, (Coeff)1);
                    break;
            }
            // factor = base^bsz
            R.assign(factor, R.one);
            for (size_t i = 0; i < bsz; ++i)
                R.mulin(factor, base);

            expectedFactors.push_back(factor);
        }
    }

    M.finalize();
}

// ============================================================
// Spray: apply sqrt(n) random elementary similarity transforms
// Each spray: add alpha*row[j] to row[i], subtract alpha*col[i] from col[j]
// This is conjugation by elementary matrix E_ij(alpha), preserving invariants.
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
        // Pick two distinct random row/col indices
        size_t i = (size_t)(rand() % n);
        size_t j;
        do { j = (size_t)(rand() % n); } while (j == i);

        Element alpha;
        do { RI.random(alpha); } while (F.isZero(alpha));

        // --- Row operation: row[i] += alpha * row[j] ---
        // Collect row j entries first to avoid aliasing
        std::vector<std::pair<size_t,Element>> rowJ;
        for (size_t col = 0; col < n; ++col) {
            Element val; F.init(val, 0);
            // SparseMatrix getEntry
            M.getEntry(val, j, col);
            if (!F.isZero(val))
                rowJ.push_back({col, val});
        }
        for (auto &[col, val] : rowJ) {
            Element cur; F.init(cur, 0);
            M.getEntry(cur, i, col);
            Element contrib; F.mul(contrib, alpha, val);
            F.addin(cur, contrib);
            M.setEntry(i, col, cur);
        }

        // --- Col operation: col[j] -= alpha * col[i] ---
        // Collect col i entries first
        std::vector<std::pair<size_t,Element>> colI;
        for (size_t row = 0; row < n; ++row) {
            Element val; F.init(val, 0);
            M.getEntry(val, row, i);
            if (!F.isZero(val))
                colI.push_back({row, val});
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
// Result struct for one (algorithm, test case) run
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
// Run one algorithm on one test case, return result
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

    bool pass = true;
    if (expected.size() > computed.size()) {
        pass = false;
    } else {
        for (size_t i = 0; i < expected.size(); ++i)
            if (expected[i] != computed[i]) { pass = false; break; }
    }
    if (!pass) {
        std::cout << "\n[FAIL] " << algoName << " on " << caseName << "\n";
        std::cout << "  computed (" << computed.size() << "):\n";
        for (size_t i = 0; i < computed.size(); ++i) {
            std::cout << "    [" << i << "] ";
            R.write(std::cout, computed[i]);
            std::cout << "\n";
        }
        std::cout << "  expected (" << expected.size() << "):\n";
        for (size_t i = 0; i < expected.size(); ++i) {
            std::cout << "    [" << i << "] ";
            R.write(std::cout, expected[i]);
            std::cout << "\n";
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
              << std::setw(24) << "Algorithm"
              << std::setw(28) << "Test Case"
              << std::setw(8)  << "n"
              << std::setw(10) << "nnz"
              << std::setw(12) << "Time (s)"
              << std::setw(8)  << "Pass"
              << "\n";
    std::cout << std::string(90, '-') << "\n";
    for (auto &r : results) {
        std::cout << std::left
                  << std::setw(24) << r.algoName
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
// Default hardcoded test cases (tall, flat, triangular, combo)
// ============================================================
std::vector<TestParams> defaultTestCases() {
    return {
        // Tall: p=x, s=3, w=1, h=10
        {3, 1, 10, 0, "tall_x"},
        // Flat: p=x-1, s=3, w=5, h=1
        {3, 5,  1, 1, "flat_xm1"},
        // Triangular: p=x, s=4, w=4, h=4
        {4, 4,  4, 0, "tri_x"},
        // Tall with x+1
        {3, 1, 10, 2, "tall_xp1"},
        // Flat with x
        {3, 5,  1, 0, "flat_x"},
    };
}

// ============================================================
// Main test suite runner
// Algorithms bitmask: bit0=Toeplitz, bit1=Butterfly, bit2=Dense, bit3=Search
// ============================================================
template<
    class FrobeniusToeplitz,
    class FrobeniusButterfly,
    class FrobeniusDense,
    class FrobeniusSearch,
    class Field,
    class PolyRing>
bool runSuite(
    FrobeniusToeplitz  &FT,
    FrobeniusButterfly &FB,
    FrobeniusDense     &FD,
    FrobeniusSearch    &FS,
    const Field        &F,
    const PolyRing     &R,
    size_t k,
    int    algoMask,   // bitmask: bit0=Toeplitz, bit1=Butterfly, bit2=Dense, bit3=Search
    int    seed,
    // Optional custom params (0 means use defaults)
    size_t custom_s = 0,
    size_t custom_w = 0,
    size_t custom_h = 0,
    int    custom_poly = -1)
{
    typedef typename PolyRing::Element Polynomial;
    typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;

    // Build test case list
    std::vector<TestParams> cases = defaultTestCases();

    // If custom params provided, prepend a custom case
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
        // Build Jordan matrix and expected invariant factors
        size_t n = tc.w * tc.h * tc.s * (tc.s + 1) / 2; // d=1 for all current polys
        SparseMat M(F, n, n);
        std::vector<Polynomial> expected;
        buildJordanMatrix(M, expected, F, R, tc.s, tc.w, tc.h, tc.poly);

        // Apply similarity sprays
        sprayMatrix(M, F, seed);

        // Run selected algorithms
        if (algoMask & 1) {
            auto res = runOne(FT, "Toeplitz", F, R, M, expected, tc.name, k);
            allPass &= res.pass;
            results.push_back(res);
        }
        if (algoMask & 2) {
            auto res = runOne(FB, "Butterfly", F, R, M, expected, tc.name, k);
            allPass &= res.pass;
            results.push_back(res);
        }
        if (algoMask & 4) {
            auto res = runOne(FD, "Dense", F, R, M, expected, tc.name, k);
            allPass &= res.pass;
            results.push_back(res);
        }
        if (algoMask & 8) {
            auto res = runOne(FS, "Search", F, R, M, expected, tc.name, k);
            allPass &= res.pass;
            results.push_back(res);
        }
    }

    printTable(results);
    return allPass;
}

#endif // __LINBOX_test_frobenius_suite_H
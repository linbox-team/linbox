/* tests/test-frobenius-egns-bench.C
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 * Written by Omesh Dhar Dwivedi <odd23@drexel.edu>
 *
 * Extended accuracy + timing benchmark for FrobeniusLargeRank (EGNS).
 *
 * Two test groups:
 *   1. FLAT SWEEP  — systematic grid of flat Jordan shapes (many small cases)
 *   2. LARGE CASES — a handful of larger matrices to expose timing differences
 *
 * Flat cases: blocks (k,k,...) i.e. all blocks of equal size, and also
 *   mixed shapes.  Matrix sizes range from n=4 to n≈60.
 * Large cases: n = 100, 200, 400 (single-eigenvalue, flat structure).
 *
 * Output: a table matching test-frobenius-suite.h format, with timing.
 *
 * No dependency on test-frobenius-suite.h.
 * 
 * 
 * Include order: algorithm headers before NTL ring wrapper.
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>

#include "linbox/algorithms/frobenius-large-generic.h"
#include "linbox/algorithms/frobenius-large-rank.h"

#include "linbox/ring/ntl/ntl-lzz_px.h"
#include "linbox/matrix/sparse-matrix.h"
#include "givaro/givtimer.h"

using namespace LinBox;

typedef NTL_zz_pX                                    PolynomialRing;
typedef PolynomialRing::CoeffField                   Field;
typedef PolynomialRing::Coeff                        Coeff;
typedef PolynomialRing::Element                      Polynomial;
typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;

// ---- Build (x - eigenval)^exp ----
static Polynomial makeXminusLambdaPow(
    const PolynomialRing &R,
    const Field          &F,
    int                   eigenval,
    size_t                exp)
{
    Coeff c_one, c_neg;
    F.assign(c_one, F.one);
    F.init(c_neg, -eigenval);

    Polynomial base;
    R.assign(base, R.zero);
    R.setCoeff(base, 0, c_neg);
    R.setCoeff(base, 1, c_one);

    Polynomial result;
    R.assign(result, R.one);
    for (size_t i = 0; i < exp; i++) R.mulin(result, base);
    return result;
}

// ---- Build flat Jordan matrix from arbitrary block sizes ----
//   Sorts blockSizes descending, lays out on diagonal.
//   expected[i] = (x - eigenval)^blockSizes[i]  (after sort)
static void buildFlatMatrix(
    SparseMat               &A,
    std::vector<Polynomial> &expected,
    const Field             &F,
    const PolynomialRing    &R,
    std::vector<size_t>      blockSizes,   // by value, sorted here
    int                      eigenval)
{
    std::sort(blockSizes.begin(), blockSizes.end(), std::greater<size_t>());

    Field::Element lam, one;
    F.init(lam, eigenval);
    F.assign(one, F.one);

    size_t row = 0;
    for (size_t bsz : blockSizes) {
        for (size_t j = 0; j < bsz; j++) {
            if (!F.isZero(lam))   A.setEntry(row+j, row+j,   lam);
            if (j + 1 < bsz)     A.setEntry(row+j, row+j+1, one);
        }
        row += bsz;
    }
    A.finalize();

    expected.clear();
    for (size_t bsz : blockSizes)
        expected.push_back(makeXminusLambdaPow(R, F, eigenval, bsz));
}

// ---- Strip trailing 1s ----
static std::vector<Polynomial>
stripOnes(const std::vector<Polynomial> &fs, const PolynomialRing &R)
{
    std::vector<Polynomial> out = fs;
    while (!out.empty() && R.isOne(out.back()))
        out.pop_back();
    return out;
}

// ---- Compare invariant lists ----
static bool sameInvariants(
    const std::vector<Polynomial> &a,
    const std::vector<Polynomial> &b,
    const PolynomialRing           &R)
{
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++)
        if (!R.areEqual(a[i], b[i])) return false;
    return true;
}

// ---- Result record ----
struct BenchResult {
    std::string algoName;
    std::string caseName;
    size_t      n;
    double      time_s;
    bool        pass;
};

// ---- Run one algorithm, time it, return result ----
template<class Alg>
static BenchResult runOne(
    Alg                            &alg,
    const std::string              &algName,
    const std::string              &caseName,
    const SparseMat                &A,
    const std::vector<Polynomial>  &expected,
    const PolynomialRing           &R)
{
    std::vector<Polynomial> fs;
    Givaro::Timer T;
    T.clear(); T.start();
    alg.frobeniusInvariants(fs, A, 0);
    T.stop();

    std::vector<Polynomial> got = stripOnes(fs, R);
    bool pass = sameInvariants(got, expected, R);

    if (!pass) {
        std::cout << "  [FAIL] " << algName << " on " << caseName << "\n";
        std::cout << "    expected (" << expected.size() << "): ";
        for (auto &f : expected) { R.write(std::cout, f); std::cout << " "; }
        std::cout << "\n";
        std::cout << "    got      (" << got.size() << "): ";
        for (auto &f : got)      { R.write(std::cout, f); std::cout << " "; }
        std::cout << "\n";
    }

    return { algName, caseName, A.rowdim(), T.usertime(), pass };
}

// ---- Print results table ----
static void printTable(const std::vector<BenchResult> &results)
{
    std::cout << "\n";
    std::cout << std::left
              << std::setw(16) << "Algorithm"
              << std::setw(36) << "Case"
              << std::setw(8)  << "n"
              << std::setw(12) << "Time (s)"
              << std::setw(8)  << "Pass"
              << "\n";
    std::cout << std::string(80, '-') << "\n";
    for (auto &r : results) {
        std::cout << std::left
                  << std::setw(16) << r.algoName
                  << std::setw(36) << r.caseName
                  << std::setw(8)  << r.n
                  << std::setw(12) << std::fixed << std::setprecision(4) << r.time_s
                  << std::setw(8)  << (r.pass ? "PASS" : "FAIL")
                  << "\n";
    }
    std::cout << "\n";
}

int main()
{
    int prime = 65521;
    Field F(prime);
    PolynomialRing R(F);

    FrobeniusLarge<PolynomialRing>          FT(R);
    FrobeniusLargeButterfly<PolynomialRing> FB(R);
    FrobeniusLargeRank<PolynomialRing>      FR(R);

    std::vector<BenchResult> results;
    int total = 0, passed = 0;

    // Helper lambda: run all three algorithms on one case
    auto bench = [&](const std::string    &desc,
                     SparseMat            &A,
                     std::vector<Polynomial> &expected)
    {
        for (auto [alg_ptr, name] : std::initializer_list<std::pair<void*, const char*>>{
                {(void*)&FT, "Toeplitz"},
                {(void*)&FB, "Butterfly"},
                {(void*)&FR, "EGNS"}})
        {
            BenchResult r;
            if      (name == std::string("Toeplitz"))  r = runOne(FT, name, desc, A, expected, R);
            else if (name == std::string("Butterfly")) r = runOne(FB, name, desc, A, expected, R);
            else                                       r = runOne(FR, name, desc, A, expected, R);
            results.push_back(r);
            total++;
            passed += r.pass;
        }
    };

    // ================================================================
    // GROUP 1: FLAT SWEEP  (systematic, many cases, small-to-medium n)
    //
    // Shape family: b copies of a single block size s  → n = b*s
    // Eigenvalues: 0 and 1
    // ================================================================
    std::cout << "=== GROUP 1: Flat sweep (uniform block sizes) ===\n\n";

    // (block_size, num_blocks) pairs — n = block_size * num_blocks
    std::vector<std::pair<size_t,size_t>> uniformShapes = {
        {1, 4},  // n=4,  e=1, trivial
        {2, 2},  // n=4,  e=2
        {1, 6},  // n=6
        {2, 3},  // n=6,  e=2
        {3, 2},  // n=6,  e=3
        {1, 8},  // n=8
        {2, 4},  // n=8,  e=2
        {4, 2},  // n=8,  e=4
        {1,10},  // n=10
        {2, 5},  // n=10, e=2
        {5, 2},  // n=10, e=5
        {3, 4},  // n=12, e=3
        {4, 3},  // n=12, e=4
        {6, 2},  // n=12, e=6
        {4, 4},  // n=16, e=4
        {5, 4},  // n=20, e=5
        {6, 4},  // n=24, e=6
        {5, 5},  // n=25, e=5
        {6, 5},  // n=30, e=6
        {8, 4},  // n=32, e=8
        {7, 5},  // n=35, e=7
        {8, 5},  // n=40, e=8
        {9, 5},  // n=45, e=9
       {10, 5},  // n=50, e=10
       {12, 5},  // n=60, e=12
    };

    for (auto [bsz, nb] : uniformShapes) {
        std::vector<size_t> blocks(nb, bsz);
        size_t n = bsz * nb;

        for (int ev : {0, 1}) {
            std::string desc = "unif_b" + std::to_string(bsz)
                             + "x" + std::to_string(nb)
                             + "_ev" + std::to_string(ev)
                             + "_n"  + std::to_string(n);
            SparseMat A(F, n, n);
            std::vector<Polynomial> expected;
            buildFlatMatrix(A, expected, F, R, blocks, ev);
            bench(desc, A, expected);
        }
    }

    // ================================================================
    // Mixed-block flat sweep:  blocks (s, s, s-1, s-1, ..., 1, 1)
    // "staircase" — tests the triangular-solve path more thoroughly
    // ================================================================
    std::cout << "=== GROUP 1b: Flat sweep (staircase block sizes) ===\n\n";

    for (size_t s = 2; s <= 6; s++) {
        for (int ev : {0, 2}) {
            std::vector<size_t> blocks;
            for (size_t k = s; k >= 1; k--) { blocks.push_back(k); blocks.push_back(k); }
            size_t n = std::accumulate(blocks.begin(), blocks.end(), (size_t)0);

            std::string desc = "stair_s" + std::to_string(s)
                             + "_ev" + std::to_string(ev)
                             + "_n"  + std::to_string(n);
            SparseMat A(F, n, n);
            std::vector<Polynomial> expected;
            buildFlatMatrix(A, expected, F, R, blocks, ev);
            bench(desc, A, expected);
        }
    }

    // ================================================================
    // GROUP 2: LARGE CASES  (fewer, bigger matrices, timing focus)
    //
    // Shape: (b copies of block size s), n = b*s
    // Chosen so EGNS stays feasible (e = s small) but n is large.
    // EGNS cost ~ O(eta * e^2 * n); Toeplitz/Butterfly cost ~ O(n^2 deg).
    // Expectation: EGNS faster when e is small and n is large.
    // ================================================================
    std::cout << "=== GROUP 2: Large matrices (timing) ===\n\n";

    // {block_size, num_blocks, eigenval}
    std::vector<std::tuple<size_t,size_t,int>> largeCases = {
        {2,  50, 0},   // n=100,  e=2  (many small blocks)
        {2, 100, 1},   // n=200,  e=2
        {2, 200, 0},   // n=400,  e=2
        {3,  33, 0},   // n=99,   e=3
        {3,  67, 1},   // n=201,  e=3
        {4,  25, 0},   // n=100,  e=4
        {4,  50, 2},   // n=200,  e=4
        {5,  20, 0},   // n=100,  e=5
        {5,  40, 1},   // n=200,  e=5
        {6,  17, 0},   // n=102,  e=6
        {8,  25, 0},   // n=200,  e=8  (taller blocks — harder for EGNS)
        {10, 20, 1},   // n=200,  e=10
        {2, 300, 0},   // n=600,  e=2  (stress test — many small blocks)
    };

    for (auto [bsz, nb, ev] : largeCases) {
        std::vector<size_t> blocks(nb, bsz);
        size_t n = bsz * nb;

        std::string desc = "large_b" + std::to_string(bsz)
                         + "x"  + std::to_string(nb)
                         + "_ev" + std::to_string(ev)
                         + "_n"  + std::to_string(n);
        SparseMat A(F, n, n);
        std::vector<Polynomial> expected;
        buildFlatMatrix(A, expected, F, R, blocks, ev);
        bench(desc, A, expected);
    }

    // ================================================================
    // Summary table
    // ================================================================
    printTable(results);

    std::cout << "==============================\n";
    std::cout << "Result: " << passed << "/" << total << " passed\n";
    return (passed == total) ? 0 : 1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\\:0,t0,+0,=s
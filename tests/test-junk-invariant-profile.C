/* linbox/tests/test-junk-invariant-profile.C
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 *
 * One-shot structural diagnostic for the junk invariant factors created by
 * a rank-(k-1) additive preconditioner.
 *
 * Convention:
 *
 *     f_1 = minpoly(A),   f_{i+1} | f_i,
 *
 * and the preconditioner used for f_k has requested rank r=k-1.
 * In standard Smith order s_1 | ... | s_n, the Villard-shape reconstruction is
 *
 *     sigma_i = t_i                         (1 <= i <= r),
 *     sigma_i = s_{i-r} t_i                 (r < i <= n).
 *
 * The test chooses one k (randomly unless -k is supplied), computes the full
 * invariant-factor list of A+P for each selected preconditioner family, uses the
 * exact known invariant factors of the constructed A, reconstructs every t_i,
 * prints the three indexed lists, and reports:
 *
 *   - exact-rank check for P;
 *   - exact shape divisions;
 *   - t_1 | ... | t_n;
 *   - coprimality of every t_i with charpoly(A);
 *   - lower-slot junk versus concentration in t_n;
 *   - the junk degree budget;
 *   - squarefreeness of T = product_i t_i;
 *   - the degree identity deg sigma_n = sum_{j<=k} deg f_j.
 *
 * The input A is a constructed Jordan or mixed-primary companion-block matrix,
 * so its invariant factors are taken from the construction rather than from a
 * second probabilistic FNF computation.
 *
 * Typical commands:
 *
 *   make test-junk-invariant-profile
 *   ./tests/test-junk-invariant-profile -p 101 -s 4 -w 2 -H 5 -a 012
 *   ./tests/test-junk-invariant-profile -p 5 -s 6 -w 2 -H 5 -k 4 -a 2
 *
 * Options:
 *   -p P       base prime, default 10000019
 *   -e E       extension degree, default 1
 *   -r R       random seed, default time(NULL)
 *   -k K       selected kth largest invariant; 0=random nontrivial k >= 2
 *   -a A       0=Toeplitz, 1=Butterfly, 2=Dense; default 012
 *   -s S       custom Jordan: number of distinct block sizes
 *   -w W       custom Jordan: repetitions per block size
 *   -H H       custom Jordan: height/step between block sizes
 *   -q Q       custom Jordan polynomial: 0=x, 1=x-1, 2=x+1
 *   -u U       random similarity spray: 0=no, 1=yes; default 1
 *   -c C       built-in Jordan case index when no custom case; default 0
 *   -m M       mixed case: 0=off, 1..5 selects a built-in mixed case
 *   -z Z       maximum full-FNF attempts on the same A+P; default 3
 */

#ifndef DISABLE_COMMENTATOR
#define DISABLE_COMMENTATOR
#endif

#include "linbox/linbox-config.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "linbox/ring/ntl.h"
#include "linbox/util/commentator.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/vector/blas-vector.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/sum.h"
#include "linbox/blackbox/toeplitz.h"
#include "linbox/blackbox/butterfly.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/algorithms/frobenius-large-generic.h"

#include "test-frobenius-suite.h"

using namespace LinBox;

typedef NTL_zz_p  Field;
typedef NTL_zz_pX PolyRing;
typedef PolyRing::Element Polynomial;
typedef PolyRing::Coeff   Coeff;
typedef Field::Element    Element;
typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;
typedef MatrixDomain<Field> MatrixDom;
typedef MatrixDom::OwnMatrix DenseMat;

static int parseAlgoMask(const char *s)
{
    if (s == nullptr || std::strlen(s) == 0) return 7; // all three
    int mask = 0;
    for (size_t i = 0; i < std::strlen(s); ++i) {
        if (s[i] == '0') mask |= 1;
        if (s[i] == '1') mask |= 2;
        if (s[i] == '2') mask |= 4;
    }
    return mask ? mask : 7;
}

static const char *algoName(int bit)
{
    switch (bit) {
        case 1: return "Toeplitz(A+T1T2)";
        case 2: return "Butterfly(A+VtDU)";
        case 4: return "Dense(A+UV)";
        default: return "Unknown";
    }
}

static long degPoly(PolyRing &R, const Polynomial &f)
{
    if (R.isZero(f)) return -1;
    return (long) R.deg(f);
}

// In the Jordan-suite expected list, entries are largest-first and omit trailing 1s.
// For k past the nontrivial range, f_k is 1.
static Polynomial expectedAt(PolyRing &R,
                             const std::vector<Polynomial> &expected,
                             size_t k)
{
    Polynomial one;
    R.assign(one, R.one);
    if (k == 0) return one;
    if (k <= expected.size()) return expected[k - 1];
    return one;
}

static size_t prefixDegree(PolyRing &R,
                           const std::vector<Polynomial> &expected,
                           size_t k)
{
    size_t ans = 0;
    size_t upto = std::min(k, expected.size());
    for (size_t i = 0; i < upto; ++i) ans += (size_t) R.deg(expected[i]);
    return ans;
}

// ============================================================
// Mixed-primary stress cases
// ============================================================
// The Jordan-suite builder only gives one primary at a time: powers of
// x, x-1, or x+1.  For mixed-primary tests we instead build rational
// canonical blocks directly.  If the expected factors are supplied in
// largest-first divisibility order
//
//     f_1, f_2, ..., f_m,     f_{i+1} | f_i,
//
// then diag(C(f_1),...,C(f_m)) has exactly those nontrivial invariant
// factors.  This lets us test examples such as
//
//     x^30 (x-1)^5,
//     x^20 (x-1)^5,
//     x^10 (x-1)^3,
//     (x-1)^3,
//
// where the x-primary and (x-1)-primary exponent profiles are both
// nonincreasing down the largest-first list.

struct PrimaryPowers {
    size_t x;
    size_t xm1;
    size_t xp1;
};

struct MixedCase {
    std::string name;
    std::vector<PrimaryPowers> profile; // largest-first exponent triples
};

static Polynomial makeLinear(PolyRing &R, int poly)
{
    Polynomial base;
    R.assign(base, R.zero);

    // poly=0: x, poly=1: x-1, poly=2: x+1.
    switch (poly) {
        case 0: R.setCoeff(base, 0, (Coeff)0);  R.setCoeff(base, 1, (Coeff)1); break;
        case 1: R.setCoeff(base, 0, (Coeff)-1); R.setCoeff(base, 1, (Coeff)1); break;
        case 2: R.setCoeff(base, 0, (Coeff)1);  R.setCoeff(base, 1, (Coeff)1); break;
        default: R.setCoeff(base, 0, (Coeff)0); R.setCoeff(base, 1, (Coeff)1); break;
    }
    return base;
}

static Polynomial powPoly(PolyRing &R, const Polynomial &base, size_t e)
{
    Polynomial out;
    R.assign(out, R.one);
    for (size_t i = 0; i < e; ++i) R.mulin(out, base);
    return out;
}

static Polynomial makeMixedFactor(PolyRing &R, const PrimaryPowers &pp)
{
    Polynomial f, lin, part;
    R.assign(f, R.one);

    lin = makeLinear(R, 0);
    part = powPoly(R, lin, pp.x);
    R.mulin(f, part);

    lin = makeLinear(R, 1);
    part = powPoly(R, lin, pp.xm1);
    R.mulin(f, part);

    lin = makeLinear(R, 2);
    part = powPoly(R, lin, pp.xp1);
    R.mulin(f, part);

    return f;
}

static bool dividesPoly(PolyRing &R, const Polynomial &small, const Polynomial &large)
{
    Polynomial g;
    R.assign(g, small);
    R.gcdin(g, large);
    return R.areEqual(g, small);
}

static void assertLargestFirstChain(PolyRing &R, const std::vector<Polynomial> &factors)
{
    for (size_t i = 1; i < factors.size(); ++i) {
        if (!dividesPoly(R, factors[i], factors[i-1])) {
            std::cerr << "Internal mixed-primary case is not in largest-first divisibility order at i="
                      << i << "\n  f_i=";
            R.write(std::cerr, factors[i-1]);
            std::cerr << "\n  f_{i+1}=";
            R.write(std::cerr, factors[i]);
            std::cerr << "\n";
            std::abort();
        }
    }
}

static void writeCompanionBlock(SparseMat &M,
                                const Field &F,
                                PolyRing &R,
                                const Polynomial &f,
                                size_t row0)
{
    const long dLong = R.deg(f);
    if (dLong <= 0) return;
    const size_t d = (size_t) dLong;

    Element one;
    F.assign(one, F.one);

    // Standard companion orientation:
    // rows 1..d-1 have a 1 in the previous column, and the last column is
    // -a_0,...,-a_{d-1} for f=x^d+a_{d-1}x^{d-1}+...+a_0.
    for (size_t i = 1; i < d; ++i)
        M.setEntry(row0 + i, row0 + i - 1, one);

    for (size_t j = 0; j < d; ++j) {
        Coeff c;
        R.getCoeff(c, f, j);
        F.negin(c);
        if (!F.isZero(c))
            M.setEntry(row0 + j, row0 + d - 1, c);
    }
}

static void buildCompanionInvariantMatrix(SparseMat &A,
                                          std::vector<Polynomial> &expected,
                                          const Field &F,
                                          PolyRing &R,
                                          const std::vector<Polynomial> &factors)
{
    expected = factors;
    assertLargestFirstChain(R, expected);

    size_t n = 0;
    for (const auto &f : expected) n += (size_t) R.deg(f);

    A.resize(n, n);
    size_t row0 = 0;
    for (const auto &f : expected) {
        writeCompanionBlock(A, F, R, f, row0);
        row0 += (size_t) R.deg(f);
    }
    A.finalize();
}

static std::vector<MixedCase> defaultMixedCases()
{
    std::vector<MixedCase> out;

    // Prompt-style mixed primary.  The literal list
    //     x^30(x-1)^5, x^20(x-1)^5, x^10, (x-1)^3
    // is not itself an invariant-factor chain, since (x-1)^3 does not
    // divide x^10.  The coprime x^10 and (x-1)^3 pieces fold into
    // x^10(x-1)^3 as an invariant factor.
    out.push_back({"mixed_prompt_folded", {
        {30, 5, 0},
        {20, 5, 0},
        {10, 3, 0}
    }});

    // Same prompt-style profile, but with an additional pure (x-1)^3 tail.
    // Exponent profiles:
    //   x:   30,20,10,0
    //   x-1:  5, 5, 3,3
    out.push_back({"mixed_prompt_tail", {
        {30, 5, 0},
        {20, 5, 0},
        {10, 3, 0},
        { 0, 3, 0}
    }});

    // Three-primary staircase using x, x-1, and x+1.
    out.push_back({"mixed_3primary", {
        {16, 8, 4},
        {12, 5, 4},
        { 7, 5, 2},
        { 3, 1, 0}
    }});

    // Flat mixed case: many equal invariant factors, each with two primaries.
    out.push_back({"mixed_flat", {
        {3, 2, 0}, {3, 2, 0}, {3, 2, 0}, {3, 2, 0},
        {3, 2, 0}, {3, 2, 0}, {3, 2, 0}, {3, 2, 0}
    }});

    // Longer mixed staircase to stress repeated degree plateaus.
    out.push_back({"mixed_stair", {
        {8, 5, 0},
        {8, 5, 0},
        {6, 4, 0},
        {6, 4, 0},
        {4, 3, 0},
        {4, 2, 0},
        {2, 2, 0},
        {1, 0, 0}
    }});

    return out;
}


// ============================================================
// Polynomial arithmetic helpers for squarefree/junk-product tests
// ============================================================

static Polynomial productFactors(PolyRing &R,
                                 const std::vector<Polynomial> &v,
                                 size_t begin,
                                 size_t end)
{
    Polynomial out;
    R.assign(out, R.one);
    end = std::min(end, v.size());
    for (size_t i = begin; i < end; ++i)
        R.mulin(out, v[i]);
    return out;
}

static Polynomial productFactors(PolyRing &R,
                                 const std::vector<Polynomial> &v)
{
    return productFactors(R, v, 0, v.size());
}

static std::vector<Polynomial> stripTrailingOnes(PolyRing &R,
                                                 const std::vector<Polynomial> &in)
{
    std::vector<Polynomial> out = in;
    while (!out.empty() && R.isOne(out.back()))
        out.pop_back();
    return out;
}

// Convert a largest-first nontrivial list into a standard Smith-order list
// of length n: s_1 | ... | s_n, including leading 1s.
static std::vector<Polynomial> toStandardLength(PolyRing &R,
                                                const std::vector<Polynomial> &largestFirstRaw,
                                                size_t n)
{
    std::vector<Polynomial> largestFirst = stripTrailingOnes(R, largestFirstRaw);
    std::vector<Polynomial> stdList;
    Polynomial one;
    R.assign(one, R.one);

    const size_t nontriv = std::min(largestFirst.size(), n);
    for (size_t i = 0; i < n - nontriv; ++i)
        stdList.push_back(one);

    for (size_t j = 0; j < nontriv; ++j)
        stdList.push_back(largestFirst[nontriv - 1 - j]);

    return stdList;
}

// Polynomial derivative using only the LinBox polynomial-ring interface.
static void derivativePoly(Polynomial &df,
                           PolyRing &R,
                           const Field &F,
                           const Polynomial &f)
{
    R.assign(df, R.zero);
    const long d = degPoly(R, f);
    if (d <= 0) return;

    for (long i = 1; i <= d; ++i) {
        Coeff c, scalar, prod;
        R.getCoeff(c, f, (size_t)i);
        if (F.isZero(c)) continue;
        F.init(scalar, (int)i);
        F.mul(prod, c, scalar);
        if (!F.isZero(prod))
            R.setCoeff(df, (size_t)(i - 1), prod);
    }
}

// Exact quotient by a monic denominator, implemented by long division.
// All invariant factors here are monic, so no leading-coefficient inverse is needed.
static bool exactQuotientMonic(Polynomial &q,
                               PolyRing &R,
                               const Field &F,
                               const Polynomial &num,
                               const Polynomial &den)
{
    R.assign(q, R.zero);
    if (R.isZero(den)) return false;
    if (R.isOne(den)) {
        R.assign(q, num);
        return true;
    }
    if (R.isZero(num)) return true;

    const long denDeg = degPoly(R, den);
    if (denDeg < 0) return false;

    Polynomial rem;
    R.assign(rem, num);

    while (!R.isZero(rem) && degPoly(R, rem) >= denDeg) {
        const long remDeg = degPoly(R, rem);
        const size_t shift = (size_t)(remDeg - denDeg);

        Coeff lead;
        R.getCoeff(lead, rem, (size_t)remDeg);
        if (F.isZero(lead)) break;

        R.setCoeff(q, shift, lead);

        for (long j = 0; j <= denDeg; ++j) {
            Coeff dj, prod, cur;
            R.getCoeff(dj, den, (size_t)j);
            if (F.isZero(dj)) continue;
            F.mul(prod, lead, dj);
            R.getCoeff(cur, rem, shift + (size_t)j);
            F.subin(cur, prod);
            R.setCoeff(rem, shift + (size_t)j, cur);
        }
    }

    return R.isZero(rem);
}

static bool isSquarefreePoly(PolyRing &R,
                             const Field &F,
                             const Polynomial &f)
{
    if (R.isZero(f)) return false;
    if (R.isOne(f)) return true;

    Polynomial df, g;
    derivativePoly(df, R, F, f);
    R.assign(g, f);
    R.gcdin(g, df);
    return R.isOne(g);
}

static bool dividesExact(Polynomial &quot,
                         PolyRing &R,
                         const Field &F,
                         const Polynomial &num,
                         const Polynomial &den)
{
    return exactQuotientMonic(quot, R, F, num, den);
}


struct ProbeResult {
    std::vector<Polynomial> sigmaLargest; // largest-first, nontrivial entries
    size_t actualRank = 0;
    size_t fnfAttemptsUsed = 0;
    bool fullDegreeOK = false;
};

// Materialize a black box by applying it to the standard basis, then perform
// exact Gaussian elimination over F.  This is intentionally dense and cubic:
// the test runs only one k and is meant as a structural diagnostic.
template<class Blackbox>
static size_t exactBlackboxRank(const Blackbox &P, const Field &F)
{
    const size_t rows = P.rowdim();
    const size_t cols = P.coldim();

    std::vector<Element> M(rows * cols);
    for (size_t i = 0; i < rows * cols; ++i)
        F.assign(M[i], F.zero);

    for (size_t j = 0; j < cols; ++j) {
        BlasVector<Field> x(F, cols), y(F, rows);
        for (size_t c = 0; c < cols; ++c) F.assign(x[c], F.zero);
        F.assign(x[j], F.one);
        P.apply(y, x);
        for (size_t i = 0; i < rows; ++i)
            F.assign(M[i * cols + j], y[i]);
    }

    size_t rank = 0;
    for (size_t col = 0; col < cols && rank < rows; ++col) {
        size_t pivot = rank;
        while (pivot < rows && F.isZero(M[pivot * cols + col])) ++pivot;
        if (pivot == rows) continue;

        if (pivot != rank) {
            for (size_t j = col; j < cols; ++j)
                std::swap(M[pivot * cols + j], M[rank * cols + j]);
        }

        Element invPivot;
        F.inv(invPivot, M[rank * cols + col]);

        for (size_t i = rank + 1; i < rows; ++i) {
            if (F.isZero(M[i * cols + col])) continue;

            Element factor;
            F.mul(factor, M[i * cols + col], invPivot);
            for (size_t j = col; j < cols; ++j) {
                Element prod;
                F.mul(prod, factor, M[rank * cols + j]);
                F.subin(M[i * cols + j], prod);
            }
        }
        ++rank;
    }
    return rank;
}

template<class Blackbox>
static void computeFullInvariants(ProbeResult &res,
                                  PolyRing &R,
                                  const Blackbox &Ap,
                                  size_t n,
                                  size_t maxAttempts)
{
    res.sigmaLargest.clear();
    res.fullDegreeOK = false;
    res.fnfAttemptsUsed = 0;

    for (size_t attempt = 1; attempt <= std::max<size_t>(1, maxAttempts); ++attempt) {
        std::vector<Polynomial> candidate;
        FrobeniusLarge<PolyRing> solver(R);
        solver.frobeniusInvariants(candidate, Ap, 0);
        candidate = stripTrailingOnes(R, candidate);

        Polynomial charAp = productFactors(R, candidate);
        res.sigmaLargest = candidate;
        res.fnfAttemptsUsed = attempt;
        res.fullDegreeOK = (degPoly(R, charAp) == (long)n);
        if (res.fullDegreeOK) break;
    }
}

template<class Blackbox>
static void probeToeplitz(ProbeResult &res,
                          FrobeniusLarge<PolyRing> &helper,
                          const Field &F,
                          PolyRing &R,
                          const Blackbox &A,
                          size_t k,
                          size_t fnfAttempts)
{
    if (k <= 1) {
        res.actualRank = 0;
        computeFullInvariants(res, R, A, A.rowdim(), fnfAttempts);
        return;
    }

    const size_t n = A.rowdim();
    const size_t rank = k - 1;

    Polynomial u, v;
    helper.randomPolynomial(u, n + k - 3);
    helper.randomPolynomial(v, n + k - 3);

    typedef Toeplitz<Field, PolyRing> Toep;
    Toep T1(R, u, n, rank);
    Toep T2(R, v, rank, n);

    typedef Compose<Toep, Toep> Prec;
    Prec P(T1, T2);
    res.actualRank = exactBlackboxRank(P, F);

    typedef Sum<Blackbox, Prec> Preconditioned;
    Preconditioned Ap(A, P);
    computeFullInvariants(res, R, Ap, n, fnfAttempts);
}

template<class Blackbox>
static void probeButterfly(ProbeResult &res,
                           FrobeniusLarge<PolyRing> &helper,
                           const Field &F,
                           PolyRing &R,
                           const Blackbox &A,
                           size_t k,
                           size_t fnfAttempts)
{
    (void)helper;
    if (k <= 1) {
        res.actualRank = 0;
        computeFullInvariants(res, R, A, A.rowdim(), fnfAttempts);
        return;
    }

    const size_t n = A.rowdim();
    const size_t rank = k - 1;

    typedef CekstvSwitch<Field> Switch;
    typedef Butterfly<Field, Switch> BB;
    typedef TransposeOwner<BB> BBt;

    Field::RandIter RI(F);
    typename Switch::Factory facU(RI);
    typename Switch::Factory facV(RI);

    BB U(F, n, facU);
    BB V(F, n, facV);
    BBt Vt(V);

    SparseMatrix<Field> D(F, n, n);
    for (size_t i = 0; i < rank; ++i) {
        Element e;
        do { RI.random(e); } while (F.isZero(e));
        D.setEntry(n - i - 1, n - i - 1, e);
    }

    typedef Compose<SparseMatrix<Field>, BB> DU_t;
    DU_t DU(D, U);
    typedef Compose<BBt, DU_t> Prec;
    Prec P(Vt, DU);
    res.actualRank = exactBlackboxRank(P, F);

    typedef Sum<Blackbox, Prec> Preconditioned;
    Preconditioned Ap(A, P);
    computeFullInvariants(res, R, Ap, n, fnfAttempts);
}

template<class Blackbox>
static void probeDense(ProbeResult &res,
                       FrobeniusLarge<PolyRing> &helper,
                       const Field &F,
                       PolyRing &R,
                       const Blackbox &A,
                       size_t k,
                       size_t fnfAttempts)
{
    (void)helper;
    if (k <= 1) {
        res.actualRank = 0;
        computeFullInvariants(res, R, A, A.rowdim(), fnfAttempts);
        return;
    }

    const size_t n = A.rowdim();
    const size_t rank = k - 1;

    Field::RandIter RI(F);
    DenseMat U(F, n, rank);
    DenseMat V(F, rank, n);

    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < rank; ++j) {
            Element e; RI.random(e); U.setEntry(i, j, e);
        }

    for (size_t i = 0; i < rank; ++i)
        for (size_t j = 0; j < n; ++j) {
            Element e; RI.random(e); V.setEntry(i, j, e);
        }

    typedef Compose<DenseMat, DenseMat> Prec;
    Prec P(U, V);
    res.actualRank = exactBlackboxRank(P, F);

    typedef Sum<Prec, Blackbox> Preconditioned;
    Preconditioned Ap(P, A);
    computeFullInvariants(res, R, Ap, n, fnfAttempts);
}

struct JunkProfile {
    std::vector<Polynomial> sStd;
    std::vector<Polynomial> sigmaStd;
    std::vector<Polynomial> tStd;
    std::vector<bool> tDefined;

    Polynomial charA;
    Polynomial charAp;
    Polynomial preservedOldProduct;
    Polynomial junkProduct;
    Polynomial junkQuotient;
    Polynomial repeatedFactorGCD;

    bool fullDegreeOK = false;
    bool rankOK = false;
    bool shapeOK = false;
    bool charQuotientOK = false;
    bool productMatchesQuotient = false;
    bool coprimeOK = false;
    bool chainOK = false;
    bool squarefreeDefined = false;
    bool squarefree = false;
    bool budgetOK = false;
    bool allLowerOne = false;
    bool degreeIdentityOK = false;
    bool lemmaContradiction = false;

    size_t firstNonunit = 0; // n means none
    size_t nonunitCount = 0;
    size_t lowerNonunitCount = 0;
    long lowerJunkDegree = 0;
    long topJunkDegree = -1;
    long totalJunkDegree = 0;
    long expectedJunkBudget = 0;
};

static JunkProfile buildJunkProfile(PolyRing &R,
                                    const Field &F,
                                    const std::vector<Polynomial> &expectedLargest,
                                    const ProbeResult &res,
                                    size_t n,
                                    size_t k)
{
    JunkProfile out;
    const size_t rank = k - 1;

    out.sStd = toStandardLength(R, expectedLargest, n);
    out.sigmaStd = toStandardLength(R, res.sigmaLargest, n);
    out.tStd.resize(n);
    out.tDefined.assign(n, false);

    R.assign(out.charA, productFactors(R, expectedLargest));
    R.assign(out.charAp, productFactors(R, res.sigmaLargest));
    R.assign(out.preservedOldProduct,
             productFactors(R, out.sStd, 0, n - rank));
    R.assign(out.junkProduct, R.one);
    R.assign(out.junkQuotient, R.zero);
    R.assign(out.repeatedFactorGCD, R.zero);

    out.fullDegreeOK = res.fullDegreeOK && (degPoly(R, out.charAp) == (long)n);
    out.rankOK = (res.actualRank == rank);
    out.shapeOK = true;
    out.coprimeOK = true;

    for (size_t i = 0; i < n; ++i) {
        Polynomial ti;
        bool defined = true;

        if (i < rank) {
            R.assign(ti, out.sigmaStd[i]);
        } else {
            const size_t oldIndex = i - rank;
            if (!dividesExact(ti, R, F, out.sigmaStd[i], out.sStd[oldIndex])) {
                defined = false;
                out.shapeOK = false;
                R.assign(ti, R.zero);
            }
        }

        out.tStd[i] = ti;
        out.tDefined[i] = defined;

        if (!defined) {
            out.coprimeOK = false;
            continue;
        }

        R.mulin(out.junkProduct, ti);
        const long d = degPoly(R, ti);
        if (d > 0) {
            if (out.nonunitCount == 0) out.firstNonunit = i;
            ++out.nonunitCount;
            if (i + 1 < n) ++out.lowerNonunitCount;
        }
        if (d >= 0) {
            out.totalJunkDegree += d;
            if (i + 1 < n) out.lowerJunkDegree += d;
            else out.topJunkDegree = d;
        }

        Polynomial g;
        R.assign(g, ti);
        R.gcdin(g, out.charA);
        if (!R.isOne(g)) out.coprimeOK = false;
    }

    if (out.nonunitCount == 0) out.firstNonunit = n;

    out.chainOK = out.shapeOK;
    if (out.chainOK) {
        for (size_t i = 0; i + 1 < n; ++i) {
            Polynomial q;
            if (!dividesExact(q, R, F, out.tStd[i + 1], out.tStd[i])) {
                out.chainOK = false;
                break;
            }
        }
    }

    out.charQuotientOK = exactQuotientMonic(out.junkQuotient, R, F,
                                             out.charAp,
                                             out.preservedOldProduct);
    out.productMatchesQuotient = out.shapeOK && out.charQuotientOK
                              && R.areEqual(out.junkProduct, out.junkQuotient);

    out.squarefreeDefined = out.shapeOK;
    if (out.squarefreeDefined) {
        Polynomial derivative;
        derivativePoly(derivative, R, F, out.junkProduct);
        R.assign(out.repeatedFactorGCD, out.junkProduct);
        R.gcdin(out.repeatedFactorGCD, derivative);
        out.squarefree = R.isOne(out.repeatedFactorGCD);
    }

    out.expectedJunkBudget = (long)prefixDegree(R, expectedLargest, rank);
    out.budgetOK = out.shapeOK && (out.totalJunkDegree == out.expectedJunkBudget);

    out.allLowerOne = out.shapeOK;
    if (out.allLowerOne) {
        for (size_t i = 0; i + 1 < n; ++i) {
            if (!R.isOne(out.tStd[i])) {
                out.allLowerOne = false;
                break;
            }
        }
    }

    const long sigmaTopDegree = degPoly(R, out.sigmaStd[n - 1]);
    const long expectedTopDegree = (long)prefixDegree(R, expectedLargest, k);
    out.degreeIdentityOK = (sigmaTopDegree == expectedTopDegree);

    out.lemmaContradiction = out.fullDegreeOK && out.rankOK && out.shapeOK
                          && out.coprimeOK && out.chainOK && out.squarefree
                          && !out.degreeIdentityOK;
    return out;
}

static void printPolyLine(PolyRing &R,
                          const std::string &name,
                          const Polynomial &f)
{
    std::cout << "  " << name << " = ";
    R.write(std::cout, f);
    std::cout << "\n";
}

static void printIndexedList(PolyRing &R,
                             const std::string &title,
                             const std::string &symbol,
                             const std::vector<Polynomial> &v,
                             size_t begin,
                             size_t end,
                             const std::vector<bool> *defined = nullptr)
{
    std::cout << "\n" << title << "\n";
    if (begin >= end || begin >= v.size()) {
        std::cout << "  (empty)\n";
        return;
    }
    end = std::min(end, v.size());

    size_t i = begin;
    while (i < end) {
        if (defined && !(*defined)[i]) {
            std::cout << "  " << symbol << "_" << (i + 1) << " = <undefined: exact division failed>\n";
            ++i;
            continue;
        }

        if (R.isOne(v[i])) {
            size_t j = i + 1;
            while (j < end && (!defined || (*defined)[j]) && R.isOne(v[j])) ++j;
            if (j - i >= 3) {
                std::cout << "  " << symbol << "_" << (i + 1)
                          << " ... " << symbol << "_" << j << " = 1\n";
            } else {
                for (size_t q = i; q < j; ++q)
                    std::cout << "  " << symbol << "_" << (q + 1) << " = 1\n";
            }
            i = j;
            continue;
        }

        std::cout << "  " << symbol << "_" << (i + 1) << " = ";
        R.write(std::cout, v[i]);
        std::cout << "    [deg " << degPoly(R, v[i]) << "]\n";
        ++i;
    }
}

static std::string classification(const JunkProfile &p)
{
    if (!p.rankOK) return "RANK_MISMATCH: resample the preconditioner";
    if (!p.fullDegreeOK) return "FULL_FNF_UNRELIABLE: charpoly degree is not n";
    if (!p.shapeOK) return "SHAPE_DIVISION_FAILURE";
    if (!p.charQuotientOK) return "CHARPOLY_JUNK_QUOTIENT_FAILURE";
    if (p.lemmaContradiction)
        return "DETERMINISTIC_LEMMA_CONTRADICTION: save this seed and instance";
    if (p.squarefree && p.allLowerOne && p.degreeIdentityOK)
        return "SQUAREFREE_CERTIFICATE_SUCCEEDS";
    if (!p.squarefree && p.allLowerOne && p.degreeIdentityOK)
        return "NONSQUAREFREE_BUT_HARMLESS: all junk is still in t_n";
    if (!p.allLowerOne && !p.degreeIdentityOK)
        return "LOWER_SLOT_JUNK_AND_DEGREE_IDENTITY_FAILURE";
    if (!p.allLowerOne && p.degreeIdentityOK)
        return "LOWER_SLOT_JUNK_BUT_IDENTITY_HOLDS: inspect budget/FNF";
    return "MIXED_OR_UNRESOLVED_DIAGNOSTIC";
}

static void printDiagnostic(const std::string &kind,
                            const std::string &caseName,
                            const char *algorithm,
                            PolyRing &R,
                            const std::vector<Polynomial> &expectedLargest,
                            size_t n,
                            size_t k,
                            const ProbeResult &res,
                            const JunkProfile &profile)
{
    const size_t rank = k - 1;
    const size_t survivingOld = n - rank;

    std::cout << "\n\n============================================================\n";
    std::cout << "JUNK INVARIANT PROFILE\n";
    std::cout << "kind=" << kind << "  case=" << caseName
              << "  algorithm=" << algorithm << "\n";
    std::cout << "n=" << n << "  selected k=" << k
              << "  requested rank=" << rank
              << "  actual rank(P)=" << res.actualRank << "\n";
    std::cout << "full-FNF attempts used=" << res.fnfAttemptsUsed
              << "  degree-n result=" << (res.fullDegreeOK ? "yes" : "no") << "\n";
    std::cout << "largest-first target: f_" << k
              << " = s_" << (n - k + 1) << "\n";
    std::cout << "============================================================\n";

    printIndexedList(R,
        "KNOWN INVARIANTS OF A (standard order; exact from construction)",
        "s", profile.sStd, 0, n);

    std::cout << "\n  Shape divisions use s_1,...,s_" << survivingOld << ".\n";
    std::cout << "  Junk budget uses s_" << (survivingOld + 1)
              << ",...,s_" << n << ".\n";

    printIndexedList(R,
        "ALL INVARIANTS OF A+P (standard order)",
        "sigma", profile.sigmaStd, 0, n);

    printIndexedList(R,
        "RECONSTRUCTED JUNK INVARIANTS",
        "t", profile.tStd, 0, n, &profile.tDefined);

    std::cout << "\nNONUNIT JUNK SUPPORT\n  ";
    bool first = true;
    for (size_t i = 0; i < n; ++i) {
        if (profile.tDefined[i] && !R.isOne(profile.tStd[i])) {
            if (!first) std::cout << ", ";
            std::cout << (i + 1);
            first = false;
        }
    }
    if (first) std::cout << "(none)";
    std::cout << "\n";

    std::cout << "\nFINAL REPORT\n";
    std::cout << "  rankOK                         = " << profile.rankOK << "\n";
    std::cout << "  fullFNF_degreeOK               = " << profile.fullDegreeOK << "\n";
    std::cout << "  shapeExactDivisionsOK          = " << profile.shapeOK << "\n";
    std::cout << "  charpolyJunkQuotientOK         = " << profile.charQuotientOK << "\n";
    std::cout << "  product(t_i)==charpolyQuotient = " << profile.productMatchesQuotient << "\n";
    std::cout << "  every t_i coprime to char(A)   = " << profile.coprimeOK << "\n";
    std::cout << "  junkChainOK                    = " << profile.chainOK << "\n";
    std::cout << "  junkProductSquarefree          = "
              << (profile.squarefreeDefined ? (profile.squarefree ? "1" : "0") : "NA") << "\n";
    std::cout << "  firstNonunitT                  = ";
    if (profile.firstNonunit == n) std::cout << "none\n";
    else std::cout << (profile.firstNonunit + 1) << "\n";
    std::cout << "  nonunitTCount                  = " << profile.nonunitCount << "\n";
    std::cout << "  lowerNonunitTCount             = " << profile.lowerNonunitCount << "\n";
    std::cout << "  lowerJunkDegree                = " << profile.lowerJunkDegree << "\n";
    std::cout << "  topJunkDegree                  = " << profile.topJunkDegree << "\n";
    std::cout << "  totalJunkDegree                = " << profile.totalJunkDegree << "\n";
    std::cout << "  expectedJunkBudget             = " << profile.expectedJunkBudget << "\n";
    std::cout << "  junkBudgetOK                   = " << profile.budgetOK << "\n";
    std::cout << "  t_1=...=t_{n-1}=1             = " << profile.allLowerOne << "\n";
    std::cout << "  degreeIdentityOK               = " << profile.degreeIdentityOK << "\n";
    std::cout << "  lemmaContradiction             = " << profile.lemmaContradiction << "\n";
    std::cout << "  classification                 = " << classification(profile) << "\n";

    if (profile.shapeOK) {
        printPolyLine(R, "T = product_i t_i", profile.junkProduct);
        printPolyLine(R, "gcd(T,T')", profile.repeatedFactorGCD);
    }
    if (profile.charQuotientOK)
        printPolyLine(R, "char(A+P) / product_{i<=n-r} s_i", profile.junkQuotient);
}

struct FinalRow {
    std::string algorithm;
    size_t actualRank = 0;
    bool fullOK = false;
    bool shapeOK = false;
    bool chainOK = false;
    bool coprimeOK = false;
    bool squarefree = false;
    long lowerJunkDegree = 0;
    bool degreeIdentityOK = false;
    bool lemmaContradiction = false;
    std::string label;
};

static void appendFinalRow(std::vector<FinalRow> &rows,
                           const char *algorithm,
                           const ProbeResult &res,
                           const JunkProfile &p)
{
    FinalRow row;
    row.algorithm = algorithm;
    row.actualRank = res.actualRank;
    row.fullOK = p.fullDegreeOK;
    row.shapeOK = p.shapeOK;
    row.chainOK = p.chainOK;
    row.coprimeOK = p.coprimeOK;
    row.squarefree = p.squarefreeDefined && p.squarefree;
    row.lowerJunkDegree = p.lowerJunkDegree;
    row.degreeIdentityOK = p.degreeIdentityOK;
    row.lemmaContradiction = p.lemmaContradiction;
    row.label = classification(p);
    rows.push_back(row);
}

int main(int argc, char **argv)
{
    commentator().setMaxDetailLevel(-1);
    commentator().setMaxDepth(-1);
    commentator().setReportStream(std::clog);

    int p = 10000019;
    int e = 1;
    int seed = (int)std::time(nullptr);
    int kArg = 0;
    int sArg = 0;
    int wArg = 0;
    int hArg = 0;
    int poly = 0;
    int spray = 1;
    int caseIndex = 0;
    int mixedIndex = 0;
    int fnfAttemptsArg = 3;
    char algoStr[16] = "";

    static Argument args[] = {
        { 'p', "-p P", "Field characteristic",                         TYPE_INT, &p },
        { 'e', "-e E", "Field extension degree",                       TYPE_INT, &e },
        { 'r', "-r R", "Random seed",                                  TYPE_INT, &seed },
        { 'k', "-k K", "Selected kth largest invariant (0=random)",     TYPE_INT, &kArg },
        { 'a', "-a A", "Algorithms: 0=Toeplitz 1=Butterfly 2=Dense",    TYPE_STR, algoStr },
        { 's', "-s S", "Custom: number of distinct block sizes",        TYPE_INT, &sArg },
        { 'w', "-w W", "Custom: repetitions per block size",            TYPE_INT, &wArg },
        { 'H', "-H H", "Custom: height/step between block sizes",       TYPE_INT, &hArg },
        { 'q', "-q Q", "Custom polynomial: 0=x, 1=x-1, 2=x+1",          TYPE_INT, &poly },
        { 'u', "-u U", "Apply random similarity spray: 0=no, 1=yes",    TYPE_INT, &spray },
        { 'c', "-c C", "Built-in Jordan case index",                    TYPE_INT, &caseIndex },
        { 'm', "-m M", "Mixed case: 0=off, 1..5=case index",            TYPE_INT, &mixedIndex },
        { 'z', "-z Z", "Maximum full-FNF attempts on same A+P",         TYPE_INT, &fnfAttemptsArg },
        END_OF_ARGUMENTS
    };

    parseArguments(argc, argv, args);
    std::srand(seed);

    const size_t s = sArg > 0 ? (size_t)sArg : 0;
    const size_t w = wArg > 0 ? (size_t)wArg : 0;
    const size_t h = hArg > 0 ? (size_t)hArg : 0;
    const size_t fnfAttempts = fnfAttemptsArg > 0 ? (size_t)fnfAttemptsArg : 1;
    const int algoMask = parseAlgoMask(algoStr[0] ? algoStr : nullptr);

    Field F(p, e);
    PolyRing R(F);

    SparseMat A(F, 0, 0);
    std::vector<Polynomial> expected;
    std::string kind;
    std::string caseName;

    const bool custom = (s > 0 && w > 0 && h > 0);
    if (custom) {
        const size_t n = w * h * s * (s + 1) / 2;
        A.resize(n, n);
        buildJordanMatrix(A, expected, F, R, s, w, h, poly);
        kind = "jordan";
        caseName = "custom_s" + std::to_string(s)
                 + "_w" + std::to_string(w)
                 + "_h" + std::to_string(h);
    } else if (mixedIndex > 0) {
        std::vector<MixedCase> mixedCases = defaultMixedCases();
        if ((size_t)mixedIndex > mixedCases.size()) {
            std::cerr << "Invalid -m index; valid range is 1.." << mixedCases.size() << "\n";
            return 2;
        }
        const MixedCase &tc = mixedCases[(size_t)mixedIndex - 1];
        std::vector<Polynomial> factors;
        for (const auto &pp : tc.profile) factors.push_back(makeMixedFactor(R, pp));
        size_t n = 0;
        for (const auto &f : factors) n += (size_t)R.deg(f);
        A.resize(n, n);
        buildCompanionInvariantMatrix(A, expected, F, R, factors);
        kind = "mixed";
        caseName = tc.name;
    } else {
        std::vector<TestParams> cases = defaultTestCases();
        if (cases.empty()) {
            std::cerr << "No built-in test cases available.\n";
            return 2;
        }
        if (caseIndex < 0 || (size_t)caseIndex >= cases.size()) {
            std::cerr << "Invalid -c index; valid range is 0.." << (cases.size() - 1) << "\n";
            return 2;
        }
        const TestParams &tc = cases[(size_t)caseIndex];
        const size_t n = tc.w * tc.h * tc.s * (tc.s + 1) / 2;
        A.resize(n, n);
        buildJordanMatrix(A, expected, F, R, tc.s, tc.w, tc.h, tc.poly);
        kind = "jordan";
        caseName = tc.name;
    }

    if (spray != 0) sprayMatrix(A, F, seed);

    if (expected.empty()) {
        std::cerr << "The selected matrix has no nontrivial invariant factors.\n";
        return 2;
    }

    size_t k = 1;
    if (kArg > 0) {
        k = (size_t)kArg;
        if (k > expected.size()) {
            std::cerr << "Selected k=" << k << " exceeds the known nontrivial count "
                      << expected.size() << ".\n";
            return 2;
        }
    } else if (expected.size() >= 2) {
        k = 2 + (size_t)(std::rand() % (int)(expected.size() - 1));
    }

    const size_t n = A.rowdim();
    const size_t requestedRank = k - 1;

    std::cout << "=== One-Shot Junk Invariant Profile ===\n";
    std::cout << "Field GF(" << p << "^" << e << ")  seed=" << seed << "\n";
    std::cout << "kind=" << kind << "  case=" << caseName << "  n=" << n << "\n";
    std::cout << "known nontrivial invariants=" << expected.size()
              << "  selected k=" << k
              << "  requested rank=" << requestedRank << "\n";
    std::cout << "algorithms:"
              << ((algoMask & 1) ? " Toeplitz" : "")
              << ((algoMask & 2) ? " Butterfly" : "")
              << ((algoMask & 4) ? " Dense" : "") << "\n";

    FrobeniusLarge<PolyRing> helper(R);
    std::vector<FinalRow> finalRows;

    if (algoMask & 1) {
        ProbeResult res;
        probeToeplitz(res, helper, F, R, A, k, fnfAttempts);
        JunkProfile profile = buildJunkProfile(R, F, expected, res, n, k);
        printDiagnostic(kind, caseName, algoName(1), R, expected, n, k, res, profile);
        appendFinalRow(finalRows, algoName(1), res, profile);
    }
    if (algoMask & 2) {
        ProbeResult res;
        probeButterfly(res, helper, F, R, A, k, fnfAttempts);
        JunkProfile profile = buildJunkProfile(R, F, expected, res, n, k);
        printDiagnostic(kind, caseName, algoName(2), R, expected, n, k, res, profile);
        appendFinalRow(finalRows, algoName(2), res, profile);
    }
    if (algoMask & 4) {
        ProbeResult res;
        probeDense(res, helper, F, R, A, k, fnfAttempts);
        JunkProfile profile = buildJunkProfile(R, F, expected, res, n, k);
        printDiagnostic(kind, caseName, algoName(4), R, expected, n, k, res, profile);
        appendFinalRow(finalRows, algoName(4), res, profile);
    }

    std::cout << "\n\n================ COMPARATIVE SUMMARY ================\n";
    std::cout << std::left
              << std::setw(22) << "algorithm"
              << std::setw(10) << "rank(P)"
              << std::setw(9)  << "fullOK"
              << std::setw(9)  << "shapeOK"
              << std::setw(9)  << "chainOK"
              << std::setw(10) << "coprime"
              << std::setw(10) << "sqfree"
              << std::setw(12) << "lowerDeg"
              << std::setw(10) << "degreeEq"
              << "classification\n";
    std::cout << std::string(130, '-') << "\n";

    bool contradiction = false;
    for (const auto &row : finalRows) {
        std::cout << std::left
                  << std::setw(22) << row.algorithm
                  << std::setw(10) << row.actualRank
                  << std::setw(9)  << row.fullOK
                  << std::setw(9)  << row.shapeOK
                  << std::setw(9)  << row.chainOK
                  << std::setw(10) << row.coprimeOK
                  << std::setw(10) << row.squarefree
                  << std::setw(12) << row.lowerJunkDegree
                  << std::setw(10) << row.degreeIdentityOK
                  << row.label << "\n";
        if (row.lemmaContradiction) contradiction = true;
    }

    std::cout << "=======================================================\n";
    return contradiction ? 3 : 0;
}
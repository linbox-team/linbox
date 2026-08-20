/* linbox/tests/test-degree-probe-lemma.C
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 * Written for probing the Villard/FrobeniusLarge degree claim.
 *
 * IMPORTANT INDEXING CONVENTION USED BY THIS TEST
 * ------------------------------------------------
 * This follows the LinBox/Gavin/FrobeniusLarge output convention:
 *
 *     f_1 = minpoly(A),   f_{i+1} | f_i,
 *
 * so f_k means the k-th largest invariant factor.  The algorithm call
 * kthInvariantFactor(..., k) uses a rank-(k-1) additive perturbation P.
 *
 * This test does NOT call kthInvariantFactorImpl directly, because that method
 * computes
 *
 *     g_raw = minpoly(A + P)
 *
 * and then immediately replaces it by
 *
 *     gcd(g_raw, minpoly(A)).
 *
 * Instead, this test builds the same three rank-(k-1) additive preconditioners
 * used by frobenius-large-generic.h and inspects the raw minimal polynomial
 *
 *     g = minpoly(A + P_{k-1})
 *
 * BEFORE the gcd step.
 *
 * For the k-th largest invariant f_k, this checks:
 *
 *     gcdOK:       gcd(g, f_1) == f_k
 *     prefixEq:    deg(g) == deg(f_1)+...+deg(f_k)
 *
 * To make the "extra factor" visible, it also prints:
 *
 *     extraDeg = deg(g) - deg(f_k)
 *     extraExp = deg(f_1)+...+deg(f_{k-1})
 *
 * When gcdOK=1, f_k divides g, so extraDeg is the degree of the quotient
 * g/f_k, i.e. the part killed by the later gcd step.  Thus the clean
 * degree-probe claim is equivalently:
 *
 *     extraDeg == extraExp.
 *
 * Useful failure line:
 *
 *     gcdOK=1 but extraEq=0
 *
 * means the kth-invariant probe succeeded but the degree-probe equality failed.
 *
 * Usage examples:
 *   make test-degree-probe-lemma
 *   ./tests/test-degree-probe-lemma -p 10000019 -k 8 -t 5 -a 012
 *   ./tests/test-degree-probe-lemma -s 5 -w 3 -H 2 -q 0 -k 10 -t 20
 *
 * Options:
 *   -p P       base prime, default 10000019
 *   -e E       extension degree, default 1
 *   -r R       random seed, default time(NULL)
 *   -k K       largest k to test, default = all nontrivial expected factors
 *   -t T       random trials per (case, algorithm, k), default 1
 *   -a A       algorithms: 0=Toeplitz(A+T1T2), 1=Butterfly(A+V^TDU), 2=Dense(A+UV), default 012
 *   -s S       custom matrix: number of distinct block sizes
 *   -w W       custom matrix: repetitions per block size
 *   -H H       custom matrix: height/step between block sizes
 *   -q Q       polynomial: 0=x, 1=x-1, 2=x+1, default 0
 *   -u U       apply random similarity spray: 0=no, 1=yes, default 1
 *   -m M       include mixed-primary companion-block stress cases: 0=no, 1=yes, default 1
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

template<class Blackbox>
static void rawProbeToeplitz(Polynomial &g,
                             FrobeniusLarge<PolyRing> &helper,
                             PolyRing &R,
                             const Blackbox &A,
                             size_t k)
{
    if (k <= 1) {
        helper.minpoly(g, A);
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

    typedef Sum<Blackbox, Prec> Preconditioned;
    Preconditioned Ap(A, P);

    helper.minpoly(g, Ap);
}

template<class Blackbox>
static void rawProbeButterfly(Polynomial &g,
                              FrobeniusLarge<PolyRing> &helper,
                              const Field &F,
                              const Blackbox &A,
                              size_t k)
{
    if (k <= 1) {
        helper.minpoly(g, A);
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

    // Same shape as frobenius-large-generic.h: D has rank nonzeros at the end.
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

    typedef Sum<Blackbox, Prec> Preconditioned;
    Preconditioned Ap(A, P);

    helper.minpoly(g, Ap);
}

template<class Blackbox>
static void rawProbeDense(Polynomial &g,
                          FrobeniusLarge<PolyRing> &helper,
                          const Field &F,
                          const Blackbox &A,
                          size_t k)
{
    if (k <= 1) {
        helper.minpoly(g, A);
        return;
    }

    const size_t n = A.rowdim();
    const size_t rank = k - 1;

    Field::RandIter RI(F);
    DenseMat U(F, n, rank);
    DenseMat V(F, rank, n);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < rank; ++j) {
            Element e;
            RI.random(e);
            U.setEntry(i, j, e);
        }
    }
    for (size_t i = 0; i < rank; ++i) {
        for (size_t j = 0; j < n; ++j) {
            Element e;
            RI.random(e);
            V.setEntry(i, j, e);
        }
    }

    typedef Compose<DenseMat, DenseMat> Prec;
    Prec P(U, V);

    typedef Sum<Prec, Blackbox> Preconditioned;
    Preconditioned Ap(P, A);

    helper.minpoly(g, Ap);
}

struct Counts {
    size_t total = 0;
    size_t jordanTotal = 0;
    size_t mixedTotal = 0;
    size_t gcdOK = 0;
    size_t prefixLe = 0;
    size_t prefixEq = 0;
    size_t extraEq = 0;
    size_t successButPrefixNe = 0;
    size_t successButExtraNe = 0;
    size_t successButViolatesLe = 0;
};

static void updateAndPrint(const std::string &kind,
                           const std::string &caseName,
                           const char *alg,
                           size_t trial,
                           size_t k,
                           PolyRing &R,
                           const Polynomial &mA,
                           const std::vector<Polynomial> &expected,
                           const Polynomial &g,
                           Counts &counts)
{
    Polynomial gcd;
    R.assign(gcd, g);
    R.gcdin(gcd, mA);

    Polynomial fk = expectedAt(R, expected, k);

    const long rawDeg    = degPoly(R, g);       // deg minpoly(A+P)
    const long fkDeg     = degPoly(R, fk);      // deg f_k, with deg(1)=0 in practice
    const long prefixDeg = (long) prefixDegree(R, expected, k);
    const long extraDeg  = rawDeg - fkDeg;      // degree of g/f_k when gcdOK=1
    const long extraExp  = prefixDeg - fkDeg;   // sum_{j<k} deg f_j
    const size_t rank    = (k == 0) ? 0 : k - 1;

    const bool invOK = R.areEqual(gcd, fk);
    const bool leOK  = (rawDeg >= 0 && rawDeg <= prefixDeg);
    const bool eqOK  = (rawDeg >= 0 && rawDeg == prefixDeg);
    const bool exOK  = (rawDeg >= 0 && extraDeg == extraExp);

    counts.total++;
    if (kind == "mixed") counts.mixedTotal++;
    else                 counts.jordanTotal++;
    if (invOK) counts.gcdOK++;
    if (leOK)  counts.prefixLe++;
    if (eqOK)  counts.prefixEq++;
    if (exOK)  counts.extraEq++;
    if (invOK && !eqOK) counts.successButPrefixNe++;
    if (invOK && !exOK) counts.successButExtraNe++;
    if (invOK && !leOK) counts.successButViolatesLe++;

    std::cout << std::left
              << std::setw(9)  << kind
              << std::setw(22) << caseName
              << std::setw(21) << alg
              << std::setw(6)  << k
              << std::setw(6)  << rank
              << std::setw(7)  << trial
              << std::setw(8)  << rawDeg
              << std::setw(8)  << fkDeg
              << std::setw(10) << extraDeg
              << std::setw(10) << extraExp
              << std::setw(10) << prefixDeg
              << std::setw(8)  << (leOK  ? "1" : "0")
              << std::setw(9)  << (eqOK  ? "1" : "0")
              << std::setw(9)  << (exOK  ? "1" : "0")
              << std::setw(8)  << (invOK ? "1" : "0");

    if (invOK && !exOK)       std::cout << "  SUCCESS_BUT_EXTRA_NE_EXPECTED";
    else if (invOK && !eqOK)  std::cout << "  SUCCESS_BUT_DEGREE_NE_PREFIX";
    else if (!invOK)          std::cout << "  unlucky_or_wrong_probe";
    std::cout << "\n";

    if (invOK && (!eqOK || !exOK)) {
        std::cout << "    g=minpoly(A+P): "; R.write(std::cout, g); std::cout << "\n";
        std::cout << "    gcd(g,f1):      "; R.write(std::cout, gcd); std::cout << "\n";
        std::cout << "    expected f_k:   "; R.write(std::cout, fk); std::cout << "\n";
        std::cout << "    degrees: raw=" << rawDeg
                  << " fk=" << fkDeg
                  << " extra=" << extraDeg
                  << " expectedExtra=" << extraExp
                  << " prefix=" << prefixDeg << "\n";
    }
}

static void runOneCase(const TestParams &tc,
                       const Field &F,
                       PolyRing &R,
                       int seed,
                       bool doSpray,
                       size_t kmaxArg,
                       size_t trials,
                       int algoMask,
                       Counts &counts)
{
    const size_t n = tc.w * tc.h * tc.s * (tc.s + 1) / 2;

    SparseMat A(F, n, n);
    std::vector<Polynomial> expected;
    buildJordanMatrix(A, expected, F, R, tc.s, tc.w, tc.h, tc.poly);
    if (doSpray) sprayMatrix(A, F, seed);

    FrobeniusLarge<PolyRing> helper(R);

    Polynomial mA;
    helper.minpoly(mA, A);

    // Default: all known nontrivial factors.  Pass -k larger to probe trailing 1s too.
    const size_t kmax = (kmaxArg == 0) ? expected.size() : kmaxArg;

    std::cout << "\nCase " << tc.name
              << "  n=" << n
              << "  expected_nontrivial=" << expected.size()
              << "  kmax=" << kmax
              << "  minpoly_deg=" << R.deg(mA)
              << "\n";

    for (size_t k = 1; k <= kmax; ++k) {
        for (size_t t = 1; t <= trials; ++t) {
            Polynomial g;

            if (algoMask & 1) {
                rawProbeToeplitz(g, helper, R, A, k);
                updateAndPrint("jordan", tc.name, algoName(1), t, k, R, mA, expected, g, counts);
            }
            if (algoMask & 2) {
                rawProbeButterfly(g, helper, F, A, k);
                updateAndPrint("jordan", tc.name, algoName(2), t, k, R, mA, expected, g, counts);
            }
            if (algoMask & 4) {
                rawProbeDense(g, helper, F, A, k);
                updateAndPrint("jordan", tc.name, algoName(4), t, k, R, mA, expected, g, counts);
            }
        }
    }
}


static void runMixedCase(const MixedCase &tc,
                         const Field &F,
                         PolyRing &R,
                         int seed,
                         bool doSpray,
                         size_t kmaxArg,
                         size_t trials,
                         int algoMask,
                         Counts &counts)
{
    std::vector<Polynomial> factors;
    for (const auto &pp : tc.profile)
        factors.push_back(makeMixedFactor(R, pp));

    size_t n = 0;
    for (const auto &f : factors) n += (size_t) R.deg(f);

    SparseMat A(F, n, n);
    std::vector<Polynomial> expected;
    buildCompanionInvariantMatrix(A, expected, F, R, factors);
    if (doSpray) sprayMatrix(A, F, seed);

    FrobeniusLarge<PolyRing> helper(R);

    Polynomial mA;
    helper.minpoly(mA, A);

    const size_t kmax = (kmaxArg == 0) ? expected.size() : kmaxArg;

    std::cout << "\nCase " << tc.name
              << "  n=" << n
              << "  expected_nontrivial=" << expected.size()
              << "  kmax=" << kmax
              << "  minpoly_deg=" << R.deg(mA)
              << "  type=mixed-primary companion blocks"
              << "\n";

    std::cout << "  expected factors, largest-first:\n";
    for (size_t i = 0; i < expected.size(); ++i) {
        std::cout << "    f_" << (i + 1) << " = ";
        R.write(std::cout, expected[i]);
        std::cout << "  deg=" << R.deg(expected[i]) << "\n";
    }

    for (size_t k = 1; k <= kmax; ++k) {
        for (size_t t = 1; t <= trials; ++t) {
            Polynomial g;

            if (algoMask & 1) {
                rawProbeToeplitz(g, helper, R, A, k);
                updateAndPrint("mixed", tc.name, algoName(1), t, k, R, mA, expected, g, counts);
            }
            if (algoMask & 2) {
                rawProbeButterfly(g, helper, F, A, k);
                updateAndPrint("mixed", tc.name, algoName(2), t, k, R, mA, expected, g, counts);
            }
            if (algoMask & 4) {
                rawProbeDense(g, helper, F, A, k);
                updateAndPrint("mixed", tc.name, algoName(4), t, k, R, mA, expected, g, counts);
            }
        }
    }
}

int main(int argc, char **argv)
{
    commentator().setMaxDetailLevel(-1);
    commentator().setMaxDepth(-1);
    commentator().setReportStream(std::clog);

    int    p        = 10000019;
    int    e        = 1;
    int    seed     = (int) std::time(nullptr);
    // LinBox parseArguments(TYPE_INT) expects int*, not size_t*.
    // Keep the command-line storage as int and convert to size_t after parsing.
    // Passing size_t* here can silently corrupt neighboring option values and
    // was the reason mixed-primary cases could appear to be skipped.
    int    kmaxArg  = 0;
    int    trialsArg= 1;
    int    sArg     = 0;
    int    wArg     = 0;
    int    hArg     = 0;
    int    poly     = 0;
    int    spray    = 1;
    int    mixed    = 1;
    char   algoStr[16] = "";

    static Argument args[] = {
        { 'p', "-p P", "Field characteristic",                         TYPE_INT, &p },
        { 'e', "-e E", "Field extension degree",                       TYPE_INT, &e },
        { 'r', "-r R", "Random seed",                                  TYPE_INT, &seed },
        { 'k', "-k K", "Max k to test (0 = expected nontrivial count)", TYPE_INT, &kmaxArg },
        { 't', "-t T", "Trials per (case, algorithm, k)",              TYPE_INT, &trialsArg },
        { 'a', "-a A", "Algorithms: 0=Toeplitz 1=Butterfly 2=Dense",    TYPE_STR, algoStr },
        { 's', "-s S", "Custom: number of distinct block sizes",        TYPE_INT, &sArg },
        { 'w', "-w W", "Custom: repetitions per block size",            TYPE_INT, &wArg },
        { 'H', "-H H", "Custom: height/step between block sizes",       TYPE_INT, &hArg },
        { 'q', "-q Q", "Custom polynomial: 0=x, 1=x-1, 2=x+1",          TYPE_INT, &poly },
        { 'u', "-u U", "Apply random similarity spray: 0=no, 1=yes",    TYPE_INT, &spray },
        { 'm', "-m M", "Include mixed-primary stress cases: 0=no, 1=yes", TYPE_INT, &mixed },
        END_OF_ARGUMENTS
    };

    parseArguments(argc, argv, args);
    std::srand(seed);

    const size_t kmax   = (kmaxArg   <= 0) ? 0 : (size_t) kmaxArg;
    const size_t trials = (trialsArg <= 0) ? 1 : (size_t) trialsArg;
    const size_t s      = (sArg      <= 0) ? 0 : (size_t) sArg;
    const size_t w      = (wArg      <= 0) ? 0 : (size_t) wArg;
    const size_t h      = (hArg      <= 0) ? 0 : (size_t) hArg;

    Field F(p, e);
    PolyRing R(F);

    int algoMask = parseAlgoMask(algoStr[0] ? algoStr : nullptr);

    std::vector<TestParams> cases;
    const bool customCase = (s > 0 && w > 0 && h > 0);
    if (customCase) {
        cases.push_back({s, w, h, poly,
                         "custom_s" + std::to_string(s)
                       + "_w" + std::to_string(w)
                       + "_h" + std::to_string(h)});
    } else {
        cases = defaultTestCases();
    }

    std::vector<MixedCase> mixedCases = defaultMixedCases();
    const bool runMixed = (mixed != 0);

    std::cout << "=== Degree Probe Lemma Test: raw minpoly + explicit extra degree ===\n";
    std::cout << "Convention: f_1=minpoly(A), k=index of kth largest factor, rank(P)=k-1\n";
    std::cout << "Field GF(" << p << "^" << e << ")"
              << "  seed=" << seed
              << "  trials=" << trials
              << "  spray=" << spray
              << "  mixed=" << mixed
              << "\n";
    std::cout << "Case plan: jordan_cases=" << cases.size()
              << "  mixed_cases=" << (runMixed ? mixedCases.size() : 0)
              << "  custom_jordan=" << (customCase ? 1 : 0)
              << "\n";
    std::cout << "Algorithms:"
              << ((algoMask & 1) ? " Toeplitz(A+T1T2)" : "")
              << ((algoMask & 2) ? " Butterfly(A+VtDU)" : "")
              << ((algoMask & 4) ? " Dense(A+UV)" : "")
              << "\n\n";

    std::cout << std::left
              << std::setw(9)  << "kind"
              << std::setw(22) << "case"
              << std::setw(21) << "algorithm"
              << std::setw(6)  << "k"
              << std::setw(6)  << "rank"
              << std::setw(7)  << "trial"
              << std::setw(8)  << "deg(g)"
              << std::setw(8)  << "degfk"
              << std::setw(10) << "extraDeg"
              << std::setw(10) << "extraExp"
              << std::setw(10) << "prefixDeg"
              << std::setw(8)  << "deg<="
              << std::setw(9)  << "prefixEq"
              << std::setw(9)  << "extraEq"
              << std::setw(8)  << "gcdOK"
              << "note\n";
    std::cout << std::string(160, '-') << "\n";

    Counts counts;
    for (const auto &tc : cases)
        runOneCase(tc, F, R, seed, spray != 0, kmax, trials, algoMask, counts);

    // Mixed-primary companion-block stress cases are controlled only by -m.
    // They now also run alongside a custom -s/-w/-H case unless -m 0 is passed.
    if (runMixed) {
        for (const auto &tc : mixedCases)
            runMixedCase(tc, F, R, seed, spray != 0, kmax, trials, algoMask, counts);
    }

    std::cout << "\n=== Summary ===\n";
    std::cout << "total probes:                         " << counts.total << "\n";
    std::cout << "jordan-suite probes:                  " << counts.jordanTotal << "\n";
    std::cout << "mixed-primary probes:                 " << counts.mixedTotal << "\n";
    std::cout << "gcd success probes:                   " << counts.gcdOK << "\n";
    std::cout << "degree <= prefix probes:              " << counts.prefixLe << "\n";
    std::cout << "degree == prefix probes:              " << counts.prefixEq << "\n";
    std::cout << "extraDeg == expected extra probes:    " << counts.extraEq << "\n";
    std::cout << "success but degree != prefix:         " << counts.successButPrefixNe << "\n";
    std::cout << "success but extraDeg != expected:     " << counts.successButExtraNe << "\n";
    std::cout << "success but degree > prefix:          " << counts.successButViolatesLe << "\n";

    // Exit 2 means: the kth-invariant probe succeeded somewhere but the equality
    // version of the degree-probe lemma failed there.  Exit 1 means no equality
    // counterexample was seen, but some probes were unlucky or violated the bound.
    if (counts.successButPrefixNe > 0 || counts.successButExtraNe > 0) return 2;
    if (counts.gcdOK != counts.total || counts.prefixLe != counts.total) return 1;
    return 0;
}

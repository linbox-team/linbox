/* linbox/tests/test-frobenius-suite-input.C
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 * Written by Omesh Dhar Dwivedi <odd23@drexel.edu>
 *
 * Speed benchmark: read a matrix from a given path in SMS format and compare
 * Toeplitz, Butterfly, and Dense preconditioners under binary, exponential,
 * and factor-aware threshold search, plus LIFs.
 * When Dense and LIFs are selected together, compare their leading factors
 * after timing.
 *
 * Usage:
 *   ./test-frobenius-suite-input -f /path/to/matrix.sms [-p 10007] [-k 0]
 *       [-c 0.999999] [-b 0] [-r 42] -a 0123456789 [-n 5] [-x] [-i] [-P]
 *
 * -a accepts digits:
 *   0=Toeplitz  1=ToeplitzSearch  2=Butterfly  3=ButterflySearch
 *   4=Dense     5=DenseSearch     6=LIFs
 *   7=ToeplitzFactorAware  8=ButterflyFactorAware  9=DenseFactorAware
 *   default with -x: all ten
 *
 * LIFs block selection:
 *   -c C sets the target success probability used to choose the block size.
 *   -b B forces block size B; -b 0 uses the probability-based choice.
 *
 * Safety guardrails:
 *   -a must be supplied explicitly; omitting it would run all ten methods.
 *   -k 0 and multiple runs on enormous matrices require -x.
 *   -x acknowledges risky settings; it never permits -n 0.
 *   This benchmark performs no warmup calls.
 *
 * No invariant profile is printed by default.
 * Use -i for grouped multiplicities and degrees.
 * Use -P for grouped multiplicities, degrees, and complete polynomials;
 * -P implies -i.
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

struct RunTiming {
    double wall;
    double user;
    double system;
};

struct TimingTotals {
    double wall;
    double user;
    double system;
};

void printRunTiming(
    const char *label,
    int run,
    int nruns,
    const RunTiming &timing)
{
    std::cout << "Timing"
              << "  algorithm=" << label
              << "  run=" << (run + 1) << "/" << nruns
              << "  wall_time=" << std::fixed << std::setprecision(6)
              << timing.wall << "s"
              << "  user_time=" << timing.user << "s"
              << "  system_time=" << timing.system << "s"
              << std::endl;
}

TimingTotals timingTotals(const std::vector<RunTiming> &timings)
{
    TimingTotals totals = { 0.0, 0.0, 0.0 };

    for (size_t i = 0; i < timings.size(); ++i) {
        totals.wall += timings[i].wall;
        totals.user += timings[i].user;
        totals.system += timings[i].system;
    }

    return totals;
}

void printTimingSummary(
    const char *label,
    const std::vector<RunTiming> &timings,
    const TimingTotals &totals)
{
    const double runs = static_cast<double>(timings.size());

    std::cout << "Timing"
              << "  algorithm=" << label
              << "  summary"
              << "  avg_wall_time=" << std::fixed << std::setprecision(6)
              << totals.wall / runs << "s"
              << "  total_wall_time=" << totals.wall << "s"
              << "  avg_user_time=" << totals.user / runs << "s"
              << "  total_user_time=" << totals.user << "s"
              << "  avg_system_time=" << totals.system / runs << "s"
              << "  total_system_time=" << totals.system << "s"
              << std::endl;
}

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

    std::vector<Polynomial> fs;

    std::vector<RunTiming> timings;
    timings.reserve(static_cast<size_t>(nruns));
    size_t nfactors = 0;
    for (int run = 0; run < nruns; ++run) {
        fs.clear();
        Givaro::Timer T; T.clear(); T.start();
        solver.frobeniusInvariants(fs, M, keff);
        T.stop();

        RunTiming timing;
        timing.wall = T.realtime();
        timing.user = T.usertime();
        timing.system = T.systime();
        timings.push_back(timing);
        printRunTiming(label, run, nruns, timing);
        nfactors = fs.size();
    }
    const TimingTotals totals = timingTotals(timings);
    printTimingSummary(label, timings, totals);

    std::cout << std::left << std::setw(22) << label
              << "  factors=" << std::setw(6) << nfactors;
    std::cout << "  avg_time=" << std::fixed << std::setprecision(6)
              << totals.user / nruns << "s"
              << "  total=" << totals.user << "s" << std::endl;

    // The final timed result is returned only after the timer has stopped.
    return fs;
}

// LIFs has a separate timed path because frobeniusInvariants() hardcodes
// b = k+2. largestInvariantFactors() accepts the selected block size directly
// but returns factors in the opposite order, so reverse inside the timed region
// to preserve the old frobeniusInvariants() behavior exactly.
template<typename Solver, typename SparseMat, typename PolyRing>
std::vector<typename PolyRing::Element> benchLIFs(
    Solver &solver,
    SparseMat &M,
    size_t requestedK,
    size_t blockSize,
    int nruns,
    const char *label)
{
    typedef typename PolyRing::Element Polynomial;

    std::vector<Polynomial> fs;

    std::vector<RunTiming> timings;
    timings.reserve(static_cast<size_t>(nruns));
    size_t returnedEntries = 0;
    for (int run = 0; run < nruns; ++run) {
        fs.clear();

        Givaro::Timer T; T.clear(); T.start();
        solver.largestInvariantFactors(fs, M, blockSize);
        std::reverse(fs.begin(), fs.end());
        T.stop();

        RunTiming timing;
        timing.wall = T.realtime();
        timing.user = T.usertime();
        timing.system = T.systime();
        timings.push_back(timing);
        printRunTiming(label, run, nruns, timing);
        returnedEntries = fs.size();
    }

    const TimingTotals totals = timingTotals(timings);
    printTimingSummary(label, timings, totals);

    std::cout << std::left << std::setw(22) << label
              << "  returned_entries=" << std::setw(6) << returnedEntries
              << "  requested_k=" << std::setw(6) << requestedK
              << "  block_b=" << std::setw(6) << blockSize
              << "  avg_time=" << std::fixed << std::setprecision(6)
              << totals.user / nruns << "s"
              << "  total=" << totals.user << "s" << std::endl;

    // The final timed result is returned only after the timer has stopped.
    return fs;
}

// Compare the requested leading factors. LIFs may return b >= k entries, while
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

// Print the requested leading profile in LinBox's nonincreasing invariant-
// factor order. Equal factors are contiguous, so no sorting is needed.
// Some FrobeniusLarge variants omit trailing unit factors, while LIFs may
// return more than requested; normalize both cases to exactly k entries.
template<typename PolyRing>
void printInvariantProfile(
    const std::vector<typename PolyRing::Element> &fs,
    size_t k,
    const PolyRing &R,
    const std::string &source,
    bool printPolynomials)
{
    const size_t explicitCount = std::min(k, fs.size());

    std::cout << std::string(70, '-') << std::endl;
    std::cout << "Invariant profile from " << source
              << " (first " << k << ", nonincreasing order)"
              << std::endl;

    size_t uniqueIndex = 1;
    size_t i = 0;
    while (i < explicitCount) {
        size_t j = i + 1;
        while (j < explicitCount && fs[j] == fs[i]) ++j;

        size_t copies = j - i;
        if (j == explicitCount && explicitCount < k && R.isOne(fs[i]))
            copies += k - explicitCount;

        std::cout << "unique=" << std::setw(4) << uniqueIndex
                  << "  copies=" << std::setw(8) << copies
                  << "  degree=" << R.deg(fs[i]);
        if (printPolynomials)
            R.write(std::cout << "  invariant=", fs[i]);
        std::cout << std::endl;

        ++uniqueIndex;
        i = j;
    }

    // Missing entries denote trailing unit invariant factors. Do not print a
    // second unit row if the preceding explicit group was already the unit.
    if (explicitCount < k &&
        (explicitCount == 0 || !R.isOne(fs[explicitCount - 1]))) {
        std::cout << "unique=" << std::setw(4) << uniqueIndex
                  << "  copies=" << std::setw(8) << (k - explicitCount)
                  << "  degree=0";
        if (printPolynomials)
            std::cout << "  invariant=1";
        std::cout << std::endl;
    }

    std::cout << "Profile entries: " << k
              << " (solver returned " << fs.size()
              << "; used " << explicitCount;
    if (explicitCount < k)
        std::cout << "; implicit units " << (k - explicitCount);
    else if (fs.size() > k)
        std::cout << "; ignored oversampling entries " << (fs.size() - k);
    std::cout << ")" << std::endl;

    if (!printPolynomials)
        std::cout << "Full polynomials hidden; rerun with -P to print them."
                  << std::endl;
}

// bit0=Toeplitz, bit1=ToeplitzSearch, bit2=Butterfly, bit3=ButterflySearch,
// bit4=Dense, bit5=DenseSearch, bit6=LIFs, bit7=ToeplitzFactorAware,
// bit8=ButterflyFactorAware, bit9=DenseFactorAware
int parseAlgoMask(const char *s) {
    if (s == nullptr || strlen(s) == 0) return 1023; // all ten
    int mask = 0;
    for (size_t i = 0; i < strlen(s); ++i) {
        if (s[i] == '0') mask |= 1;
        if (s[i] == '1') mask |= 2;
        if (s[i] == '2') mask |= 4;
        if (s[i] == '3') mask |= 8;
        if (s[i] == '4') mask |= 16;
        if (s[i] == '5') mask |= 32;
        if (s[i] == '6') mask |= 64;
        if (s[i] == '7') mask |= 128;
        if (s[i] == '8') mask |= 256;
        if (s[i] == '9') mask |= 512;
    }
    return mask ? mask : 1023;
}

bool validAlgoSelection(const std::string &s) {
    if (s.empty()) return false;
    for (size_t i = 0; i < s.size(); ++i)
        if (s[i] < '0' || s[i] > '9') return false;
    return true;
}

int main(int argc, char **argv) {

    commentator().setMaxDetailLevel(-1);
    commentator().setMaxDepth(-1);
    commentator().setReportStream(std::clog);

    int         p       = 10007;
    int         e       = 1;
    int         k       = 0;
    int         seed    = time(NULL);
    int         nruns   = 5;
    double      lifsConfidence = 0.999999;
    int         forcedLifsBlock = 0;
    bool        allowRisky = false;
    bool        printProfile = false;
    bool        printPolynomials = false;
    std::string filepath;
    std::string algoStr;

    static Argument args[] = {
        { 'f', "-f F", "Input matrix file path (SMS format)",
            TYPE_STR, &filepath },
        { 'p', "-p P", "Field characteristic",
            TYPE_INT, &p },
        { 'e', "-e E", "Field extension degree",
            TYPE_INT, &e },
        { 'k', "-k K", "Number of invariant factors (0 = log2(n), requires -x)",
            TYPE_INT, &k },
        { 'c', "-c C", "LIFs target success probability when -b 0 (default: 0.999999)",
            TYPE_DOUBLE, &lifsConfidence },
        { 'b', "-b B", "LIFs block size (0 = choose automatically from -c)",
            TYPE_INT, &forcedLifsBlock },
        { 'r', "-r R", "Random seed",
            TYPE_INT, &seed },
        { 'n', "-n N", "Number of timed runs",
            TYPE_INT, &nruns },
        { 'a', "-a A", "Algorithms: 0=Toeplitz 1=ToeplitzSearch 2=Butterfly 3=ButterflySearch 4=Dense 5=DenseSearch 6=LIFs 7=ToeplitzFactorAware 8=ButterflyFactorAware 9=DenseFactorAware (required unless -x)",
            TYPE_STR, &algoStr },
        { 'x', "-x", "Acknowledge and allow risky benchmark settings",
            TYPE_BOOL, &allowRisky },
        { 'i', "-i", "Print grouped invariant multiplicities and degrees",
            TYPE_BOOL, &printProfile },
        { 'P', "-P", "Print grouped invariant multiplicities, degrees, and complete polynomials (implies -i)",
            TYPE_BOOL, &printPolynomials },
        END_OF_ARGUMENTS
    };

    parseArguments(argc, argv, args);
    if (printPolynomials) printProfile = true;
    srand(seed);

    if (nruns <= 0) {
        std::cerr << "Error: number of timed runs must be positive." << std::endl;
        return -1;
    }
    if (k < 0) {
        std::cerr << "Error: requested invariant-factor count cannot be negative."
                  << std::endl;
        return -1;
    }
    if (forcedLifsBlock < 0) {
        std::cerr << "Error: LIFs block size cannot be negative." << std::endl;
        return -1;
    }
    if (algoStr.empty()) {
        if (!allowRisky) {
            std::cerr << "Safety stop: -a was omitted, which would run all ten "
                      << "algorithms. Specify -a explicitly (for example, -a 46), "
                      << "or rerun with -x to acknowledge this workload."
                      << std::endl;
            return -1;
        }
        std::cerr << "Safety warning: -x accepted omitted -a; all ten algorithms "
                  << "will run." << std::endl;
    }
    else if (!validAlgoSelection(algoStr)) {
        std::cerr << "Error: -a must contain only algorithm digits 0 through 9."
                  << std::endl;
        return -1;
    }

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

    if (n == 0) {
        std::cerr << "Error: input matrix has dimension zero." << std::endl;
        return -1;
    }
    if (M.coldim() != n) {
        std::cerr << "Error: input matrix must be square (rows=" << n
                  << ", columns=" << M.coldim() << ")." << std::endl;
        return -1;
    }

    size_t keff = (k == 0)
        ? (size_t)std::max(1.0, std::log2((double)n))
        : static_cast<size_t>(k);
    int algoMask = parseAlgoMask(algoStr.empty() ? nullptr : algoStr.c_str());

    if (k == 0) {
        if (!allowRisky) {
            std::cerr << "Safety stop: -k 0 selects k=" << keff
                      << " (approximately log2(n)) for n=" << n
                      << ". Specify -k explicitly, or rerun with -x to "
                      << "acknowledge this workload." << std::endl;
            return -1;
        }
        std::cerr << "Safety warning: -x accepted -k 0; the benchmark will request "
                  << keff << " invariant factors." << std::endl;
    }

    const size_t enormousMatrixOrder = 100000;
    const size_t enormousMatrixNnz = 10000000;
    const bool enormousMatrix =
        n >= enormousMatrixOrder || nnz >= enormousMatrixNnz;
    if (enormousMatrix && nruns > 1) {
        if (!allowRisky) {
            std::cerr << "Safety stop: this matrix is enormous (n=" << n
                      << ", nnz=" << nnz << "; threshold n>="
                      << enormousMatrixOrder << " or nnz>=" << enormousMatrixNnz
                      << ") and -n " << nruns
                      << " requests multiple timed runs. Use -n 1, or rerun "
                      << "with -x to acknowledge the potentially multi-day job."
                      << std::endl;
            return -1;
        }
        std::cerr << "Safety warning: -x accepted " << nruns
                  << " timed runs on an enormous matrix (n=" << n
                  << ", nnz=" << nnz << ")."
                  << std::endl;
    }

    if (keff == 0 || keff > n) {
        std::cerr << "Error: requested invariant-factor count must satisfy 1 <= k <= n."
                  << std::endl;
        return -1;
    }

    InvariantFactors<Field, PolyRing> IFD(F, R);
    size_t lifsBlock = 0;
    if (algoMask & 64) {
        if (forcedLifsBlock == 0 &&
            !(lifsConfidence > 0.0 && lifsConfidence < 1.0)) {
            std::cerr << "Error: LIFs target success probability must satisfy 0 < -c < 1."
                      << std::endl;
            return -1;
        }

        if (forcedLifsBlock == 0) {
            // Reproduce InvariantFactors::min_block_size here because that
            // helper assigns the field's Givaro::Integer cardinality directly
            // to size_t, which is ambiguous for NTL_zz_p on some platforms.
            // The output-parameter overload performs the conversion safely.
            size_t q = 0;
            F.cardinality(q);
            if (q < 2) {
                std::cerr << "Error: field cardinality must be at least 2."
                          << std::endl;
                return -1;
            }

            const double probabilityConstant = (q == 2) ? 3.0 : 2.0;
            const size_t oversampling = static_cast<size_t>(std::ceil(
                std::log(probabilityConstant /
                         (1.0 - std::sqrt(lifsConfidence))) /
                std::log(static_cast<double>(q))));
            if (oversampling > n - keff) {
                std::cerr << "Error: probability-selected LIFs block size exceeds "
                          << "matrix dimension n=" << n << "." << std::endl;
                return -1;
            }
            lifsBlock = keff + oversampling;
        }
        else {
            lifsBlock = static_cast<size_t>(forcedLifsBlock);
        }

        if (lifsBlock < keff) {
            std::cerr << "Error: LIFs block size b=" << lifsBlock
                      << " is smaller than requested k=" << keff << "."
                      << std::endl;
            return -1;
        }
        if (lifsBlock > n) {
            std::cerr << "Error: LIFs block size b=" << lifsBlock
                      << " exceeds matrix dimension n=" << n << "."
                      << std::endl;
            return -1;
        }
    }

    std::cout << "=== Frobenius Speed Benchmark ===" << std::endl;
    std::cout << "File:      " << filepath << std::endl;
    std::cout << "n:         " << n        << std::endl;
    std::cout << "nnz:       " << nnz      << std::endl;
    std::cout << "Field:     GF(" << p << "^" << e << ")" << std::endl;
    std::cout << "k:         " << keff << (k == 0 ? " (log2(n) default)" : "") << std::endl;
    std::cout << "seed:      " << seed     << std::endl;
    std::cout << "runs:      " << nruns << std::endl;
    std::cout << "warmups:   0 (disabled)" << std::endl;
    std::cout << "timing:    wall/real elapsed time is primary; user and system CPU "
              << "are also recorded per run and in aggregate" << std::endl;
    std::cout << "legacy:    avg_time and total remain user CPU for compatibility"
              << std::endl;
    std::cout << "profile:    "
              << (printPolynomials ? "full polynomials (-P)"
                  : (printProfile ? "multiplicity + degree (-i)" : "off"))
              << std::endl;
    std::cout << "Algorithms:"
              << ((algoMask &  1) ? " Toeplitz"        : "")
              << ((algoMask &  2) ? " ToeplitzSearch"  : "")
              << ((algoMask &  4) ? " Butterfly"       : "")
              << ((algoMask &  8) ? " ButterflySearch" : "")
              << ((algoMask & 16) ? " Dense"           : "")
              << ((algoMask & 32) ? " DenseSearch"     : "")
              << ((algoMask & 64) ? " LIFs"            : "")
              << ((algoMask & 128) ? " ToeplitzFactorAware"  : "")
              << ((algoMask & 256) ? " ButterflyFactorAware" : "")
              << ((algoMask & 512) ? " DenseFactorAware"     : "")
              << std::endl;
    if (algoMask & 64) {
        std::cout << "LIFs requested r: " << keff << std::endl;
        std::cout << "LIFs block b:     " << lifsBlock << std::endl;
        if (forcedLifsBlock == 0) {
            std::cout << "LIFs selection:   automatic from -c" << std::endl;
            std::cout << "LIFs target:      " << std::setprecision(12)
                      << lifsConfidence << std::setprecision(6) << std::endl;
        } else {
            std::cout << "LIFs selection:   forced by -b" << std::endl;
        }
    }
    std::cout << std::string(70, '-') << std::endl;

    FrobeniusLarge<PolyRing>                FT(R);
    FrobeniusLargeSearch<PolyRing>          FTS(R);
    FrobeniusLargeButterfly<PolyRing>       FB(R);
    FrobeniusLargeButterflySearch<PolyRing> FBS(R);
    FrobeniusLargeDense<PolyRing>           FD(R);
    FrobeniusLargeDenseSearch<PolyRing>     FDS(R);
    FrobeniusLargeFactorAware<PolyRing>          FTF(R);
    FrobeniusLargeButterflyFactorAware<PolyRing> FBF(R);
    FrobeniusLargeDenseFactorAware<PolyRing>     FDF(R);
    std::vector<Polynomial> denseResult;
    std::vector<Polynomial> lifsResult;
    std::vector<Polynomial> profileResult;
    std::string profileSource;

    // Select one result for the optional matrix profile. Prefer Dense, then
    // LIFs, then the first requested FrobeniusLarge variant. When profile
    // output is disabled, no extra result is retained or copied for it.
    int profileBit = 0;
    if (printProfile) {
        if      (algoMask & 16) { profileBit = 16; profileSource = "Dense"; }
        else if (algoMask & 64) { profileBit = 64; profileSource = "LIFs"; }
        else if (algoMask &  1) { profileBit =  1; profileSource = "Toeplitz"; }
        else if (algoMask &  2) { profileBit =  2; profileSource = "ToeplitzSearch"; }
        else if (algoMask &  4) { profileBit =  4; profileSource = "Butterfly"; }
        else if (algoMask &  8) { profileBit =  8; profileSource = "ButterflySearch"; }
        else if (algoMask & 32) { profileBit = 32; profileSource = "DenseSearch"; }
        else if (algoMask & 128) { profileBit = 128; profileSource = "ToeplitzFactorAware"; }
        else if (algoMask & 256) { profileBit = 256; profileSource = "ButterflyFactorAware"; }
        else if (algoMask & 512) { profileBit = 512; profileSource = "DenseFactorAware"; }
    }

    if (algoMask &  1) {
        if (profileBit == 1)
            profileResult = bench(FT, M, keff, nruns, R, "Toeplitz");
        else
            bench(FT, M, keff, nruns, R, "Toeplitz");
    }
    if (algoMask &  2) {
        if (profileBit == 2)
            profileResult = bench(FTS, M, keff, nruns, R, "ToeplitzSearch");
        else
            bench(FTS, M, keff, nruns, R, "ToeplitzSearch");
    }
    if (algoMask &  4) {
        if (profileBit == 4)
            profileResult = bench(FB, M, keff, nruns, R, "Butterfly");
        else
            bench(FB, M, keff, nruns, R, "Butterfly");
    }
    if (algoMask &  8) {
        if (profileBit == 8)
            profileResult = bench(FBS, M, keff, nruns, R, "ButterflySearch");
        else
            bench(FBS, M, keff, nruns, R, "ButterflySearch");
    }
    if (algoMask & 16)
        denseResult = bench(FD, M, keff, nruns, R, "Dense");
    if (algoMask & 32) {
        if (profileBit == 32)
            profileResult = bench(FDS, M, keff, nruns, R, "DenseSearch");
        else
            bench(FDS, M, keff, nruns, R, "DenseSearch");
    }
    if (algoMask & 64)
        lifsResult = benchLIFs<InvariantFactors<Field, PolyRing>,
                               SparseMat, PolyRing>(
            IFD, M, keff, lifsBlock, nruns, "LIFs");
    if (algoMask & 128) {
        if (profileBit == 128)
            profileResult = bench(FTF, M, keff, nruns, R,
                                  "ToeplitzFactorAware");
        else
            bench(FTF, M, keff, nruns, R, "ToeplitzFactorAware");
    }
    if (algoMask & 256) {
        if (profileBit == 256)
            profileResult = bench(FBF, M, keff, nruns, R,
                                  "ButterflyFactorAware");
        else
            bench(FBF, M, keff, nruns, R, "ButterflyFactorAware");
    }
    if (algoMask & 512) {
        if (profileBit == 512)
            profileResult = bench(FDF, M, keff, nruns, R,
                                  "DenseFactorAware");
        else
            bench(FDF, M, keff, nruns, R, "DenseFactorAware");
    }

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

    // Profile formatting and optional polynomial output happen only after all
    // timed solver calls (and the optional correctness check) have completed.
    if (printProfile) {
        if (profileBit == 16)
            printInvariantProfile(denseResult, keff, R, profileSource,
                                  printPolynomials);
        else if (profileBit == 64)
            printInvariantProfile(lifsResult, keff, R, profileSource,
                                  printPolynomials);
        else
            printInvariantProfile(profileResult, keff, R, profileSource,
                                  printPolynomials);
    }

    return 0;
}

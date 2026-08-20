/* linbox/tests/test-frobenius-mixed-butterfly.C
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 * Written by Omesh Dhar Dwivedi <odd23@drexel.edu>
 *
 * Random mixed-irreducible Frobenius test suite.  This compares the
 * binary, exponential, and factor-aware searches using only the
 * Butterfly preconditioner.
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 * ========LICENCE========
 */

#include "linbox/linbox-config.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <NTL/ZZ.h>
#include <NTL/lzz_pXFactoring.h>

#include "linbox/util/commentator.h"
#include "linbox/ring/ntl.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/algorithms/frobenius-large-generic.h"
#include "givaro/givtimer.h"

using namespace LinBox;

namespace {

enum SearchMask {
    RunBinary      = 1,
    RunExponential = 2,
    RunFactorAware = 4
};

struct Layer {
    size_t irreducible;
    size_t survival;
};

struct TestResult {
    std::string algorithm;
    std::string profile;
    size_t n;
    size_t nnz;
    double usertime;
    bool pass;
};

uint64_t splitmix64(uint64_t x)
{
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

void resetNTLSeed(uint64_t seed)
{
    const unsigned long s = static_cast<unsigned long>(seed | UINT64_C(1));
    NTL::SetSeed(NTL::to_ZZ(s));
}

int parseSearchMask(const std::string &s)
{
    if (s.empty()) return RunBinary | RunExponential | RunFactorAware;

    int mask = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '2') mask |= RunBinary;
        if (s[i] == '3') mask |= RunExponential;
        if (s[i] == '8') mask |= RunFactorAware;
    }
    return mask;
}

const char *modeName(int mode)
{
    if (mode == 1) return "early";
    if (mode == 2) return "late";
    return "uniform";
}

size_t randomSurvival(std::mt19937_64 &rng, size_t maximum, int mode)
{
    std::uniform_int_distribution<size_t> distribution(1, maximum);
    const size_t a = distribution(rng);
    if (mode == 1) return std::min(a, distribution(rng));
    if (mode == 2) return std::max(a, distribution(rng));
    return a;
}

template<class PolyRing>
void makeDistinctIrreducibles(
    std::vector<typename PolyRing::Element> &irreducibles,
    std::vector<size_t> &degrees,
    const PolyRing &R,
    size_t count,
    size_t maximumDegree,
    uint64_t seed)
{
    typedef typename PolyRing::Element Polynomial;

    irreducibles.clear();
    degrees.clear();
    resetNTLSeed(seed);

    for (size_t j = 0; j < count; ++j) {
        // p_0 is linear, which makes every requested matrix dimension
        // attainable exactly.  The remaining degrees cycle through the
        // requested range, while the irreducibles themselves are random.
        const size_t degree = (j == 0) ? 1 : 1 + (j % maximumDegree);
        Polynomial base, candidate;
        NTL::BuildIrred(base, static_cast<long>(degree));

        bool distinct = false;
        for (size_t attempt = 0; attempt < 10000 && !distinct; ++attempt) {
            NTL::BuildRandomIrred(candidate, base);
            distinct = R.deg(candidate) == degree && R.isIrreducible(candidate);
            for (size_t i = 0; i < irreducibles.size() && distinct; ++i)
                if (R.areEqual(candidate, irreducibles[i])) distinct = false;
        }
        if (!distinct)
            throw std::runtime_error("could not generate enough distinct irreducibles");

        irreducibles.push_back(candidate);
        degrees.push_back(degree);
    }
}

void makeRandomLayers(
    std::vector<Layer> &layers,
    const std::vector<size_t> &degrees,
    size_t matrixDimension,
    size_t factorCount,
    int mode,
    uint64_t seed)
{
    layers.clear();

    std::mt19937_64 rng(seed);

    // A degree-one layer surviving through every position guarantees
    // exactly factorCount nonunit invariant factors.
    layers.push_back(Layer{0, factorCount});
    size_t used = factorCount;

    // Give every requested irreducible its own random lifetime.  Reserve
    // enough room for one layer of each irreducible still to come.
    for (size_t j = 1; j < degrees.size(); ++j) {
        size_t reserve = 0;
        for (size_t z = j + 1; z < degrees.size(); ++z)
            reserve += degrees[z];
        if (used + reserve + degrees[j] > matrixDimension)
            throw std::runtime_error("target dimension is too small for the requested profile");

        const size_t budget = matrixDimension - used - reserve;
        const size_t maximum = std::min(factorCount, budget / degrees[j]);
        const size_t survival = randomSurvival(rng, maximum, mode);
        layers.push_back(Layer{j, survival});
        used += degrees[j] * survival;
    }

    if (used > matrixDimension)
        throw std::runtime_error("target dimension is too small for the requested profile");

    size_t remaining = matrixDimension - used;

    while (remaining > 0) {
        std::vector<size_t> feasible;
        for (size_t j = 0; j < degrees.size(); ++j)
            if (degrees[j] <= remaining) feasible.push_back(j);

        std::uniform_int_distribution<size_t> choose(0, feasible.size() - 1);
        const size_t j = feasible[choose(rng)];
        const size_t maximum = std::min(factorCount, remaining / degrees[j]);
        const size_t survival = randomSurvival(rng, maximum, mode);

        layers.push_back(Layer{j, survival});
        remaining -= degrees[j] * survival;
    }
}

template<class PolyRing>
void buildInvariantFactors(
    std::vector<typename PolyRing::Element> &factors,
    const PolyRing &R,
    const std::vector<typename PolyRing::Element> &irreducibles,
    const std::vector<Layer> &layers,
    size_t factorCount,
    size_t matrixDimension)
{
    typedef typename PolyRing::Element Polynomial;

    factors.resize(factorCount);
    for (size_t i = 0; i < factorCount; ++i)
        R.assign(factors[i], R.one);

    for (size_t a = 0; a < layers.size(); ++a)
        for (size_t i = 0; i < layers[a].survival; ++i)
            R.mulin(factors[i], irreducibles[layers[a].irreducible]);

    size_t degreeSum = 0;
    for (size_t i = 0; i < factors.size(); ++i) {
        if (R.isOne(factors[i]))
            throw std::runtime_error("generated a unit among the requested invariant factors");
        degreeSum += R.deg(factors[i]);
        if (i + 1 < factors.size()) {
            Polynomial remainder;
            R.rem(remainder, factors[i], factors[i + 1]);
            if (!R.isZero(remainder))
                throw std::runtime_error("generated polynomials are not a divisibility chain");
        }
    }
    if (degreeSum != matrixDimension)
        throw std::runtime_error("generated invariant-factor degrees do not sum to n");
}

template<class Field, class PolyRing>
void writeCompanionBlock(
    SparseMatrix<Field, SparseMatrixFormat::CSR> &M,
    const Field &F,
    const PolyRing &R,
    const typename PolyRing::Element &polynomial,
    size_t offset)
{
    typedef typename Field::Element Element;
    typedef typename PolyRing::Coeff Coeff;

    const size_t degree = R.deg(polynomial);
    for (size_t i = 0; i + 1 < degree; ++i)
        M.setEntry(offset + i + 1, offset + i, F.one);

    for (size_t i = 0; i < degree; ++i) {
        Coeff coefficient;
        Element negative;
        R.getCoeff(coefficient, polynomial, i);
        F.neg(negative, coefficient);
        if (!F.isZero(negative))
            M.setEntry(offset + i, offset + degree - 1, negative);
    }
}

template<class Field, class PolyRing>
void buildCompanionSum(
    SparseMatrix<Field, SparseMatrixFormat::CSR> &M,
    const Field &F,
    const PolyRing &R,
    const std::vector<typename PolyRing::Element> &factors,
    size_t matrixDimension)
{
    M.resize(matrixDimension, matrixDimension);
    size_t offset = 0;
    for (size_t i = 0; i < factors.size(); ++i) {
        writeCompanionBlock(M, F, R, factors[i], offset);
        offset += R.deg(factors[i]);
    }
    if (offset != matrixDimension)
        throw std::runtime_error("companion-block dimensions do not sum to n");
    M.finalize();
}

template<class Field>
void sprayMatrix(
    SparseMatrix<Field, SparseMatrixFormat::CSR> &M,
    const Field &F,
    uint64_t seed)
{
    typedef typename Field::Element Element;
    if (M.rowdim() < 2) return;

    std::mt19937_64 rng(seed);
    typename Field::RandIter fieldRandom(F, seed);
    const size_t n = M.rowdim();
    const size_t sprays = static_cast<size_t>(std::ceil(std::sqrt(static_cast<double>(n))));
    std::uniform_int_distribution<size_t> index(0, n - 1);

    for (size_t spray = 0; spray < sprays; ++spray) {
        const size_t i = index(rng);
        size_t j;
        do { j = index(rng); } while (j == i);

        Element alpha;
        do { fieldRandom.random(alpha); } while (F.isZero(alpha));

        std::vector<std::pair<size_t, Element> > rowJ;
        for (size_t column = 0; column < n; ++column) {
            Element value; F.init(value, 0);
            M.getEntry(value, j, column);
            if (!F.isZero(value)) rowJ.push_back(std::make_pair(column, value));
        }
        for (size_t z = 0; z < rowJ.size(); ++z) {
            Element current; F.init(current, 0);
            M.getEntry(current, i, rowJ[z].first);
            Element contribution; F.mul(contribution, alpha, rowJ[z].second);
            F.addin(current, contribution);
            M.setEntry(i, rowJ[z].first, current);
        }

        std::vector<std::pair<size_t, Element> > columnI;
        for (size_t row = 0; row < n; ++row) {
            Element value; F.init(value, 0);
            M.getEntry(value, row, i);
            if (!F.isZero(value)) columnI.push_back(std::make_pair(row, value));
        }
        Element negativeAlpha; F.neg(negativeAlpha, alpha);
        for (size_t z = 0; z < columnI.size(); ++z) {
            Element current; F.init(current, 0);
            M.getEntry(current, columnI[z].first, j);
            Element contribution; F.mul(contribution, negativeAlpha, columnI[z].second);
            F.addin(current, contribution);
            M.setEntry(columnI[z].first, j, current);
        }
    }
    M.finalize();
}

template<class FrobeniusObject, class Field, class PolyRing>
TestResult runOne(
    FrobeniusObject &algorithm,
    const std::string &algorithmName,
    const std::string &profileName,
    const Field &F,
    const PolyRing &R,
    const SparseMatrix<Field, SparseMatrixFormat::CSR> &M,
    const std::vector<typename PolyRing::Element> &expected,
    size_t limit,
    uint64_t algorithmSeed)
{
    typedef typename PolyRing::Element Polynomial;
    std::vector<Polynomial> computed;

    resetNTLSeed(algorithmSeed);
    Givaro::Timer timer;
    timer.clear(); timer.start();
    algorithm.frobeniusInvariants(computed, M, limit);
    timer.stop();

    const size_t expectedCount = limit > 0
        ? std::min(limit, expected.size()) : expected.size();
    bool pass = computed.size() >= expectedCount;
    for (size_t i = 0; i < expectedCount && pass; ++i)
        pass = R.areEqual(computed[i], expected[i]);

    // Once the full expected nonunit prefix was requested, the solver may
    // return a unit suffix.  Accept units, but reject extra nonunit factors.
    if (expectedCount == expected.size())
        for (size_t i = expectedCount; i < computed.size() && pass; ++i)
            pass = R.isOne(computed[i]);

    if (!pass) {
        std::cout << "\n[FAIL] " << algorithmName << " on " << profileName << "\n";
        std::cout << "  computed " << computed.size()
                  << ", expected nonunit prefix " << expectedCount << "\n";
        const size_t shown = std::max(computed.size(), expectedCount);
        for (size_t i = 0; i < shown; ++i) {
            std::cout << "  [" << i << "] computed=";
            if (i < computed.size()) R.write(std::cout, computed[i]);
            else std::cout << "<missing>";
            std::cout << " expected=";
            if (i < expectedCount) R.write(std::cout, expected[i]);
            else std::cout << "<unit suffix>";
            std::cout << "\n";
        }
    }

    TestResult result;
    result.algorithm = algorithmName;
    result.profile = profileName;
    result.n = M.rowdim();
    result.nnz = M.size();
    result.usertime = timer.usertime();
    result.pass = pass;
    return result;
}

void printTable(const std::vector<TestResult> &results)
{
    std::cout << "\n" << std::left
              << std::setw(27) << "Algorithm"
              << std::setw(22) << "Profile"
              << std::setw(8)  << "n"
              << std::setw(10) << "nnz"
              << std::setw(12) << "Time (s)"
              << std::setw(8)  << "Pass" << "\n";
    std::cout << std::string(87, '-') << "\n";
    for (size_t i = 0; i < results.size(); ++i)
        std::cout << std::left
                  << std::setw(27) << results[i].algorithm
                  << std::setw(22) << results[i].profile
                  << std::setw(8)  << results[i].n
                  << std::setw(10) << results[i].nnz
                  << std::setw(12) << std::fixed << std::setprecision(4)
                  << results[i].usertime
                  << std::setw(8) << (results[i].pass ? "PASS" : "FAIL")
                  << "\n";
    std::cout << "\n";
}

template<class PolyRing>
void printProfileSummary(
    const PolyRing &R,
    const std::vector<typename PolyRing::Element> &irreducibles,
    const std::vector<size_t> &degrees,
    const std::vector<Layer> &layers,
    const std::vector<typename PolyRing::Element> &factors,
    size_t trial,
    size_t matrixDimension,
    int mode,
    bool verbose)
{
    std::cout << "Profile " << trial << ": n=" << matrixDimension
              << " nonunit-factors=" << factors.size()
              << " irreducibles=" << irreducibles.size()
              << " layers=" << layers.size()
              << " mode=" << modeName(mode)
              << " deg(F1)=" << R.deg(factors[0]) << " degrees=[";
    for (size_t j = 0; j < degrees.size(); ++j) {
        if (j) std::cout << ',';
        std::cout << degrees[j];
    }
    std::cout << "]\n";

    if (!verbose) return;
    for (size_t j = 0; j < irreducibles.size(); ++j) {
        std::cout << "  p" << j << '=';
        R.write(std::cout, irreducibles[j]);
        std::cout << " survival=[";
        bool first = true;
        for (size_t a = 0; a < layers.size(); ++a)
            if (layers[a].irreducible == j) {
                if (!first) std::cout << ',';
                std::cout << layers[a].survival;
                first = false;
            }
        std::cout << "]\n";
    }
}

} // anonymous namespace

int main(int argc, char **argv)
{
    int characteristic = 10000019;
    int matrixDimension = 1000;
    int factorCount = 20;
    int irreducibleCount = 6;
    int maximumDegree = 3;
    int profileCount = 1;
    int mode = 0;
    int limit = 0;
    int seed = static_cast<int>(std::time(NULL));
    bool verbose = false;
    std::string algorithms = "238";

    static Argument args[] = {
        { 'n', "-n N", "Exact matrix dimension", TYPE_INT, &matrixDimension },
        { 'm', "-m M", "Number of nonunit invariant factors", TYPE_INT, &factorCount },
        { 'i', "-i I", "Number of distinct irreducibles", TYPE_INT, &irreducibleCount },
        { 'd', "-d D", "Maximum irreducible degree", TYPE_INT, &maximumDegree },
        { 't', "-t T", "Number of random profiles", TYPE_INT, &profileCount },
        { 'g', "-g G", "Layer survival: 0=uniform, 1=early drops, 2=late drops", TYPE_INT, &mode },
        { 'k', "-k K", "Invariant-factor limit (0 = all)", TYPE_INT, &limit },
        { 'p', "-p P", "Prime field characteristic", TYPE_INT, &characteristic },
        { 'r', "-r R", "Random seed", TYPE_INT, &seed },
        { 'a', "-a A", "Butterfly searches: 2=binary 3=exponential 8=factor-aware", TYPE_STR, &algorithms },
        { 'v', "-v", "Print irreducibles and every layer survival", TYPE_BOOL, &verbose },
        END_OF_ARGUMENTS
    };

    parseArguments(argc, argv, args);

    const int searchMask = parseSearchMask(algorithms);
    if (characteristic <= 2 || matrixDimension <= 0 || factorCount <= 0 ||
        irreducibleCount <= 0 || maximumDegree <= 0 || profileCount <= 0 ||
        limit < 0 || mode < 0 || mode > 2 || searchMask == 0) {
        std::cerr << "invalid test parameters\n";
        return -1;
    }

    typedef NTL_zz_p Field;
    typedef NTL_zz_pX PolyRing;
    typedef PolyRing::Element Polynomial;
    typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;

    Field F(characteristic);
    PolyRing R(F);

    FrobeniusLargeButterfly<PolyRing> binary(R);
    FrobeniusLargeButterflySearch<PolyRing> exponential(R);
    FrobeniusLargeButterflyFactorAware<PolyRing> factorAware(R);

    std::cout << "=== Mixed-Irreducible Frobenius Suite (Butterfly) ===\n"
              << "Field: GF(" << characteristic << ") seed=" << seed
              << " k=" << limit << " profiles=" << profileCount << "\n"
              << "Searches: "
              << ((searchMask & RunBinary) ? "Binary " : "")
              << ((searchMask & RunExponential) ? "Exponential " : "")
              << ((searchMask & RunFactorAware) ? "FactorAware " : "")
              << "\n";

    std::vector<TestResult> results;
    bool allPass = true;

    try {
        for (int trial = 0; trial < profileCount; ++trial) {
            const uint64_t trialSeed = splitmix64(static_cast<uint64_t>(
                static_cast<uint32_t>(seed)) + static_cast<uint64_t>(trial));

            std::vector<Polynomial> irreducibles;
            std::vector<size_t> degrees;
            std::vector<Layer> layers;
            std::vector<Polynomial> expected;

            makeDistinctIrreducibles(irreducibles, degrees, R,
                static_cast<size_t>(irreducibleCount),
                static_cast<size_t>(maximumDegree), splitmix64(trialSeed));
            makeRandomLayers(layers, degrees,
                static_cast<size_t>(matrixDimension),
                static_cast<size_t>(factorCount), mode,
                splitmix64(trialSeed + 1));
            buildInvariantFactors(expected, R, irreducibles, layers,
                static_cast<size_t>(factorCount),
                static_cast<size_t>(matrixDimension));

            SparseMat matrix(F,
                static_cast<size_t>(matrixDimension),
                static_cast<size_t>(matrixDimension));
            buildCompanionSum(matrix, F, R, expected,
                static_cast<size_t>(matrixDimension));
            sprayMatrix(matrix, F, splitmix64(trialSeed + 2));

            printProfileSummary(R, irreducibles, degrees, layers, expected,
                static_cast<size_t>(trial), static_cast<size_t>(matrixDimension),
                mode, verbose);

            const std::string profileName = std::string("mixed_")
                + modeName(mode) + "_" + std::to_string(trial);
            const uint64_t algorithmSeed = splitmix64(trialSeed + 3);

            if (searchMask & RunBinary) {
                TestResult r = runOne(binary, "ButterflyBinary", profileName,
                    F, R, matrix, expected, static_cast<size_t>(limit), algorithmSeed);
                allPass = allPass && r.pass; results.push_back(r);
            }
            if (searchMask & RunExponential) {
                TestResult r = runOne(exponential, "ButterflyExponential", profileName,
                    F, R, matrix, expected, static_cast<size_t>(limit), algorithmSeed);
                allPass = allPass && r.pass; results.push_back(r);
            }
            if (searchMask & RunFactorAware) {
                TestResult r = runOne(factorAware, "ButterflyFactorAware", profileName,
                    F, R, matrix, expected, static_cast<size_t>(limit), algorithmSeed);
                allPass = allPass && r.pass; results.push_back(r);
            }
        }
    }
    catch (const std::exception &error) {
        std::cerr << "generation error: " << error.what() << "\n";
        return -1;
    }

    printTable(results);
    return allPass ? 0 : -1;
}
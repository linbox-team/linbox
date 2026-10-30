/* linbox/tests/test-frobenius-suite.C
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 * Written by Omesh Dhar Dwivedi <odd23@drexel.edu>
 *
 * Constructed-matrix Frobenius correctness and timing suite.
 *
 * Timing policy:
 *   - Every timed invocation reports wall, user, and system time.
 *   - The final table treats wall time as the primary observed runtime.
 *   - There are no warmup invocations.
 *
 * LIFs block policy:
 *   - For k > 0, -c chooses the smallest probability-certified block b.
 *   - -b forces b and overrides -c.
 *   - For k = 0 (all factors), LIFs uses b=n separately for each case.
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 * ========LICENCE========
 */

#include "linbox/linbox-config.h"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "linbox/ring/modular.h"
#include "linbox/util/commentator.h"
#include "linbox/ring/ntl.h"

#include "linbox/algorithms/frobenius-large-generic.h"
#include "linbox/algorithms/invariant-factors.h"

#include "test-frobenius-suite.h"

using namespace LinBox;

// bit0=Toeplitz, bit1=ToeplitzSearch, bit2=Butterfly,
// bit3=ButterflySearch, bit4=Dense, bit5=DenseSearch, bit6=LIFs,
// bit7=ToeplitzFactorAware, bit8=ButterflyFactorAware,
// bit9=DenseFactorAware.
int parseAlgoMask(const std::string &selection)
{
    if (selection.empty()) return 1023;

    int mask = 0;
    for (size_t i = 0; i < selection.size(); ++i)
        mask |= 1 << (selection[i] - '0');
    return mask;
}

bool validAlgoSelection(const std::string &selection)
{
    if (selection.empty()) return true;
    for (size_t i = 0; i < selection.size(); ++i) {
        if (selection[i] < '0' || selection[i] > '9') return false;
    }
    return true;
}

int main(int argc, char **argv)
{
    int         p = 10000019;
    int         e = 1;
    int         kInput = 0;
    int         seed = static_cast<int>(std::time(NULL));
    int         sInput = 0;
    int         wInput = 0;
    int         hInput = 0;
    int         poly = 0;
    int         nruns = 1;
    double      lifsConfidence = 0.999999;
    int         forcedLifsBlockInput = 0;
    bool        allowRisky = false;
    std::string algoStr;

    static Argument args[] = {
        { 'k', "-k K", "Number of invariant factors (0 = all)",
            TYPE_INT, &kInput },
        { 'p', "-p P", "Characteristic of field GF(p)",
            TYPE_INT, &p },
        { 'e', "-e E", "Extension degree (NTL_zz_p requires 1)",
            TYPE_INT, &e },
        { 'r', "-r R", "Random seed",
            TYPE_INT, &seed },
        { 'n', "-n N", "Number of timed runs (default: 1)",
            TYPE_INT, &nruns },
        { 'c', "-c C", "LIFs target success probability per call (default: 0.999999)",
            TYPE_DOUBLE, &lifsConfidence },
        { 'b', "-b B", "LIFs block size (0 = choose automatically from -c)",
            TYPE_INT, &forcedLifsBlockInput },
        { 's', "-s S", "Custom: number of distinct block sizes",
            TYPE_INT, &sInput },
        { 'w', "-w W", "Custom: repetitions per block size",
            TYPE_INT, &wInput },
        { 'H', "-H H", "Custom: step between block sizes",
            TYPE_INT, &hInput },
        { 'q', "-q Q", "Custom polynomial (0=x, 1=x-1, 2=x+1)",
            TYPE_INT, &poly },
        { 'a', "-a A", "Algorithms: 0=Toeplitz 1=ToeplitzSearch 2=Butterfly 3=ButterflySearch 4=Dense 5=DenseSearch 6=LIFs 7=ToeplitzFactorAware 8=ButterflyFactorAware 9=DenseFactorAware (default: all)",
            TYPE_STR, &algoStr },
        { 'x', "-x", "Acknowledge risky settings for an enormous custom case",
            TYPE_BOOL, &allowRisky },
        END_OF_ARGUMENTS
    };

    parseArguments(argc, argv, args);

    if (p <= 1) {
        std::cerr << "Error: field characteristic must be greater than 1."
                  << std::endl;
        return 1;
    }
    if (e != 1) {
        std::cerr << "Error: this suite uses NTL_zz_p, which supports only "
                  << "extension degree -e 1." << std::endl;
        return 1;
    }
    if (kInput < 0) {
        std::cerr << "Error: invariant-factor count cannot be negative."
                  << std::endl;
        return 1;
    }
    if (nruns <= 0) {
        std::cerr << "Error: number of timed runs must be positive."
                  << std::endl;
        return 1;
    }
    if (forcedLifsBlockInput < 0) {
        std::cerr << "Error: LIFs block size cannot be negative."
                  << std::endl;
        return 1;
    }
    if (sInput < 0 || wInput < 0 || hInput < 0) {
        std::cerr << "Error: custom dimensions cannot be negative."
                  << std::endl;
        return 1;
    }
    if (poly < 0 || poly > 2) {
        std::cerr << "Error: custom polynomial -q must be 0, 1, or 2."
                  << std::endl;
        return 1;
    }
    if (!validAlgoSelection(algoStr)) {
        std::cerr << "Error: -a must contain only algorithm digits 0 through 9."
                  << std::endl;
        return 1;
    }

    const bool anyCustomDimension =
        sInput != 0 || wInput != 0 || hInput != 0;
    const bool completeCustomDimensions =
        sInput > 0 && wInput > 0 && hInput > 0;
    if (anyCustomDimension && !completeCustomDimensions) {
        std::cerr << "Error: a custom case requires positive -s, -w, and -H "
                  << "together." << std::endl;
        return 1;
    }

    const size_t k = static_cast<size_t>(kInput);
    const size_t customS = static_cast<size_t>(sInput);
    const size_t customW = static_cast<size_t>(wInput);
    const size_t customH = static_cast<size_t>(hInput);
    const size_t forcedLifsBlock =
        static_cast<size_t>(forcedLifsBlockInput);
    const int algoMask = parseAlgoMask(algoStr);

    // The built-in cases are tiny and intentionally run all algorithms. Apply
    // multi-day safeguards only to an enormous user-constructed case.
    size_t customN = 0;
    const size_t enormousOrder = 100000;
    if (completeCustomDimensions) {
        if (!checkedCaseOrder(customS, customW, customH, customN)) {
            std::cerr << "Error: custom dimensions overflow matrix order."
                      << std::endl;
            return 1;
        }

        const bool enormousCustom = customN >= enormousOrder;
        const bool riskyAllAlgorithms = enormousCustom && algoStr.empty();
        const bool riskyRepetitions = enormousCustom && nruns > 1;
        const bool riskyAllFactors = enormousCustom && k == 0;
        if ((riskyAllAlgorithms || riskyRepetitions || riskyAllFactors) &&
            !allowRisky) {
            std::cerr << "Safety stop: custom case n=" << customN
                      << " is enormous (threshold n>=" << enormousOrder
                      << ")." << std::endl;
            if (riskyAllAlgorithms)
                std::cerr << "  -a was omitted, so all ten algorithms would run."
                          << std::endl;
            if (riskyRepetitions)
                std::cerr << "  -n " << nruns
                          << " requests multiple timed runs." << std::endl;
            if (riskyAllFactors)
                std::cerr << "  -k 0 requests all invariant factors."
                          << std::endl;
            std::cerr << "Specify safer settings or add -x to acknowledge the "
                      << "workload." << std::endl;
            return 1;
        }
        if ((riskyAllAlgorithms || riskyRepetitions || riskyAllFactors) &&
            allowRisky)
            std::cerr << "Safety warning: -x accepted risky settings for custom "
                      << "case n=" << customN << "." << std::endl;
    }

    typedef NTL_zz_p  Field;
    typedef NTL_zz_pX PolyRing;

    Field    F(p, e);
    PolyRing R(F);

    size_t lifsBlock = 0;
    if (algoMask & 64) {
        if (k == 0 && forcedLifsBlock != 0) {
            std::cerr << "Error: -b cannot be forced with -k 0; computing all "
                      << "factors requires LIFs block b=n for each case."
                      << std::endl;
            return 1;
        }

        if (k > 0 && forcedLifsBlock == 0) {
            if (!(lifsConfidence > 0.0 && lifsConfidence < 1.0)) {
                std::cerr << "Error: LIFs target probability must satisfy "
                          << "0 < -c < 1." << std::endl;
                return 1;
            }

            // Do not call InvariantFactors::min_block_size(): its direct
            // Givaro::Integer-to-size_t conversion is ambiguous for NTL_zz_p
            // on some platforms. This is the same formula with the safe
            // cardinality output-parameter overload.
            size_t q = 0;
            F.cardinality(q);
            if (q < 2) {
                std::cerr << "Error: field cardinality must be at least 2."
                          << std::endl;
                return 1;
            }

            const double probabilityConstant = (q == 2) ? 3.0 : 2.0;
            const double rawOversampling =
                std::log(probabilityConstant /
                         (1.0 - std::sqrt(lifsConfidence))) /
                std::log(static_cast<double>(q));
            if (!std::isfinite(rawOversampling) || rawOversampling < 0.0 ||
                rawOversampling >
                    static_cast<double>(std::numeric_limits<size_t>::max())) {
                std::cerr << "Error: could not represent the selected LIFs "
                          << "oversampling size." << std::endl;
                return 1;
            }

            const size_t oversampling = static_cast<size_t>(
                std::ceil(rawOversampling));
            if (k > std::numeric_limits<size_t>::max() - oversampling) {
                std::cerr << "Error: selected LIFs block size overflows size_t."
                          << std::endl;
                return 1;
            }
            lifsBlock = k + oversampling;
        } else if (k > 0) {
            lifsBlock = forcedLifsBlock;
        }

        if (k > 0 && lifsBlock < k) {
            std::cerr << "Error: LIFs block size b=" << lifsBlock
                      << " is smaller than requested k=" << k << "."
                      << std::endl;
            return 1;
        }
    }

    FrobeniusLarge<PolyRing>                    FT(R);
    FrobeniusLargeSearch<PolyRing>              FTS(R);
    FrobeniusLargeButterfly<PolyRing>           FB(R);
    FrobeniusLargeButterflySearch<PolyRing>     FBS(R);
    FrobeniusLargeDense<PolyRing>               FD(R);
    FrobeniusLargeDenseSearch<PolyRing>         FDS(R);
    InvariantFactors<Field, PolyRing>            IFD(F, R);
    FrobeniusLargeFactorAware<PolyRing>          FTF(R);
    FrobeniusLargeButterflyFactorAware<PolyRing> FBF(R);
    FrobeniusLargeDenseFactorAware<PolyRing>     FDF(R);

    std::srand(seed);

    std::cout << "=== Frobenius Test Suite ===" << std::endl;
    std::cout << "Field:     GF(" << p << ")" << std::endl;
    std::cout << "seed:      " << seed << std::endl;
    std::cout << "k:         ";
    if (k == 0) std::cout << "all";
    else        std::cout << k;
    std::cout << std::endl;
    std::cout << "runs:      " << nruns << std::endl;
    std::cout << "warmups:   0 (disabled)" << std::endl;
    std::cout << "timing:    wall/real elapsed time is primary; user and system "
              << "CPU are also recorded per run and in aggregate" << std::endl;
    if (completeCustomDimensions)
        std::cout << "custom n:  " << customN << std::endl;
    std::cout << "Algorithms:"
              << ((algoMask &   1) ? " Toeplitz"                 : "")
              << ((algoMask &   2) ? " ToeplitzSearch"           : "")
              << ((algoMask &   4) ? " Butterfly"                : "")
              << ((algoMask &   8) ? " ButterflySearch"          : "")
              << ((algoMask &  16) ? " Dense"                    : "")
              << ((algoMask &  32) ? " DenseSearch"              : "")
              << ((algoMask &  64) ? " LIFs"                     : "")
              << ((algoMask & 128) ? " ToeplitzFactorAware"      : "")
              << ((algoMask & 256) ? " ButterflyFactorAware"     : "")
              << ((algoMask & 512) ? " DenseFactorAware"         : "")
              << std::endl;

    if (algoMask & 64) {
        std::cout << "LIFs requested r: ";
        if (k == 0) std::cout << "all";
        else        std::cout << k;
        std::cout << std::endl;
        if (k == 0) {
            std::cout << "LIFs block b:     n (selected per case)" << std::endl;
            std::cout << "LIFs selection:   full-order for -k 0" << std::endl;
        } else {
            std::cout << "LIFs block b:     " << lifsBlock << std::endl;
            if (forcedLifsBlock == 0) {
                std::cout << "LIFs selection:   automatic from -c" << std::endl;
                std::cout << "LIFs target/call: " << std::setprecision(12)
                          << lifsConfidence << std::setprecision(6)
                          << std::endl;
            } else {
                std::cout << "LIFs selection:   forced by -b" << std::endl;
            }
        }
    }

    const bool pass = runSuite(FT, FTS, FB, FBS, FD, FDS, IFD,
                               FTF, FBF, FDF, F, R, k, algoMask, seed,
                               nruns, lifsBlock, customS, customW,
                               customH, poly);
    return pass ? 0 : 1;
}

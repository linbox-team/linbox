/*
 * Copyright(C) LinBox
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#pragma once

#include <linbox/algorithms/rns.h>
#include <linbox/solutions/methods.h>

namespace LinBox {
    /**
     * The algorithm find out the p-adic writing of A^{-1} * b.
     * So that A^{-1} * b = c0 + c1 * p + c2 * p^2 + ... + c{k-1} * p^{k-1}.
     * The chosen p is multi-modular.
     *
     * It is based on Chen/Storjohann RNS-based p-adic lifting.
     * Based on https://cs.uwaterloo.ca/~astorjoh/p92-chen.pdf
     * A BLAS Based C Library for Exact Linear Algebra on Integer Matrices (ISSAC 2009)
     * But it has been slightly modified in order to use BLAS3 multiplication within the main loop.
     *
     *  RNS Dixon algorithm goes this way:
     *      (i)     Use (p1, ..., pl) primes with an arbitrary l.
     *      (ii)    Algorithm goes:
     *                  for i = 1 .. l:
     *                  |   Bi = A^{-1} mod pi                      < Pre-computing
     *                  [r1|...|rl] = [b|...|b]
     *                  [y1|...|yl] = [0|...|0]
     *                  for j = 1 .. k:
     *                  |   for i = 1 .. l:
     *                  |   |   (Qi, Ri) = such that ri = pi Qi + Ri with |Ri| < pi
     *                  |   |   ci = Bi ri mod pi                   < Matrix-vector in Z/pZ
     *                  |   |   yi = yi + ci * pi^(i-1)             < Done over ZZ
     *                  |   V = [R1|...|Rl] - A [c1|...|cl]         < Matrix-matrix in ZZ
     *                  |   for i = 1 .. l:
     *                  |   |   ri = Qi + (Vi / pi)
     *              @note The computation of V can be done in a RNS system such that each RNS
     * base-prime is bigger than each (p1, ..., pl). This way, [R1|...|Rl] and [c1|...|cl] are
     * zero-cost to get in the RNS system.
     *      (iii)   y = CRT_Reconstruct(y1, ..., yl)
     *      (iv)    x = Rational_Reconstruct(y)
     *
     * One can configure how many primes are used with `Method::DixonRNS.primeBaseLength`.
     * According to the paper, a value of lp = 2 (ln(n) + log2(||A||)) or without the factor 2
     * can be used, but it depends on the problem, really.
     */
    template <class _Field, class _Ring, class _PrimeGenerator>
    class MultiModLiftingContainer final : public LiftingContainer<_Ring> {
        using BaseClass = LiftingContainer<_Ring>;

    public:
        using Ring = _Ring;
        using Field = _Field;
        using PrimeGenerator = _PrimeGenerator;

        using RNSSystem = FFPACK::rns_double;
        using RNSDomain = FFPACK::RNSInteger<FFPACK::rns_double>;
        using RNSElement = typename RNSDomain::Element;
        using RNSElementPtr = typename RNSDomain::Element_ptr;

        using IElement = typename Ring::Element;
        using IMatrix = DenseMatrix<_Ring>;
        using IVector = DenseVector<_Ring>;
        using FElement = typename Field::Element;
        using FMatrix = DenseMatrix<_Field>;
        using FVector = DenseVector<_Field>;

    public:
        // -------------------
        // ----- Main behavior

        // @fixme Split to inline file
        MultiModLiftingContainer(const Ring& ring, PrimeGenerator primeGenerator, const IMatrix& A,
                                 const IVector& b, const Method::DixonRNS& m)
            : _ring(ring)
            , _method(m)
            , _A(A)
            , _b(b)
            , _n(A.rowdim())
            , _rMatrix(_ring)
            , _qMatrix(_ring)
        {
            linbox_check(A.rowdim() == A.coldim());

            // std::cout << "----------" << std::endl;
            // A.write(std::cout << "A: ", Tag::FileFormat::Maple) << std::endl;
            // std::cout << "b: " << b << std::endl;

            // This will contain the primes or our MultiMod basis
            _primesCount = m.primesCount;
            _primes.resize(_primesCount);

            // Some preparation work
            Integer infinityNormA;
            InfinityNorm(infinityNormA, A);
            double logInfinityNormA = Givaro::logtwo(infinityNormA);

            {
                // Based on Chen-Storjohann's paper, this is the bit size
                // of the needed RNS basis for the residue computation
                double rnsBasisBitSize = std::ceil(1.0 + Givaro::logtwo(1 + infinityNormA * _n));
                _rnsPrimesCount = std::ceil(rnsBasisBitSize / (primeGenerator.getBits() - 1));
                _rnsPrimes.resize(_rnsPrimesCount);
                // std::cout << "_rnsPrimesCount: " << _rnsPrimesCount << std::endl;

                std::vector<double> primes;
                for (auto j = 0u; j < _primesCount + _rnsPrimesCount; ++j) {
                    auto p = *primeGenerator;
                    ++primeGenerator;

                    // @note std::lower_bound finds the iterator where to put p in the sorted
                    // container. The name of the routine might be strange, but, hey, that's not my
                    // fault.
                    auto lb = std::lower_bound(primes.begin(), primes.end(), p);
                    if (lb != primes.end() && *lb == p) {
                        --j;
                        continue;
                    }

                    // Inserting the primes at the right place to keep the array sorted
                    primes.insert(lb, p);
                }

                // We take the smallest primes for our MultiMod basis
                std::copy(primes.begin(), primes.begin() + _primesCount, _primes.begin());

                // And the others for our RNS basis
                std::copy(primes.begin() + _primesCount, primes.end(), _rnsPrimes.begin());

                // We check that we really need all the primes within the RNS basis,
                // as the first count was just an upper estimation.
                double bitSize = 0.0;
                for (int h = _rnsPrimes.size() - 1; h >= 0; --h) {
                    bitSize += Givaro::logtwo(_rnsPrimes[h]);

                    if (bitSize > rnsBasisBitSize && h > 0) {
                        _rnsPrimes.erase(_rnsPrimes.begin(), _rnsPrimes.begin() + h);
                        _rnsPrimesCount -= h;
                        break;
                    }
                }
            }

            // Setting fields up
            for (auto& pj : _primes) {
                _fields.emplace_back(pj);
            }

            // Initialize all inverses
            // @note An inverse mod some p within DixonSolver<Dense> was already computed,
            // and pass through to the lifting container. Here, we could use that, but we have
            // to keep control of generated primes, so that the RNS base has bigger primes
            // than the .
            // commentator().start("[MMLifting][Init] A^{-1} mod pj precomputations");
            {
                _B.reserve(_primesCount);
                for (auto& F : _fields) {
                    _B.emplace_back(A, F);
                }

// @fixme To be replaced with Paladin
#pragma omp parallel for
                for (auto j = 0u; j < _primesCount; ++j) {
                    int nullity = 0;
                    auto& F = _fields[j];
                    BlasMatrixDomain<Field> bmd(F);
                    bmd.invin(_B[j], nullity);
                    if (nullity > 0) {
                        // @fixme Should redraw another prime!
                        throw LinBoxError("Wrong prime, sorry.");
                    }
                }
            }
            // commentator().stop("[MMLifting][Init] A^{-1} mod pj precomputations");

            // Making A into the RNS domain
            {
                _rnsSystem = new RNSSystem(_rnsPrimes);
                _rnsDomain = new RNSDomain(*_rnsSystem);
                _rnsA = FFLAS::fflas_new(*_rnsDomain, _n, _n);
                _rnsc = FFLAS::fflas_new(*_rnsDomain, _n, _primesCount);
                _rnsR = FFLAS::fflas_new(*_rnsDomain, _n, _primesCount);
                _rnsPrimesInverses = FFLAS::fflas_new(*_rnsDomain, _primesCount);

                // @note So that 2^(16*cmax) is the max element of A.
                double cmax = logInfinityNormA / 16.;
                FFLAS::finit_rns(*_rnsDomain, _n, _n, std::ceil(cmax), A.getPointer(), A.stride(),
                                 _rnsA);
            }

            // Compute the inverses of pj for each RNS prime
            {
                for (auto j = 0u; j < _primesCount; ++j) {
                    auto prime = _primes[j];

                    auto& rnsPrimeInverse = _rnsPrimesInverses[j];
                    auto stride = rnsPrimeInverse._stride;

                    for (auto h = 0u; h < _rnsPrimesCount; ++h) {
                        auto& rnsF = _rnsSystem->_field_rns[h];
                        auto& primeInverse = rnsPrimeInverse._ptr[h * stride];
                        rnsF.inv(primeInverse, prime);
                    }
                }
            }

            // Compute how many iterations are needed
            {
                double log2PrimesProduct = 0.0;
                for (auto& pj : _primes) {
                    log2PrimesProduct += Givaro::logtwo(pj);
                }

                auto hb = RationalSolveHadamardBound(A, b);
                _log2Bound = hb.solutionLogBound;
                _numBound = hb.numBound;
                _denBound = hb.denBound;

                // _iterationsCount = log2(2 * N * D) / log2(p1 * p2 * ...)
                _iterationsCount = std::ceil(_log2Bound / log2PrimesProduct);
                std::cout << "_iterationsCount " << _iterationsCount << std::endl;
            }

            //----- Locals setup

            _rMatrix = IMatrix(_ring, _n, _primesCount);
            _qMatrix = IMatrix(_ring, _n, _primesCount);

            _FR.reserve(_primesCount);
            for (auto j = 0u; j < _primesCount; ++j) {
                auto& F = _fields[j];

                _FR.emplace_back(F, _n);

                // Initialize all residues to b
                for (auto i = 0u; i < _n; ++i) {
                    _rMatrix.refEntry(i, j) = _b[i];
                }
            }
        }

        ~MultiModLiftingContainer()
        {
            FFLAS::fflas_delete(_rnsR);
            FFLAS::fflas_delete(_rnsc);
            FFLAS::fflas_delete(_rnsA);
            delete _rnsDomain;
            delete _rnsSystem;
        }

        // --------------------------
        // ----- LiftingContainer API

        const Ring& ring() const final { return _ring; }

        /// The length of the container.
        size_t length() const final { return _iterationsCount; }

        /// The dimension of the problem/solution.
        size_t size() const final { return _n; }

        /// @note Useless, but in the API.
        const IElement& prime() const final { return _ring.one; }

        // ------------------------------
        // ----- NOT LiftingContainer API
        // ----- but still needed

        double log2Bound() const { return _log2Bound; }
        Integer numBound() const { return _numBound; }
        Integer denBound() const { return _denBound; }

        uint32_t primesCount() const { return _primesCount; }
        const FElement& prime(uint32_t index) const { return _primes.at(index); }
        const std::vector<Field>& primesFields() const { return _fields; }

        // --------------
        // ----- Iterator

        /**
         * Returns false if the next digit cannot be computed (bad modulus).
         * c is a vector of integers but all element are below p = p1 * ... * pl
         */
        bool next(std::vector<FVector>& digits)
        {
            VectorDomain<Ring> IVD(_ring);
            BlasMatrixDomain<Ring> IMD(_ring);

// commentator().start("[MultiModLifting] c = A^{-1} r mod p");
#pragma omp parallel for
            for (auto j = 0u; j < _primesCount; ++j) {
                auto pj = _primes[j];
                auto& FR = _FR[j];
                uint64_t upj = pj;

                // @note There is no VectorDomain::divmod yet.
                // Euclidian division so that rj = pj Qj + Rj
                uint64_t uR;
                for (auto i = 0u; i < _n; ++i) {
                    Integer::divmod(_qMatrix.refEntry(i, j), uR, _rMatrix.getEntry(i, j), upj);
                    // @note No need to init, because we know that uR < pj,
                    // so that would do an extra check.
                    FR[i] = static_cast<FElement>(uR);
                }

                // digit = A^{-1} * R mod pj
                auto& digit = digits[j];
                auto& B = _B[j];
                B.apply(digit, FR);

                // Store the very same result in an RNS system,
                // but fact is all the primes of the RNS system are bigger
                // than the modulus used to compute the digit, we just copy the result for
                // everybody.
                for (auto i = 0u; i < _n; ++i) {
                    setRNSMatrixElementAllResidues(_rnsR, _primesCount, i, j, FR[i]);
                    setRNSMatrixElementAllResidues(_rnsc, _primesCount, i, j, digit[i]);
                }
            }
            // commentator().stop("[MultiModLifting] c = A^{-1} r mod p");

            // ----- Compute the next residues!

            // r <= Q + (R - A c) / p

            using RNSParallel = FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::RNSModulus,
                                                              FFLAS::StrategyParameter::Threads>;
            using FGEMMParallel = FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,
                                                                FFLAS::StrategyParameter::Threads>;

            // @fixme @cpernet @jgdumas Should we move that PAR_BLOCK outside of the function
            // and let the user do it?
            // commentator().start("[MultiModLifting] FGEMM R <= R - Ac");
            // Firstly compute R <= R - A c as a fgemm within the RNS domain.
            if (_method.rnsFgemmType == RnsFgemmType::BothSequential) {
                rns_fgemm<FFLAS::ParSeqHelper::Sequential, FFLAS::ParSeqHelper::Sequential>();
            }
            else if (_method.rnsFgemmType == RnsFgemmType::BothParallel) {
                rns_fgemm<RNSParallel, FGEMMParallel>();
            }
            else if (_method.rnsFgemmType == RnsFgemmType::ParallelFgemmOnly) {
                rns_fgemm<FFLAS::ParSeqHelper::Sequential, FGEMMParallel>();
            }
            else if (_method.rnsFgemmType == RnsFgemmType::ParallelRnsOnly) {
                rns_fgemm<RNSParallel, FFLAS::ParSeqHelper::Sequential>();
            }

            // commentator().stop("[MultiModLifting] FGEMM R <= R - Ac");

            // We divide each residues by the according pj, which is done by multiplying.
            // @note The matrix _rnsR is RNS-major, meaning that it is stored
            // as [R mod q0][R mod q1][...] where [R mod qh] represents a full matrix.
            // We use this fact to keep better cache coherency.
            // commentator().start("[MultiModLifting] MUL FOR INV R <= R / p");
            auto rnsStride = 0u;
            for (auto h = 0u; h < _rnsPrimesCount; ++h) {
                auto& rnsF = _rnsSystem->_field_rns[h];

#pragma omp parallel for
                for (auto j = 0u; j < _primesCount; ++j) {
                    auto& rnsPrimeInverse = _rnsPrimesInverses[j];
                    auto stridePrimeInverse = rnsPrimeInverse._stride;
                    auto rnsPrimeInverseForRnsPrimeH = rnsPrimeInverse._ptr[h * stridePrimeInverse];

                    for (auto i = 0u; i < _n; ++i) {
                        rnsF.mulin(_rnsR._ptr[rnsStride + (i * _primesCount + j)],
                                   rnsPrimeInverseForRnsPrimeH);
                    }
                }

                rnsStride += _rnsR._stride;
            }
            // commentator().stop("[MultiModLifting] MUL FOR INV R <= R / p");

            // commentator().start("[MultiModLifting] CONVERT TO INTEGER r <= Q + R");
            FFLAS::fconvert_rns(*_rnsDomain, _n, _primesCount, 0, _rMatrix.getWritePointer(),
                                _primesCount, _rnsR + 0);
            IMD.addin(_rMatrix, _qMatrix);
            // commentator().stop("[MultiModLifting] CONVERT TO INTEGER r <= Q + R");

            return true;
        }

    private:
        // Helper function, setting all residues of a matrix element to the very same value.
        // This doesn't check the moduli.
        inline void setRNSMatrixElementAllResidues(RNSElementPtr& A, size_t lda, size_t i, size_t j,
                                                   double value)
        {
            auto& Aij = A[i * lda + j];
            auto stride = Aij._stride;
            for (auto h = 0u; h < _rnsPrimesCount; ++h) {
                Aij._ptr[h * stride] = value;
            }
        }

        // @note This allows us to factor out some of the rns fgemm variants common code.
        template <class RnsParSeq, class FgemmParSeq>
        inline void rns_fgemm()
        {
            PAR_BLOCK
            {
                using ComposedParSeqHelper = FFLAS::ParSeqHelper::Compose<RnsParSeq, FgemmParSeq>;
                using MMHelper =
                    FFLAS::MMHelper<RNSDomain, FFLAS::MMHelperAlgo::Classic,
                                    FFLAS::ModeCategories::DefaultTag, ComposedParSeqHelper>;
                ComposedParSeqHelper composedParSeqHelper(NUM_THREADS, NUM_THREADS);
                MMHelper mmHelper(*_rnsDomain, -1, composedParSeqHelper);

                FFLAS::fgemm(*_rnsDomain, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, _n,
                             _primesCount, _n, _rnsDomain->mOne, _rnsA, _n, _rnsc, _primesCount,
                             _rnsDomain->one, _rnsR, _primesCount, mmHelper);
            }
        }

    public: // @fixme BACK TO PRIVATE!
        const Ring& _ring;
        Method::DixonRNS _method; // A copy of the user-provided method.

        // The problem: A^{-1} * b
        const IMatrix& _A;
        const IVector& _b;

        double _log2Bound;
        Integer _numBound;
        Integer _denBound;

        RNSSystem* _rnsSystem = nullptr;
        RNSDomain* _rnsDomain = nullptr;
        RNSElementPtr _rnsA; // The matrix A, but in the RNS system
        // A matrix of digits c[j], being the current digits mod pj, in the RNS system
        RNSElementPtr _rnsc;
        RNSElementPtr _rnsR;
        size_t _rnsPrimesCount = 0u;
        // Stores the inverse of pj within the RNS base prime into _rnsPrimesInverses[j]
        RNSElementPtr _rnsPrimesInverses;

        std::vector<double> _primes;
        std::vector<double> _rnsPrimes;
        // Length of the ci sequence. So that p^{k-1} > 2ND (Hadamard bound).
        size_t _iterationsCount = 0u;
        size_t _n = 0u;           // Row/column dimension of A.
        size_t _primesCount = 0u; // How many primes. Equal to _primes.size().

        std::vector<FMatrix> _B;    // Inverses of A mod p[i]
        std::vector<Field> _fields; // All fields Modular<p[i]>

        //----- Iteration

        std::vector<FVector> _FR;
        IMatrix _rMatrix;
        IMatrix _qMatrix;
    };
}

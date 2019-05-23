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
            , _A(A)
            , _b(b)
            , _n(A.rowdim())
        {
            linbox_check(A.rowdim() == A.coldim());

            A.write(std::cout << "A: ", Tag::FileFormat::Maple) << std::endl;
            std::cout << "b: " << b << std::endl;

            // @fixme Pass it through Method::DixonRNS (and rename it Method::DixonMultiMod?)
            _l = 2;
            std::cout << "l: " << _l << std::endl;

            // Generating primes
            IElement iTmp;
            _ring.assign(_p, _ring.one);
            for (auto j = 0u; j < _l; ++j) {
                // @fixme Ensure that all primes are different
                // @fixme Take into account bestBitSize!
                _primes.emplace_back(*primeGenerator);
                _fields.emplace_back(_primes.back());
                _ring.init(iTmp, _primes.back());
                _ring.mulin(_p, iTmp);

                std::cout << "primes[" << j << "]: " << Integer(_primes.back()) << std::endl;

                ++primeGenerator;
            }

            std::cout << "p: " << _p << std::endl;

            // Compute how many iterations are needed
            auto hb = RationalSolveHadamardBound(A, b);
            double pLog = Givaro::logtwo(_p);
            // _k = log2(2 * N * D) / log2(p)
            _k = std::ceil((1.0 + hb.numLogBound + hb.denLogBound) / pLog);
            std::cout << "k: " << _k << std::endl;

            // @fixme Fact is RationalReconstruction which needs numbound and denbound
            // expects them to be in non-log...
            _ring.init(_numbound, Integer(1) << static_cast<uint64_t>(std::ceil(hb.numLogBound)));
            _ring.init(_denbound, Integer(1) << static_cast<uint64_t>(std::ceil(hb.denLogBound)));

            // Initialize all inverses
            // @note An inverse mod some p within DixonSolver<Dense> was already computed,
            // and pass through to the lifting container. Here, we could use that, but we have
            // to keep control of generated primes, so that the RNS base has bigger primes
            // than the .
            {
                _B.reserve(_l);

                for (const auto& F : _fields) {
                    BlasMatrixDomain<Field> bmd(F);
                    _B.emplace_back(F, _n, _n);
                    auto& Bpi = _B.back();

                    // @fixme Taken for rational-solver.inl. BETTER USE REBIND!!!
                    for (size_t i = 0; i < _n; ++i) {
                        for (size_t j = 0; j < _n; ++j) {
                            F.init(Bpi.refEntry(i, j), A.getEntry(i, j));
                        }
                    }

                    // @fixme @cpernet Use FFLAS directly, so that we can have a REAL in place inv.
                    bmd.invin(Bpi);

                    Bpi.write(std::cout << "B mod " << Integer(F.characteristic()) << ": ",
                              Tag::FileFormat::Maple)
                        << std::endl;
                }
            }
        }

        // --------------------------
        // ----- LiftingContainer API

        const Ring& ring() const final { return _ring; }

        /// The length of the container.
        size_t length() const final { return _k; }

        /// The dimension of the problem/solution.
        size_t size() const final { return _n; }

        /**
         * We are compliant to the interface even though
         * p is multi-modular and thus not a prime per se.
         */
        const IElement& prime() const final { return _p; }

        // ------------------------------
        // ----- NOT LiftingContainer API
        // ----- but still needed

        const IElement numbound() const { return _numbound; }

        const IElement denbound() const { return _denbound; }

        // --------------
        // ----- Iterator

        /**
         * Needed API for rational reconstruction.
         * Each call to next() will update
         */
        class const_iterator {
        private:
            const MultiModLiftingContainer& _lc;
            std::vector<IVector> _r; // @todo Could be a matrix? Might not be useful, as it is never
                                     // used directly in computations.
            std::vector<IVector> _Q;
            std::vector<IVector> _R; // @fixme This one should be expressed in a RNS system q, and
                                     // HAS TO BE A MATRIX for gemm.
            std::vector<FVector>
                _Fc; // @note No need to be a matrix, as we will embed it into an RNS system later.
            size_t _position;

            // @fixme Better use Givaro::RNSSystem?
            RNS<true> _pRns; // RNS system for primes

        public:
            const_iterator(const MultiModLiftingContainer& lc, size_t position = 0)
                : _lc(lc)
                , _position(position)
            {
                VectorDomain<Ring> IVD(_lc._ring);

                _pRns.init(_lc._primes);

                _r.reserve(_lc._l);
                _Q.reserve(_lc._l);
                _R.reserve(_lc._l);
                _Fc.reserve(_lc._l);
                for (auto j = 0u; j < _lc._l; ++j) {
                    auto& F = _lc._fields[j];

                    _r.emplace_back(_lc._ring, _lc._n);
                    _Q.emplace_back(_lc._ring, _lc._n);
                    _R.emplace_back(_lc._ring, _lc._n);
                    _Fc.emplace_back(F, _lc._n);

                    // Initialize all residues to b
                    _r.back() = _lc._b; // Copying data
                }

                // @fixme Allocate c

                // @todo Set up an RNS system
            }

            /**
             * Returns false if the next digit cannot be computed (bad modulus).
             * c is a vector of integers but all element are below p = p1 * ... * pl
             */
            bool next(IVector& c)
            {
                std::cout << "----- NEXT" << std::endl;

                VectorDomain<Ring> IVD(_lc._ring);

                // @fixme Should be done in parallel!
                for (auto j = 0u; j < _lc._l; ++j) {
                    auto pj = _lc._primes[j];
                    auto& r = _r[j];
                    auto& Q = _Q[j];
                    auto& R = _R[j];

                    // @todo @cpernet Is there a VectorDomain::divmod somewhere?
                    // Euclidian division so that rj = pj Qj + Rj
                    for (auto i = 0u; i < _lc._n; ++i) {
                        // @fixme @cpernet Is this OK for any Ring or should we be sure we are using
                        // Integers?
                        _lc._ring.quoRem(Q[i], R[i], r[i], pj);
                    }

                    std::cout << "--- FOR " << Integer(pj) << std::endl;
                    std::cout << "r: " << r << std::endl;
                    std::cout << "Q: " << Q << std::endl;
                    std::cout << "R: " << R << std::endl;

                    // Convert R to the field
                    // @fixme @cpernet Could this step be ignored?
                    // If not, put that in already allocated memory, and not use a temporary here.
                    auto& F = _lc._fields[j];
                    FVector FR(F, R); // rebind

                    auto& B = _lc._B[j];
                    auto& Fc = _Fc[j];
                    B.apply(Fc, FR);

                    std::cout << "Fc: " << Fc << std::endl;

                    // @todo Convert _c[i] to RNS
                }

                // ----- CRT reconstruct c from (cj)

                std::cout << "--- CRT reconstruction" << std::endl;

                // @cpernet Is that RNS system what I should use? I tweaked it so that I can use it.
                std::vector<FElement> fElements(_lc._l);
                for (auto i = 0u; i < _lc._n; ++i) {
                    for (auto j = 0u; j < _lc._l; ++j) {
                        fElements[j] = _Fc[j][i];
                    }
                    // @fixme This cra function should be called reconstruct or such.
                    _pRns.cra(c[i], fElements);
                }

                std::cout << "c: " << c << std::endl;

                // ----- Compute the next residue!

                std::cout << "--- Residue update" << std::endl;

                // @note This is a dummy implementation, for now.

                // r <= (rj - A c) / pj
                for (auto j = 0u; j < _lc._l; ++j) {
                    auto pj = _lc._primes[j];
                    auto& r = _r[j];
                    auto& Q = _Q[j];
                    auto& R = _R[j];

                    auto& Fc = _Fc[j];
                    // @fixme For now, we convert cj to integer,
                    // but it should be converted into a RNS system, on pre-allocated memory.
                    IVector Ic(_lc._ring, Fc);

                    // @fixme Should become a matrix-matrix multiplication!
                    // @fixme Should be able to do a gemv
                    _lc._A.apply(r, Ic); // r = A c
                    IVD.negin(r);        // r = - A c
                    IVD.addin(r, R);     // r = R - A c

                    // r = (R - A c) / pj
                    IElement Ipj;
                    _lc._ring.init(Ipj, pj);
                    for (auto i = 0u; i < _lc._n; ++i) {
                        _lc._ring.divin(r[i], Ipj); // @fixme Is there a divin in VectorDomain?
                    }

                    IVD.addin(r, Q); // r = Q + (R - A c) / pj
                }

                ++_position;
                return true;
            }

            bool operator!=(const const_iterator& iterator) const
            {
                return _position != iterator._position;
            }

            bool operator==(const const_iterator& iterator) const
            {
                return _position == iterator._position;
            }
        };

        const_iterator begin() const { return const_iterator(*this); }
        const_iterator end() const { return const_iterator(*this, _k); }

    private:
        const Ring& _ring;

        // The problem: A^{-1} * b
        const IMatrix& _A;
        const IVector& _b;

        IElement _numbound;
        IElement _denbound;

        IElement _p;                   // The global modulus for lifting: a multiple of all _primes.
        std::vector<FElement> _primes; // @fixme We might want something else as a type!
        size_t _k; // Length of the ci sequence. So that p^{k-1} > 2ND (Hadamard bound).
        size_t _n; // Row/column dimension of A.
        size_t _l; // How many primes. Equal to _primes.size().

        std::vector<FMatrix> _B;    // Inverses of A mod p[i]
        std::vector<Field> _fields; // All fields Modular<p[i]>
    };
}

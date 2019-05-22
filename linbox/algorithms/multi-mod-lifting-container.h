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
            , _n(A.rowdim())
        {
            linbox_check(A.rowdim() == A.coldim());

            A.write(std::cout << "A: ", Tag::FileFormat::Maple) << std::endl;
            std::cout << "b: " << b << std::endl;

            // @fixme Have l = log(||A||) + log(n) or so
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

                std::cout << "primes[" << i << "]: " << Integer(_primes.back()) << std::endl;

                ++primeGenerator;
            }

            std::cout << "p: " << _p << std::endl;

            // Compute how many iterations are needed
            auto hb = RationalSolveHadamardBound(A, b);
            double pLog = Givaro::logtwo(_p);
            _k = std::ceil((1.0 + hb.numLogBound + hb.denLogBound)
                           / pLog); // log2(2 * N * D) / log2(p)
            std::cout << "k: " << _k << std::endl;

            // @fixme Fact is RationalReconstruction which needs numbound and denbound
            // expects them to be in non-log...
            _ring.init(_numbound, Integer(1) << static_cast<uint64_t>(std::ceil(hb.numLogBound)));
            _ring.init(_denbound, Integer(1) << static_cast<uint64_t>(std::ceil(hb.denLogBound)));

            // Initialize all inverses
            // @fixme Somehow, the inverse mod p within DixonSolver<Dense> was already computed,
            // and pass through to the lifting container. Here, we can't do that, because p is
            // bigger than what DixonSolver<Dense> thought about it. So there might be a lot of
            // computation done there that is completely useless when using this container. Meaning
            // that we need a RNSDixonSolver<Dense>.
            {
                for (const auto& F : _fields) {
                    BlasMatrixDomain<Field> bmd(F);
                    auto Bpi = std::make_unique<FMatrix>(F, _n, _n);

                    // @fixme Taken for rational-solver.inl. BETTER USE REBIND!!!
                    for (size_t i = 0; i < _n; ++i) {
                        for (size_t j = 0; j < _n; ++j) {
                            F.init(Bpi->refEntry(i, j), A.getEntry(i, j));
                        }
                    }

                    bmd.invin(*Bpi);
                    Bpi->write(std::cout << "B mod " << Integer(F.characteristic()) << ": ", Tag::FileFormat::Maple) << std::endl;
                    _B.emplace_back(std::move(Bpi));
                }
            }

            // @note As _r is row major, we store each ri on each row.
            // So that r[i] = current residue for p[i].
            // _r.init(_ring, l, b.size());
            // for (auto i = 0u; i < l; ++i) {
            //     // @fixme Is there a vector domain to copy to a matrix?
            //     for (auto j = 0u; j < b.size(); ++j) {
            //         _ring.assign(_r[i][j], b[j]);
            //     }
            // }

            // @fixme Allocate Q and R
            // @fixme Allocate c

            // @todo Set up an RNS system
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
         * p is multi-modular and thus not a prime.
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
            size_t _position;

        public:
            const_iterator(const MultiModLiftingContainer& lc, size_t position = 0)
                : _lc(lc)
                , _position(position)
            {
                // @fixme Initialize reisdue _r
            }

            /**
             * Returns false if the next digit cannot be computed (bad modulus).
             * ci is a vector of integers but all element are below p = p1 * ... * pl
             */
            bool next(IVector& ci)
            {
                /*  for i = 1 .. l:
                 *  |   (Qi, Ri) = such that ri = pi Qi + Ri with |Ri| < pi
                 *  |   ci = Bi Ri mod pi                   < Matrix-vector in Z/pZ
                 *  V = [R1|...|Rl] - A [c1|...|cl]         < Matrix-matrix in ZZ
                 *  for i = 1 .. l:
                 *  |   ri = Qi + (Vi / pi)
                 */

                // @fixme Could be parallel!
                for (auto j = 0u; j < _l; ++j) {
                    // @fixme How to do euclidian division?
                    // ri = pi Qi + Ri

                    // @todo If R might already be a field element
                    // @cpernet!!!
                    // @fixme We will probably need a low-level API
                    // so that we can say that the j-th row of _ci takes
                    // the result of B * R mod pj
                    // _B[j]->apply(*_ci[j], *_R[j]);

                    // @todo Convert _c[i] to RNS
                }

                // @fixme CRT reconstruct ci from (cij)

                std::cout << "ci: " << ci << std::endl;

                // @fixme How can we do A [c1|...|cl] in ZZ if the ci are in the fields?

                // @fixme Compute the next residue!

                // @fixme @note For us, Aci is a matrix!

                // // compute Aci = _matA * ci
                // IVector Aci(_lc.ring(), _lc.size());
                // // @fixme _lc._MAD.applyV(Aci, ci, _res);

                // // update _res -= Aci
                // // @fixme _lc._VDR.subin(_res, Aci);
                // typename BlasVector<Ring>::iterator p0;

                // // update _res = _res / p
                // int index = 0;
                // for (p0 = _res.begin(); p0 != _res.end(); ++p0, ++index) {
                //     _lc.ring().divin(*p0, _lc._p);
                // }

                // increase position of the iterator
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

        IElement _numbound;
        IElement _denbound;

        IElement _p;                   // The global modulus for lifting: a multiple of all _primes.
        std::vector<FElement> _primes; // @fixme We might want something else as a type!
        size_t _k; // Length of the ci sequence. So that p^{k-1} > 2ND (Hadamard bound).
        size_t _n; // Row/column dimension of A.
        size_t _l; // How many primes. Equal to _primes.size().

        // @note r is a big matrix in ZZ holding all residues
        // IMatrix _r;
        FMatrix _ci; // Contains [ci mod p0 | ... | ci mod p{l-1}] on each row.
        std::vector<std::unique_ptr<FMatrix>> _B; // Inverses of A mod p[i]
        // std::vector<IVector> _Q;
        // std::vector<FVector> _R;
        std::vector<Field> _fields;
    };
}

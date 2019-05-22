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

        using IElement = typename _Ring::Element;
        using IMatrix = DenseMatrix<_Ring>;
        using IVector = DenseVector<_Ring>;
        using FMatrix = DenseMatrix<_Field>;
        using FVector = DenseVector<_Field>;

    public:
        // -------------------
        // ----- Main behavior

        // @fixme Split to inline file
        MultiModLiftingContainer(const Ring& ring, PrimeGenerator primeGenerator, const IMatrix& A,
                                 const IVector& b, const Method::DixonRNS& m)
            : _ring(ring)
            , _r(_ring)
            , _c(_field)
        {
            // @fixme Compute hadamard and such

            // @note From baseClass, we have _length = log2(2 * N * D)

            // @fixme Have l = log(||A||) + log(n) or so

            // @fixme Initialize fields _F[i]

            // Ap[0] = A mod p[0]
            // Ap[1] = A mod p[1]

            // B[0] = inv(Ap[0]) mod p[0]
            // B[1] = inv(Ap[1]) mod p[1]
            // @fixme How?

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

        // @fixme USELESS?
        IVector& nextdigit(IVector& digit, const IVector& residu) const
        {
            // @fixme The residu can't be r, here!
            // So the overall does a lot more job than it needs.
            // See below for the solution.

            // @fixme With this design, are we forced to CRT_Reconstruct each ci?
            // Is this bad?
            // If we don't want that, we need to not extent LiftingContainerBase,
            // and reimplement some of the behavior.
            // Because the only thing needed to user API (rational reconstruction)
            // is bool next (IVector& digit) from iterator.

            /*  for i = 1 .. l:
             *  |   (Qi, Ri) = such that ri = pi Qi + Ri with |Ri| < pi
             *  |   ci = Bi Ri mod pi                   < Matrix-vector in Z/pZ
             *  V = [R1|...|Rl] - A [c1|...|cl]         < Matrix-matrix in ZZ
             *  for i = 1 .. l:
             *  |   ri = Qi + (Vi / pi)
             */

            // @fixme Could be parallel!
            // for (auto i = 0u; i < l; ++i) {
            //     Hom<Ring, Field> hom(_ring, _F[i]);

            //     // @fixme How to do euclidian division?
            //     // ri = pi Qi + Ri

            //     // @todo If R might already be a field element
            //     _B[i]->apply(_c[i], hom.convert(_R[i]));

            //     // @todo Convert _c[i] to RNS
            // }

            // @fixme How can we do A [c1|...|cl] in ZZ if the ci are in the fields?

            // @fixme Compute the next residue!

            return digit;
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

		const IElement numbound() const
		{
			return _numbound;
		}

		const IElement denbound() const
		{
			return _denbound;
		}

        // --------------
        // ----- Iterator

        /**
         * Needed API for rational reconstruction.
         * Each call to next() will update
         */
        class const_iterator {
        private:
            BlasVector<Ring> _res;
            const MultiModLiftingContainer& _lc;
            size_t _position;

        public:
            const_iterator(const MultiModLiftingContainer& lc, size_t end = 0)
                : _lc(lc)
                , _position(end)
            {
                // @fixme Initialize _residue
            }

            /**
             * Returns false if the next digit cannot be computed (bad modulus).
             */
            bool next(IVector& digit)
            {
                // compute v2 = _matA * digit
                IVector v2(_lc.ring(), _lc.size());
                // @fixme _lc._MAD.applyV(v2, digit, _res);

                // update _res -= v2
                // @fixme _lc._VDR.subin(_res, v2);
                typename BlasVector<Ring>::iterator p0;

                // update _res = _res / p
                int index = 0;
                for (p0 = _res.begin(); p0 != _res.end(); ++p0, ++index) {
                    _lc.ring().divin(*p0, _lc._p);
                }

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
        Field _field;

        IElement _numbound;
        IElement _denbound;

        IElement _p;
        size_t _k; //< Length of the ci sequence. So that p^{k-1} > 2ND (Hadamard bound)
        size_t _n; //< Row/column dimension of A.

        // @note r is a big matrix in ZZ holding all residues
        IMatrix _r;
        FMatrix _c;
        std::vector<FMatrix> _B; // Inverses of A mod p[i]
        std::vector<IVector> _Q;
        std::vector<FVector> _R;
        std::vector<Field> _F;
    };
}

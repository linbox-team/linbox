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
     * @fixme Should this just be a different LiftingContainer?
     *
     * Chen/Storjohann RNS-based p-adic lifting.
     * The algorithm solves Ax = b over the integers.
     *
     * Based on https://cs.uwaterloo.ca/~astorjoh/p92-chen.pdf
     * A BLAS Based C Library for Exact Linear Algebra on Integer Matrices (ISSAC 2009)
     *
     *  Dixon algorithm goes this way:
     *      (i)     Compute B := A^{-1} mod p
     *              (with p a random number which is hopefully orthogonal to det(A))
     *      (ii)    Compute (ci) such that A^{-1} b = c0 + c1 p + ... + ci p^i mod p^{i+1}
     *              Which means:    r = b
     *                              for i = 0 .. k-1:
     *                              |   ci = B r mod p
     *                              |   r = (r - A ci) / p
     *              (stop when p^k > 2ND given by Hadamard bound)
     *      (iii)   Rational reconstruct with c = c0 + c1 p + ... + ck p^{k-1} (over the integers)
     *
     * The RNS part:
     *      (i)     Use p = p0p1...p{lp-1} with an arbitrary lp and (q0, q1, ..., q{lq-1}) also primes.
     *      (ii)    We now do our computation in a RNS system (p0, ..., p{lp-1}, q0, ..., q{lq-1}):
    *                               r = b
     *                              for i = 0 .. k-1:
     *                              |   for j = 0 .. lq-1:
     *                              |   |   ci[qj] = Bj r mod qj
     *                              |   (Q, R) = such that r = pQ + R with |R| < p
     *                              |   r = Q + (R - A ci) / p      < Matrix-vector multiplication done in RNS domain
     *                              |                                 and final addition over ZZ
     *              @note (R - A ci) / p can be computed in a RNS system.
     *              We know that (R - A ci) is divisible by p,
     *              so its representation is 0 on all lp first terms of the representation, meaning
     *              we just need representation of ci mod (q0, ..., q{lq-1}).
     *              For the division part, we just have to multiply the RNS representation of (R - A ci) by
     *              (1/p) mod (q0, ..., q{lq-1}).
     *              /!\ @fixme The paper does not talk about matrix-matrix multiplication,
     *              but instead about exploiting RNS.
     *      (iii)   We first RNS-reconstruct the solution before rational reconstruction.
     *
     * One can configure how many primes are used with `Method::DixonRNS.primeBaseLength`.
     * According to the paper, a value of lp = 2 (ln(n) + log2(||A||)) or without the factor 2
     * can be used, but it depends on the problem, really.
     */
    template <class Field, class Ring, class PrimeGenerator>
    class DixonRNSSolver {
    public:
        DixonRNSSolver(const Ring& ring, PrimeGenerator primeGenerator);

        /**
         * Dense solving.
         */
        template <class IntVector, class Vector>
        void solve(IntVector& xNum, typename IntVector::Element& xDen, const DenseMatrix<Ring>& A,
                   const Vector& b, const Method::DixonRNS& m);
    };
}

#include "./dixon-rns-solver.inl"
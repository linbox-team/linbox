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
    template <class Field, class Ring, class PrimeGenerator>
    inline DixonRNSSolver<Field, Ring, PrimeGenerator>::DixonRNSSolver(
        const Ring& ring, PrimeGenerator primeGenerator)
    {
    }

    /**
     * Dense solving.
     */
    template <class Field, class Ring, class PrimeGenerator>
    template <class IntVector, class Vector>
    inline void DixonRNSSolver<Field, Ring, PrimeGenerator>::solve(
        IntVector& xNum, typename IntVector::Element& xDen, const DenseMatrix<Ring>& A,
        const Vector& b, const Method::DixonRNS& m)
    {
    }
}
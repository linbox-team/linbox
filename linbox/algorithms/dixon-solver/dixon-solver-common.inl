/*
 * Copyright (C) LinBox Team
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */
#pragma once

namespace LinBox {
    /*! @brief NO DOC !
     * @bug why is this hard coded ?
     */
    template <class Prime>
    inline bool checkBlasPrime(const Prime& p)
    {
        return p < Prime(67108863);
    }

	// Check that all primes within the vector are valid for blas.
    template <>
    inline bool checkBlasPrime(const BlasVector<Givaro::ZRing<Integer>>& p)
    {
        for (size_t i = 0; i < p.size(); ++i) {
            if (!checkBlasPrime(p[i])) {
                return false;
            }
        }

        return true;
    }
}
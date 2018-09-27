/* Copyright (C) 2018 The LinBox group
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

#include <vector>
#include <gmpxx.h>

/**
 * Provides way to pack matrices and vectors of Integer into std::vectors.
 */

namespace LinBox {

    template <class Matrix> void gmp_unpackMat(Matrix& M, int* A_mp_alloc, int* A_a_size, std::vector<mp_limb_t>& A_mp_data)
    {
        size_t ni = M.rowdim(), nj = M.rowdim();

        __mpz_struct* ptr;
        for (size_t i = 0; i < ni; ++i) {
            for (size_t j = 0; j < nj; ++j) {

                ptr = const_cast<__mpz_struct*>(M.getEntry(i, j).get_mpz());
                A_mp_alloc[j + i * nj] = ptr->_mp_alloc;
                A_a_size[j + i * nj] = ptr->_mp_size;
                mp_limb_t* a_array = ptr->_mp_d;
                for (long k = 0; k < ptr->_mp_alloc; ++k) A_mp_data.push_back(a_array[k]);
            }
        }
    }

    template <class Matrix> void gmp_packMat(Matrix& M, int* A_mp_alloc, int* A_a_size, std::vector<mp_limb_t>& A_mp_data)
    {
        size_t ni = M.rowdim(), nj = M.rowdim();
        __mpz_struct* ptr2;
        size_t count = 0;
        Givaro::Integer temp;

        for (size_t i = 0; i < ni; ++i) {
            for (size_t j = 0; j < nj; ++j) {

                ptr2 = const_cast<__mpz_struct*>(temp.get_mpz());
                ptr2->_mp_alloc = A_mp_alloc[j + i * nj];
                ptr2->_mp_size = A_a_size[j + i * nj];
                _mpz_realloc(ptr2, ptr2->_mp_alloc);
                for (long k = 0; k < ptr2->_mp_alloc; ++k) {
                    ptr2->_mp_d[k] = (A_mp_data[k + count]);
                }
                count += ptr2->_mp_alloc;
                M.setEntry(i, j, temp);
            }
        }
    }

    template <class Matrix>
    void gmp_unpackSparseMat(Matrix& M, std::vector<int>& A_mp_alloc, std::vector<int>& A_a_size,
                             std::vector<mp_limb_t>& A_mp_data, std::vector<long>& A_index)
    {
        size_t ni = M.rowdim(), nj = M.rowdim();
        Givaro::ZRing<Givaro::Integer> ZZ;

        __mpz_struct* ptr;
        for (size_t i = 0; i < ni; ++i) {
            for (size_t j = 0; j < nj; ++j) {
                if (!ZZ.areEqual(M.getEntry(i, j), ZZ.zero)) {

                    ptr = const_cast<__mpz_struct*>(M.getEntry(i, j).get_mpz());
                    A_mp_alloc.push_back(ptr->_mp_alloc);
                    A_a_size.push_back(ptr->_mp_size);
                    mp_limb_t* a_array = ptr->_mp_d;
                    A_index.push_back(i);
                    A_index.push_back(j);
                    for (long k = 0; k < ptr->_mp_alloc; ++k) A_mp_data.push_back(a_array[k]);
                }
            }
        }
    }

    template <class Matrix>
    void gmp_packSparseMat(Matrix& M, std::vector<int>& A_mp_alloc, std::vector<int>& A_a_size, std::vector<mp_limb_t>& A_mp_data,
                           std::vector<long>& A_index)
    {
        size_t ni = M.rowdim(), nj = M.rowdim();
        Givaro::ZRing<Givaro::Integer> ZZ;
        __mpz_struct* ptr2;
        size_t count = 0;
        Givaro::Integer temp;

        for (size_t i = 0; i < ni; ++i) {
            for (size_t j = 0; j < nj; ++j) {

                M.setEntry(i, j, ZZ.zero);
            }
        }

        for (size_t i = 0; i < A_a_size.size(); ++i) {

            // ptr = const_cast<__mpz_struct*>(M.getEntry().get_mpz());
            ptr2 = const_cast<__mpz_struct*>(temp.get_mpz());
            ptr2->_mp_alloc = A_mp_alloc[i];
            ptr2->_mp_size = A_a_size[i];
            _mpz_realloc(ptr2, ptr2->_mp_alloc);

            for (long k = 0; k < ptr2->_mp_alloc; ++k) {
                ptr2->_mp_d[k] = (A_mp_data[k + count]);
            }
            count += ptr2->_mp_alloc;

            M.setEntry(A_index[i * 2], A_index[i * 2 + 1], temp);
        }
    }

    template <class Vector> void gmp_unpackVec(Vector& V, int* B_mp_alloc, int* B_a_size, std::vector<mp_limb_t>& B_mp_data)
    {
        size_t nj = V.size();
        Givaro::Integer temp;

        // Split vector B into arrays
        __mpz_struct* ptr;
        for (size_t j = 0; j < nj; j++) {

            ptr = const_cast<__mpz_struct*>(V.getEntry(j).get_mpz());
            B_mp_alloc[j] = ptr->_mp_alloc;
            B_a_size[j] = ptr->_mp_size;
            mp_limb_t* a_array = ptr->_mp_d;
            for (long i = 0; i < ptr->_mp_alloc; ++i) {
                B_mp_data.push_back(a_array[i]);
            }
        }
    }

    template <class Vector> void gmp_packVec(Vector& V, int* B_mp_alloc, int* B_a_size, std::vector<mp_limb_t>& B_mp_data)
    {
        size_t nj = V.size();
        // Reconstruction of vector B
        __mpz_struct* ptr2;
        size_t count = 0;
        Givaro::Integer temp;
        for (size_t j = 0; j < nj; j++) {
            ptr2 = const_cast<__mpz_struct*>(temp.get_mpz());
            ptr2->_mp_alloc = B_mp_alloc[j];
            ptr2->_mp_size = B_a_size[j];
            _mpz_realloc(ptr2, ptr2->_mp_alloc);
            for (long i = 0; i < ptr2->_mp_alloc; ++i) {
                ptr2->_mp_d[i] = (B_mp_data[i + count]);
            }
            count += ptr2->_mp_alloc;
            V.setEntry(j, temp);
        }
    }

}

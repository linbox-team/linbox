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

#include <linbox/config.h>
#include <linbox/integer.h>
#include <linbox/matrix/dense-matrix.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/vector/blas-vector.h>
#include <vector>

/**
 * Provides way to serialize any kind of data,
 * like matrices, vectors and Integer in a platform-independent way.
 *
 * Serialize functions add data to a prexisting vector of bytes,
 * they return the number of bytes written.
 * Unserialize ones read a vector of bytes starting at a specific offset,
 * the number of bytes read.
 *
 * As a convention, all numbers are written little-endian.
 *
 * @todo GMP Integers can be configured with limbs of different sizes (32 or 64 bits),
 * depending on the machine. We do not handle that right now,
 * but storing info about their dimension might be a good idea,
 * to at least emit a warning.
 */

namespace LinBox {
    // Basic serializations

    uint64_t serialize(std::vector<uint8_t>& bytes, float value);
    uint64_t serialize(std::vector<uint8_t>& bytes, double value);

    uint64_t serialize(std::vector<uint8_t>& bytes, int8_t value);
    uint64_t serialize(std::vector<uint8_t>& bytes, uint8_t value);

    uint64_t serialize(std::vector<uint8_t>& bytes, int16_t value);
    uint64_t serialize(std::vector<uint8_t>& bytes, uint16_t value);

    uint64_t serialize(std::vector<uint8_t>& bytes, int32_t value);
    uint64_t serialize(std::vector<uint8_t>& bytes, uint32_t value);

    uint64_t serialize(std::vector<uint8_t>& bytes, int64_t value);
    uint64_t serialize(std::vector<uint8_t>& bytes, uint64_t value);

    // Basic unserializations

    uint64_t unserialize(float& value, const std::vector<uint8_t>& bytes, uint64_t offset = 0u);
    uint64_t unserialize(double& value, const std::vector<uint8_t>& bytes, uint64_t offset = 0u);

    uint64_t unserialize(int8_t& value, const std::vector<uint8_t>& bytes, uint64_t offset = 0u);
    uint64_t unserialize(uint8_t& value, const std::vector<uint8_t>& bytes, uint64_t offset = 0u);

    uint64_t unserialize(int16_t& value, const std::vector<uint8_t>& bytes, uint64_t offset = 0u);
    uint64_t unserialize(uint16_t& value, const std::vector<uint8_t>& bytes, uint64_t offset = 0u);

    uint64_t unserialize(int32_t& value, const std::vector<uint8_t>& bytes, uint64_t offset = 0u);
    uint64_t unserialize(uint32_t& value, const std::vector<uint8_t>& bytes, uint64_t offset = 0u);

    uint64_t unserialize(int64_t& value, const std::vector<uint8_t>& bytes, uint64_t offset = 0u);
    uint64_t unserialize(uint64_t& value, const std::vector<uint8_t>& bytes, uint64_t offset = 0u);

    /**
     * Serializes an Integer with its underlying __mpz_struct.
     * Returns the number of bytes written.
     *
     * Format is (by bytes count):
     *  0-3  _mp_size   Number of mp_limb_t in _mp_d, can be negative to indicate negative number.
     *  4-.. _mp_d      The limbs of length (abs(_mp_size) * sizeof(mp_limb_t))
     */
    uint64_t serialize(std::vector<uint8_t>& bytes, const Integer& integer);

    /**
     * Unserializes an Integer.
     */
    uint64_t unserialize(Integer& integer, const std::vector<uint8_t>& bytes, uint64_t offset = 0u);

    /**
     * Serializes a BlasMatrix.
     *
     * Format is (by bytes count):
     *  0-7  n      Row dimension of matrix
     *  8-15 m      Column dimension of matrix
     *  16-..       Entries of the matrix, (n * m) row-majored
     */
    template <class Field>
    uint64_t serialize(std::vector<uint8_t>& bytes, const BlasMatrix<Field>& M);

    /**
     * Unserializes a BlasMatrix.
     * The matrix will be resized if necessary.
     */
    template <class Field>
    uint64_t unserialize(BlasMatrix<Field>& M, const std::vector<uint8_t>& bytes, uint64_t offset = 0u);

    /**
     * Serializes a SparseMatrix.
     *
     * Format is (by bytes count):
     *  0-7   n     Row dimension of matrix
     *  8-15  m     Column dimension of matrix
     *  23-..       Entries of the matrix, only non-zero, stored as:
     *      0-7  i      Row index
     *      8-15 j      Column index
     *      16-..       Entry value
     *  (8 bytes) End of sparse entries, with a value of 0xFFFFFFFF'FFFFFFFF
     */
    template <class Field>
    uint64_t serialize(std::vector<uint8_t>& bytes, const SparseMatrix<Field>& M);

    /**
     * Unserializes a SparseMatrix.
     * The matrix will be resized if necessary.
     */
    template <class Field>
    uint64_t unserialize(SparseMatrix<Field>& M, const std::vector<uint8_t>& bytes, uint64_t offset = 0u);

    /**
     * Serializes a BlasVector.
     *
     * Format is (by bytes count):
     *  0-7  l      Length of the vector
     *  16-..       Entries of the vector
     */
    template <class Field>
    uint64_t serialize(std::vector<uint8_t>& bytes, const BlasVector<Field>& V);

    /**
     * Unserializes a BlasVector.
     * The vector will be resized if necessary.
     */
    template <class Field>
    uint64_t unserialize(BlasVector<Field>& V, const std::vector<uint8_t>& bytes, uint64_t offset = 0u);
}

#include "serialization.inl"

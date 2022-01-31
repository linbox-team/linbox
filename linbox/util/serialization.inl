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

#include "serialization.h"

namespace LinBox {
    // ----- Basic serializations

    template <class T>
    inline uint64_t serialize_raw(std::vector<uint8_t>& bytes, const T& value)
    {
        auto uValue = reinterpret_cast<const uint8_t*>(&value);
        bytes.insert(bytes.end(), uValue, uValue + sizeof(T));
        return sizeof(T);
    }

    inline uint64_t serialize(std::vector<uint8_t>& bytes, float value) { return serialize_raw(bytes, value); }
    inline uint64_t serialize(std::vector<uint8_t>& bytes, double value) { return serialize_raw(bytes, value); }

    inline uint64_t serialize(std::vector<uint8_t>& bytes, int8_t value) { return serialize_raw(bytes, value); }
    inline uint64_t serialize(std::vector<uint8_t>& bytes, uint8_t value) { return serialize_raw(bytes, value); }

    inline uint64_t serialize(std::vector<uint8_t>& bytes, int16_t value)
    {
#if defined(__LINBOX_HAVE_BIG_ENDIAN)
        value = __builtin_bswap16(value);
#endif
        auto uValue = reinterpret_cast<uint8_t*>(&value);
        bytes.insert(bytes.end(), uValue, uValue + 2u);
        return 2u;
    }
    inline uint64_t serialize(std::vector<uint8_t>& bytes, uint16_t value)
    {
#if defined(__LINBOX_HAVE_BIG_ENDIAN)
        value = __builtin_bswap16(value);
#endif
        auto uValue = reinterpret_cast<uint8_t*>(&value);
        bytes.insert(bytes.end(), uValue, uValue + 2u);
        return 2u;
    }

    inline uint64_t serialize(std::vector<uint8_t>& bytes, int32_t value)
    {
#if defined(__LINBOX_HAVE_BIG_ENDIAN)
        value = __builtin_bswap32(value);
#endif
        auto uValue = reinterpret_cast<uint8_t*>(&value);
        bytes.insert(bytes.end(), uValue, uValue + 4u);
        return 4u;
    }
    inline uint64_t serialize(std::vector<uint8_t>& bytes, uint32_t value)
    {
#if defined(__LINBOX_HAVE_BIG_ENDIAN)
        value = __builtin_bswap32(value);
#endif
        auto uValue = reinterpret_cast<uint8_t*>(&value);
        bytes.insert(bytes.end(), uValue, uValue + 4u);
        return 4u;
    }

    inline uint64_t serialize(std::vector<uint8_t>& bytes, int64_t value)
    {
#if defined(__LINBOX_HAVE_BIG_ENDIAN)
        value = __builtin_bswap64(value);
#endif
        auto uValue = reinterpret_cast<uint8_t*>(&value);
        bytes.insert(bytes.end(), uValue, uValue + 8u);
        return 8u;
    }
    inline uint64_t serialize(std::vector<uint8_t>& bytes, uint64_t value)
    {
#if defined(__LINBOX_HAVE_BIG_ENDIAN)
        value = __builtin_bswap64(value);
#endif
        auto uValue = reinterpret_cast<uint8_t*>(&value);
        bytes.insert(bytes.end(), uValue, uValue + 8u);
        return 8u;
    }

    // ----- Basic unserializations

    template <class T>
    inline uint64_t unserialize_raw(T& value, const std::vector<uint8_t>& bytes, uint64_t offset)
    {
        auto uValue = &bytes.at(offset);
        memcpy(&value,uValue,sizeof(T));
        return sizeof(T);
    }

    inline uint64_t unserialize(float& value, const std::vector<uint8_t>& bytes, uint64_t offset)
    {
        return unserialize_raw(value, bytes, offset);
    }
    inline uint64_t unserialize(double& value, const std::vector<uint8_t>& bytes, uint64_t offset)
    {
        return unserialize_raw(value, bytes, offset);
    }

    inline uint64_t unserialize(int8_t& value, const std::vector<uint8_t>& bytes, uint64_t offset)
    {
        return unserialize_raw(value, bytes, offset);
    }
    inline uint64_t unserialize(uint8_t& value, const std::vector<uint8_t>& bytes, uint64_t offset)
    {
        return unserialize_raw(value, bytes, offset);
    }

    inline uint64_t unserialize(int16_t& value, const std::vector<uint8_t>& bytes, uint64_t offset)
    {
        auto uValue = &bytes.at(offset);
        memcpy(&value,uValue,sizeof(int16_t));
#if defined(__LINBOX_HAVE_BIG_ENDIAN)
        value = __builtin_bswap16(value);
#endif
        return 2u;
    }
    inline uint64_t unserialize(uint16_t& value, const std::vector<uint8_t>& bytes, uint64_t offset)
    {
        auto uValue = &bytes.at(offset);
        memcpy(&value,uValue,sizeof(uint16_t));
#if defined(__LINBOX_HAVE_BIG_ENDIAN)
        value = __builtin_bswap16(value);
#endif
        return 2u;
    }

    inline uint64_t unserialize(int32_t& value, const std::vector<uint8_t>& bytes, uint64_t offset)
    {
        auto uValue = &bytes.at(offset);
        memcpy(&value,uValue,sizeof(int32_t));
#if defined(__LINBOX_HAVE_BIG_ENDIAN)
        value = __builtin_bswap32(value);
#endif
        return 4u;
    }
    inline uint64_t unserialize(uint32_t& value, const std::vector<uint8_t>& bytes, uint64_t offset)
    {
        auto uValue = &bytes.at(offset);
        memcpy(&value,uValue,sizeof(uint32_t));
#if defined(__LINBOX_HAVE_BIG_ENDIAN)
        value = __builtin_bswap32(value);
#endif
        return 4u;
    }

    inline uint64_t unserialize(int64_t& value, const std::vector<uint8_t>& bytes, uint64_t offset)
    {
        auto uValue = &bytes.at(offset);
        memcpy(&value,uValue,sizeof(int64_t));
#if defined(__LINBOX_HAVE_BIG_ENDIAN)
        value = __builtin_bswap64(value);
#endif
        return 8u;
    }
    inline uint64_t unserialize(uint64_t& value, const std::vector<uint8_t>& bytes, uint64_t offset)
    {
        auto uValue = &bytes.at(offset);
        memcpy(&value,uValue,sizeof(uint64_t));
#if defined(__LINBOX_HAVE_BIG_ENDIAN)
        value = __builtin_bswap64(value);
#endif
        return 8u;
    }

    // ----- Integer

    inline uint64_t serialize(std::vector<uint8_t>& bytes, const Integer& integer)
    {
        const __mpz_struct* mpzStruct = integer.get_mpz();

        // We copy mpSize to a fixed-size variable
        // to ensure compatibility between machines.
        int32_t mpSize = mpzStruct->_mp_size;
        auto bytesWritten = serialize(bytes, mpSize);

        /**
         * @note As said, we're not sure of how many bytes we will write here,
         * because mp_limb_t can be either 32 or 64 bytes-long depending on configuration.
         * We force it to be 64 bit-long so that it becomes platform-independent.
         */
        for (auto i = 0, l = std::abs(mpSize); i < l; ++i) {
            bytesWritten += serialize(bytes, static_cast<uint64_t>(mpzStruct->_mp_d[i]));
        }

        return bytesWritten;
    }

    inline uint64_t unserialize(Integer& integer, const std::vector<uint8_t>& bytes, uint64_t offset)
    {
        __mpz_struct* mpzStruct = integer.get_mpz();

        int32_t mpSize;
        uint64_t bytesRead = 0u;
        bytesRead += unserialize(mpSize, bytes, offset + bytesRead);

        mpzStruct->_mp_alloc = std::abs(mpSize);
        mpzStruct->_mp_size = mpSize;
        _mpz_realloc(mpzStruct, mpzStruct->_mp_alloc);

        // @note We use this proxy limb for the very same reason
        // than above: the GMP real limb can be 64 or 32.
        uint64_t limb;
        for (auto i = 0, l = std::abs(mpSize); i < l; ++i) {
            bytesRead += unserialize(limb, bytes, offset + bytesRead);
            mpzStruct->_mp_d[i] = static_cast<mp_limb_t>(limb);
        }

        return bytesRead;
    }

    // ----- BlasMatrix

    template <class Field>
    inline uint64_t serialize(std::vector<uint8_t>& bytes, const BlasMatrix<Field>& M)
    {
        uint64_t n = M.rowdim(), m = M.coldim();
        auto bytesWritten = serialize(bytes, n);
        bytesWritten += serialize(bytes, m);

        for (uint64_t i = 0; i < n; ++i) {
            for (uint64_t j = 0; j < m; ++j) {
                bytesWritten += serialize(bytes, M.getEntry(i, j));
            }
        }

        return bytesWritten;
    }

    template <class Field>
    inline uint64_t unserialize(BlasMatrix<Field>& M, const std::vector<uint8_t>& bytes, uint64_t offset)
    {
        uint64_t n, m;
        uint64_t bytesRead = 0u;
        bytesRead += unserialize(n, bytes, offset + bytesRead);
        bytesRead += unserialize(m, bytes, offset + bytesRead);

        M.resize(n, m);
        for (uint64_t i = 0; i < n; ++i) {
            for (uint64_t j = 0; j < m; ++j) {
                typename Field::Element entry;
                bytesRead += unserialize(entry, bytes, offset + bytesRead);
                M.setEntry(i, j, entry);
            }
        }

        return bytesRead;
    }


    // ----- SparseMatrix

    template <class Field>
    inline uint64_t serialize(std::vector<uint8_t>& bytes, const SparseMatrix<Field>& M)
    {
        const auto& F = M.field();
        uint64_t n = M.rowdim(), m = M.coldim();
        auto bytesWritten = serialize(bytes, n);
        bytesWritten += serialize(bytes, m);

        for (uint64_t i = 0; i < n; ++i) {
            for (uint64_t j = 0; j < m; ++j) {
                const auto& entry = M.getEntry(i, j);
                if (!F.isZero(entry)) {
                    bytesWritten += serialize(bytes, i);
                    bytesWritten += serialize(bytes, j);
                    bytesWritten += serialize(bytes, entry);
                }
            }
        }

        constexpr const uint64_t endMarker = 0xFFFFFFFFFFFFFFFF;
        bytesWritten += serialize(bytes, endMarker);

        return bytesWritten;
    }

    template <class Field>
    inline uint64_t unserialize(SparseMatrix<Field>& M, const std::vector<uint8_t>& bytes, uint64_t offset)
    {
        uint64_t n, m;
        uint64_t bytesRead = 0u;
        bytesRead += unserialize(n, bytes, offset + bytesRead);
        bytesRead += unserialize(m, bytes, offset + bytesRead);

        M.resize(n, m);
        while (true) {
            uint64_t i, j;
            bytesRead += unserialize(i, bytes, offset + bytesRead);

            // Check if there is the mark of the end of the matrix entries
            if (i == 0xFFFFFFFFFFFFFFFF) {
                break;
            }

            bytesRead += unserialize(j, bytes, offset + bytesRead);

            typename Field::Element entry;
            bytesRead += unserialize(entry, bytes, offset + bytesRead);
            M.setEntry(i, j, entry);
        }

        return bytesRead;
    }


    // ----- BlasVector

    template <class Field>
    inline uint64_t serialize(std::vector<uint8_t>& bytes, const BlasVector<Field>& V)
    {
        uint64_t l = V.size();
        auto bytesWritten = serialize(bytes, l);

        for (uint64_t i = 0; i < l; ++i) {
            bytesWritten += serialize(bytes, V[i]);
        }

        return bytesWritten;
    }

    template <class Field>
    inline uint64_t unserialize(BlasVector<Field>& V, const std::vector<uint8_t>& bytes, uint64_t offset)
    {
        uint64_t l;
        uint64_t bytesRead = 0u;
        bytesRead += unserialize(l, bytes, offset + bytesRead);

        V.resize(l);
        for (uint64_t i = 0; i < l; ++i) {
            bytesRead += unserialize(V[i], bytes, offset + bytesRead);
        }

        return bytesRead;
    }
}

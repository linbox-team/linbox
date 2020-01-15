/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@imag.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#pragma once

#include "common_define.hpp"

#include <cstring> // memset

// Note: No AVX2 on HPAC so no add on integers (_mm256_add_epi64)
// Instead, we rely on SSE4.1
#if not defined(__USE_SIMD)
#define __USE_SIMD 1
#endif

#if __USE_SIMD > 0
    #define __SIMD_STEP 2u
    #define __SIMD_VEC __m128i
    #define __SIMD_ADD_EPI64(A, B) _mm_add_epi64(A, B)
    #define __SIMD_SUB_EPI64(A, B) _mm_sub_epi64(A, B)
#else
    #define __SIMD_STEP 1u
    #define __SIMD_VEC Element_t
    #define __SIMD_ADD_EPI64(A, B) (A + B)
    #define __SIMD_SUB_EPI64(A, B) (A - B)
#endif

// Note: We know that figue has 28 cores
#ifndef DLP_TASKS
#define DLP_TASKS 896
#endif

//----------------//
//----- SPMM -----//

//! Z <- Z + A * Y
//! This will compute modulo p
inline
void spmm_add(const Modulus_t p, Element_t* const Z, const SDataMatrix_t& A, const Element_t* const Y,
              const Index_t cols)
{
    for (auto& sData : A) {
        const Index_t i = sData.rowIndex;
        const Index_t j = sData.colIndex;
        const Element_t data = sData.data;

        for (Index_t k = 0u; k < cols; ++k) {
            Element_t value = __DLP_ABS(data) * Y[j * cols + k];
            value %= p;    
            if (data < 0)
                value = p - value;
                
            Z[i * cols + k] += value;
            if (Z[i * cols + k] >= p)
                Z[i * cols + k] -= p;
        }
    }
}

//! Z <- Z + A * Y
//! This will NOT compute modulo p
void spmm_add(Element_t* const Z, const SMatrixOnes_t& A, const Element_t* const Y,
              const Index_t offset, const Index_t length, const Index_t cols)
{
    auto sZ = reinterpret_cast<__SIMD_VEC*>(Z);
    auto sY = reinterpret_cast<const __SIMD_VEC*>(Y);

    // Note: Here we suppose that cols can perfectly divisible by 2u * __SIMD_STEP
    const Index_t simdCols = cols / __SIMD_STEP;
    
    for (Index_t i = 0u; i < length; ++i) {
        for (const auto c : A.at(offset + i)) {
            for (Index_t k = 0u; k < simdCols; k += 2u) {
                const Index_t si = simdCols * i + k;
                const Index_t sc = simdCols * c + k;

                // This should be done in one cycle
                sZ[si] = __SIMD_ADD_EPI64(sZ[si], sY[sc]);
                sZ[si + 1u] = __SIMD_ADD_EPI64(sZ[si + 1u], sY[sc + 1u]);
            }
        }
    }
}

//! Z <- Z - A * Y
//! This will NOT compute modulo p
inline
void spmm_sub(Element_t* const Z, const SMatrixOnes_t& A, const Element_t* const Y,
              const Index_t offset, const Index_t length, const Index_t cols)
{
    auto sZ = reinterpret_cast<__SIMD_VEC*>(Z);
    auto sY = reinterpret_cast<const __SIMD_VEC*>(Y);

    // Note: Here we suppose that cols can perfectly divisible by 2u * __SIMD_STEP
    const Index_t simdCols = cols / __SIMD_STEP;
    
    for (Index_t i = 0u; i < length; ++i) {
        for (const auto c : A.at(offset + i)) {
            for (Index_t k = 0u; k < simdCols; k += 2u) {
                const Index_t si = simdCols * i + k;
                const Index_t sc = simdCols * c + k;

                // This should be done in one cycle
                sZ[si] = __SIMD_SUB_EPI64(sZ[si], sY[sc]);
                sZ[si + 1u] = __SIMD_SUB_EPI64(sZ[si + 1u], sY[sc + 1u]);
            }
        }
    }
}

inline size_t nnz(const SMatrixOnes_t& A) {
    size_t n(0);
    for(auto const & rows : A)
        n += rows.size();
    return n;
}

//! Z <- A * Y (Y and Z are row-majored)
inline
void spmm(const Modulus_t p, DMatrix_t& Z, const SMatrix_t& A, const DMatrix_t& Y)
{
    const Element_t* const pY = Y.data.data();
    Size_t refa = ((A.nRows / DLP_TASKS) / __SIMD_STEP) * __SIMD_STEP;
    Size_t refb = refa + __SIMD_STEP;
    Size_t rest = (A.nRows - refa*DLP_TASKS) / __SIMD_STEP;
    Size_t LaunchTasks = DLP_TASKS;

//     std::cerr << "refa       : " << refa << std::endl;
//     std::cerr << "refb       : " << refb << std::endl;
//     std::cerr << "rest       : " << rest << std::endl;
//     std::cerr << "LaunchTasks: " << LaunchTasks << std::endl;

    if (refa == 0) {
        LaunchTasks = A.nRows/__SIMD_STEP;
        refa = __SIMD_STEP;
        refb = refa + __SIMD_STEP;
        rest = (A.nRows - refa*LaunchTasks);
        if (rest) ++LaunchTasks;
    }
        
//     std::cerr << "refa       : " << refa << std::endl;
//     std::cerr << "refb       : " << refb << std::endl;
//     std::cerr << "rest       : " << rest << std::endl;
//     std::cerr << "LaunchTasks: " << LaunchTasks << std::endl;
        
    for (Index_t l = 0; l < LaunchTasks; ++l) {
    #pragma omp task shared(Z, A, Y)
    {
        Size_t offset,length;
        if (l<rest) {
            offset = l*refb;
            length = refb;
        } else {
            offset = rest*refb+(l-rest)*refa;
            if (l == LaunchTasks - 1u) {
                length = A.nRows - offset;
            } else {
                length = refa;
            }
        }

        const Size_t realLength = length * Z.nCols;

//         std::cerr << "l     : " << l << std::endl;
//         std::cerr << "length: " << length << std::endl;
//         std::cerr << "offset: " << offset << std::endl;
        
        
        Element_t* const pZ = Z.data.data() + offset * Y.nCols;
    
        // Reset
        memset(pZ, 0, realLength * sizeof(Element_t));
            
        // One
// { std::stringstream buffer; buffer << "spmm_add: (" << Z.nRows << 'x' << Z.nCols << ")=(" << A.ones.size() << '|' << nnz(A.ones) << ")(" << Y.nRows << 'x' << Y.nCols << ')' << offset << '|' << length <<  std::endl; std::cerr << buffer.str(); }
        
        spmm_add(pZ, A.ones, pY, offset, length, Y.nCols);

        // Mones
// { std::stringstream buffer; buffer << "spmm_sub: (" << Z.nRows << 'x' << Z.nCols << ")=(" << A.mOnes.size() << '|' << nnz(A.mOnes) << ")(" << Y.nRows << 'x' << Y.nCols << ')' << offset << '|' << length <<  std::endl; std::cerr << buffer.str(); }
        spmm_sub(pZ, A.mOnes, pY, offset, length, Y.nCols);

        // Note: We DO know that in our case, we can acculumate and do the modulus afterwards
        // This is clearly not generic. Our modulus is 54 bits and max number of +1/-1
        // for one line is 516.
        // So, we delayed this modulo computation to here.                           
        for (Index_t i = 0u; i < realLength; ++i) {
            pZ[i] %= p;
            if (pZ[i] < 0)
                pZ[i] += p;
        }
    }
    }
    
    #pragma omp taskwait
        
    // Note: We DO know that in our case, max "other" in absolute value is 4,
    // So, no complex field tricks required.
    // There's only < 10000 elements in that part, so, yeah, no need for parallel.
    spmm_add(p, Z.data.data(), A.others, pY, Y.nCols);
}

//----------------//
//----- SPMV -----//

//! Z <- Z - A * Y
//! This will NOT compute modulo p
inline
void spmv_add(Element_t* const Z, const SMatrixOnes_t& A, const Element_t* const Y,
              const Index_t offset, const Index_t length)
{
  for (Index_t i = 0u; i < length; ++i)
    for (const auto c : A.at(offset + i))
      Z[i] += Y[c];
}

//! Z <- Z - A * Y
//! This will NOT compute modulo p
inline
void spmv_sub(Element_t* const Z, const SMatrixOnes_t& A, const Element_t* const Y,
              const Index_t offset, const Index_t length)
{

    for (Index_t i = 0u; i < length; ++i)
    for (const auto c : A.at(offset + i))
        Z[i] -= Y[c];
}

//! Z <- A * Y (Y and Z are row-majored)
inline
void spmv(const Modulus_t p, DVector_t& Z, const SMatrix_t& A, const DVector_t& Y)
{
    const Element_t* const pY = Y.data();

    // PG: There is no SIMD used, no reason to use SIMD STEP below
    //
    // Size_t refa = ((A.nRows / DLP_TASKS) / __SIMD_STEP) * __SIMD_STEP;
    // Size_t refb = refa + __SIMD_STEP;
    // Size_t rest = (A.nRows - refa*DLP_TASKS) / __SIMD_STEP;
    // Size_t LaunchTasks = DLP_TASKS;
    // if (refa == 0) {
    //     LaunchTasks = A.nRows/__SIMD_STEP;
    //     refa = __SIMD_STEP;
    //     refb = refa + __SIMD_STEP;
    //     rest = (A.nRows - refa*LaunchTasks);
    //     if (rest) ++LaunchTasks;
    // }

    Size_t refa = (A.nRows / DLP_TASKS);
    Size_t refb = refa + 1;
    Size_t rest = (A.nRows - refa*DLP_TASKS);
    Size_t LaunchTasks = DLP_TASKS;
    if (refa == 0) {
        LaunchTasks = A.nRows;
        refa = 1;
        refb = refa + 1;
        rest = (A.nRows - refa*LaunchTasks);
        if (rest) ++LaunchTasks;
    }
    
    for (Index_t l = 0; l < LaunchTasks; ++l) {
    #pragma omp task shared(Z, A, Y)
    {
//         const Size_t offset = refLength * l;
//         const Size_t length = (l == DLP_TASKS - 1u)? A.nRows - offset : refLength;
        Size_t offset,length;
        if (l<rest) {
            offset = l*refb;
            length = refb;
        } else {
            offset = rest*refb+(l-rest)*refa;
            if (l == LaunchTasks - 1u) {
                length = A.nRows - offset;
            } else {
                length = refa;
            }
        }
        Element_t* const pZ = Z.data() + offset;
    
        // Reset
        memset(pZ, 0, length * sizeof(Element_t));

        // Ones
        spmv_add(pZ, A.ones, pY, offset, length);

        // Mones	
        spmv_sub(pZ, A.mOnes, pY, offset, length);

        // Reduce          
        for (Index_t i = 0u; i < length; ++i) {
            pZ[i] %= p;
            if (pZ[i] < 0)
                pZ[i] += p;
        }
    }
    }
    
    #pragma omp taskwait

    // This does not use SIMD, so no problem here using SPMM version
    spmm_add(p, Z.data(), A.others, pY, 1u);
}


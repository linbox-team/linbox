/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@gmail.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */


#include "common_define.hpp"
#include "common_sdm.hpp"
#include "common_read.hpp"
#include "common_spmm.hpp"
#include "common_dv.hpp"

#include <fflas-ffpack/utils/args-parser.h>
#include <fflas-ffpack/paladin/parallel.h>
#include <recint/recint.h>

#ifdef __USE_128bits
using FieldElement_t = RecInt::rmint<7u, RecInt::MG_INACTIVE>;
#else
using FieldElement_t = RecInt::rmint<6u, RecInt::MG_INACTIVE>;
#endif

// Z <- A * V (mod p)
// Note: global module from rmint have to be set before any call
inline void gemv(DVector_t& Z, const DMatrix_t& A, const DVector_t& V)
{
    memset(Z.data(), 0, Z.size() * sizeof(Element_t));
    
    auto rZ = reinterpret_cast<FieldElement_t*>(Z.data());
    auto rA = reinterpret_cast<const FieldElement_t* const>(A.data.data());
    auto rV = reinterpret_cast<const FieldElement_t* const>(V.data());
    
    const size_t AnRows = A.nRows;
    const size_t AnCols = A.nCols;
    
#if 1
    FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Row,FFLAS::StrategyParameter::Threads> H(DLP_TASKS);


    SYNCH_GROUP ( 
        FOR1D(i, AnRows, H, 
             {
                 for(Index_t j=0u; j< AnCols; ++j)
                     rZ[i] += rA[i*AnCols+j]*rV[j];
             }));

    
#else

    const Size_t refLength = A.nRows / DLP_TASKS;
    for (Index_t l = 0; l < DLP_TASKS; ++l) {
#pragma omp task shared(A)
    {
        const Size_t offset = refLength * l;
        const Size_t length = (l == DLP_TASKS - 1u)? AnRows - offset : refLength;
     
        for (Index_t i = 0u; i < length; ++i) {
            Index_t oi = offset + i;
            for (Index_t j = 0u; j < AnCols; ++j)
                rZ[oi] += rA[oi * AnCols + j] * rV[j];
        }
    }
    }
    
#pragma omp taskwait
#endif
}

// Z <- A * V (mod p)
// Note: global module from rmint have to be set before any call
inline void gemm(DMatrix_t& Z, const DMatrix_t& A, const DMatrix_t& V)
{
    memset(Z.data.data(), 0, Z.nRows*Z.nCols * sizeof(Element_t));
    
    auto rZ = reinterpret_cast<FieldElement_t*>(Z.data.data());
    auto rA = reinterpret_cast<const FieldElement_t* const>(A.data.data());
    auto rV = reinterpret_cast<const FieldElement_t* const>(V.data.data());
    
    const size_t AnRows = A.nRows;
    const size_t AnCols = A.nCols;
    const size_t VnCols = V.nCols;


#if 1
    FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Row,FFLAS::StrategyParameter::Threads> H(DLP_TASKS);


    SYNCH_GROUP ( 
        FOR1D(i, AnRows, H, 
             {
                 for(Index_t j=0u; j< AnCols; ++j)
                     for(Index_t k=0u; k< VnCols; ++k)
                         rZ[i*VnCols+k] += rA[i*AnCols+j]*rV[j*VnCols+k];
             }));
#else
    FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads> H(DLP_TASKS);


    SYNCH_GROUP ( 
        FOR2D(i, j, AnRows, AnCols, H, 
             {
                 for(Index_t k=0u; k< VnCols; ++k)
                     rZ[i*VnCols+k] += rA[i*AnCols+j]*rV[j*VnCols+k];
             }));

    
#endif
}

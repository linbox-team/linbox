/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@imag.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#pragma once

#include "common_common.hpp"

#include <set>
#include <list>
#include <vector>

//------------------//
//----- Timers -----//

#include <fflas-ffpack/utils/timer.h>

using TTimer = FFLAS::OMPTimer;

//--------------------//
//----- Typedefs -----//

// Left projection: a list of line indexes
using LeftProj_t = std::set<Index_t>;

// Dense matrix or vector
using DVector_t = std::vector<Element_t>;

struct DMatrix_t
{
    Size_t nRows = 0u;
    Size_t nCols = 0u;

    DVector_t data; // Size is nRows x nCols, row major order

    DMatrix_t() = default;
    DMatrix_t(const Size_t& rows, const Size_t& cols)
        : nRows(rows)
        , nCols(cols)
        , data(nRows * nCols)
    {}
};

// Sparse vectors
using SVectorOnes_t = std::vector<Index_t>;
using SVector_t = std::vector<std::pair<Index_t, Element_t>>;

// Sparse matrix containing only ones
using SMatrixOnes_t = std::vector<SVectorOnes_t>;

// Sparse matrix containing elements
struct SData_t
{
    Index_t rowIndex;
    Index_t colIndex;
    Element_t data;
};
bool operator==(const SData_t& u, const SData_t& v) {
    std::cerr << u.rowIndex << ' ' << v.rowIndex << std::endl;
    if (u.rowIndex != v.rowIndex) return false;
    if (u.colIndex != v.colIndex) return false;
    if (u.data != v.data) return false;
    return true;
}

using SDataMatrix_t = std::vector<SData_t>;

// Our matrix is mainly -1/1 and sometimes other things.
struct SMatrix_t
{
    Size_t nRows = 0u;
    Size_t nCols = 0u;

    uint64_t nnz = 0u;
    uint64_t nOnes = 0u;
    uint64_t nMOnes = 0u;
    uint64_t nOthers = 0u;

    SMatrixOnes_t ones;
    SMatrixOnes_t mOnes;
    SDataMatrix_t others;
};

bool operator==(const SMatrix_t& A, const SMatrix_t& B) {
    if (A.nRows != B.nRows) return false;
    if (A.nCols != B.nCols) return false;
    if (A.nnz != B.nnz) return false;
    if (A.nOnes != B.nOnes) return false;
    if (A.nMOnes != B.nMOnes) return false;
    if (A.nOthers != B.nOthers) return false;

    if (!(A.ones == B.ones)) return false;
    if (!(A.mOnes == B.mOnes)) return false;
    if (!(A.others == B.others)) return false;
    return true;
}

std::ostream& operator<< (std::ostream& out, const SMatrix_t& A) {
    out << A.nRows << ' ' << A.nCols << " M" << std::endl;
    for(size_t i=0; i<A.ones.size(); ++i)
        for(auto j: A.ones[i])
            out << (i+1) << ' ' << (j+1) << ' ' << "1" << std::endl;
    for(size_t i=0; i<A.mOnes.size(); ++i)
        for(auto j: A.mOnes[i])
            out << (i+1) << ' ' << (j+1) << ' ' << "-1" << std::endl;
    for(auto sda: A.others)
        out << (sda.rowIndex+1) << ' ' << (sda.colIndex+1) << ' ' << (sda.data) << std::endl;
    return out << "0 0 0" << std::endl;
}

void getModulus(const std::string& modulus, Modulus_t& p)
{
    if (modulus.empty()) {
        std::cerr << "/!\\ No modulus specified with option -q." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::stringstream str(modulus);
    str >> p;

    std::cerr << "[GMOD] Using modulus " << p << std::endl;
}


    

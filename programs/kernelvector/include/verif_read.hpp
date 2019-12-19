/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@imag.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#pragma once

#include "verif_define.hpp"

//-----------------//
//----- Tools -----//

//! Will sort the matrix by row, propagating swaps to col and dat
template<typename T, typename E>
void quickSortMatrix(std::vector<T>& row, std::vector<T>& col, std::vector<E>& dat,
                     T left, T right)
{
    T i = left, j = right;
    T pivot = row[ (left + right) / 2];

    while (i <= j)
    {
        while (row[i] < pivot) i++;
        while (row[j] > pivot) j--;

        if (i < j)
        {
            std::swap(row[i], row[j]);
            std::swap(col[i], col[j]);
            std::swap(dat[i], dat[j]);
            i++; j--;
        }
        else if (i == j)
        {
            ++i; if (j>0) --j;
        }
    }

    if (left < j)  quickSortMatrix(row, col, dat, left, j);
    if (i < right) quickSortMatrix(row, col, dat, i, right);
}

//------------------//
//----- Matrix -----//

//! Four elements per line
inline
void readS4OMatrix(const Field_t& field, std::istream& ifs, SMatrix_t& matrix,
                   Index_t& rowDim, Index_t& colDim, Index_t& nnz, bool transposed)
{
    std::cerr << "S4O reader not implemented" << std::endl;
    exit(EXIT_FAILURE);
}

inline
void readSRDMatrix(const Field_t& field, std::istream& ifs, SMatrix_t& matrix,
                   Index_t& rowDim, Index_t& colDim, Index_t& nnz, bool transposed)
{
    std::vector<Index_t> iRow, iCol;
    std::vector<FieldElement_t> iDat;
    
    // Header
    ifs >> rowDim >> colDim;
    nnz = 0u;
    
    // Content
    Index_t nElements;
    Index_t row = 0u, col;
    FieldElement_t data;
    
    while (ifs >> nElements) {
        for (Index_t i = 0u; i < nElements; ++i) {
            ifs >> col;
            ifs.ignore(1);   // Separator character ':'
            field.read(ifs, data);
            
            iRow.emplace_back(row);
            iCol.emplace_back(col);
            iDat.emplace_back(data);
        }

        ++row;
        nnz += nElements;
    }
    
    iRow.shrink_to_fit();
    iCol.shrink_to_fit();
    iDat.shrink_to_fit();

    // Apply transpose if needed
    if (transposed)
    {
        std::swap(iRow, iCol);
        std::swap(rowDim, colDim);

        // Re-sort by row (implicitly asked by sparse_init())
        quickSortMatrix(iRow, iCol, iDat, Index_t(0u), Index_t(nnz - 1u));
    }
    
    // Really init the matrix
    FFLAS::sparse_init(field, matrix, iRow.data(), iCol.data(), iDat.data(), rowDim, colDim, nnz);
}

template<class FE_Container>
std::ostream& mprint(std::ostream& os, const FE_Container& M, const Field_t& field) {
    os << '['; 
    for(const auto& val: M) 
        field.write(os, val) << ' '; 
    return os << ']';
}  


// Note: The SMS have to be row-sorted.
inline
void readSMSMatrix(const Field_t& field, std::istream& ifs, SMatrix_t& matrix,
                   Index_t& rowDim, Index_t& colDim, Index_t& nnz, bool transposed)
{
    std::vector<Index_t> iRow, iCol;
    std::vector<FieldElement_t> iDat;
    
    // Header
    std::string tmp;
    ifs >> rowDim >> colDim >> tmp;
    nnz = 0u;
    
    // Content
    Index_t row, col;
    FieldElement_t data;
    
    while (ifs >> row) {
        // Stop on 0 0 0
        if (row == 0u)
            break;
            
        // Read data
        ifs >> col;
        field.read(ifs, data);
    
        // Save
        iRow.emplace_back(row - 1u);
        iCol.emplace_back(col - 1u);
        iDat.emplace_back(data);
        
        ++nnz;
    }

    // Apply transpose if needed
    if (transposed)
    {
        std::swap(iRow, iCol);
        std::swap(rowDim, colDim);

        // Re-sort by row (implicitly asked by sparse_init())
        quickSortMatrix(iRow, iCol, iDat, Index_t(0u), Index_t(nnz - 1u));
    }
    
    // Really init the matrix
    FFLAS::sparse_init(field, matrix, iRow.data(), iCol.data(), iDat.data(), rowDim, colDim, nnz);
// mprint(std::cerr << "A: ", iDat, field) << std::endl;
}


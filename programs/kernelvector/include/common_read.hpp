
/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@imag.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#pragma once

#include "common_define.hpp"

#include <algorithm> // std::sort

//------------------//
//----- Matrix -----//

inline
void reinitMatrix(SMatrix_t& matrix)
{
    matrix.nnz = 0u;
    matrix.nOnes = 0u;
    matrix.nMOnes = 0u;
    matrix.nOthers = 0u;
    
    matrix.ones.clear();
    matrix.mOnes.clear();
    matrix.others.clear();
}

//! Four elements per line
// FIXME -0 is not correctly interpreted (because it is read as 0)
inline
void readS4OMatrix(Modulus_t p, std::istream& ifs, SMatrix_t& matrix)
{
    // Reinit
    reinitMatrix(matrix);
    
    // Header
    ifs >> matrix.nRows >> matrix.nCols;
    
    matrix.ones.resize(matrix.nRows);
    matrix.mOnes.resize(matrix.nRows);
    
    // S40 indexes are stored in hexidecimal
    ifs >> std::hex;

    // Content
    const uint8_t nElements = 4u;
    int64_t colID;

    for (Index_t row = 0u; row < matrix.nRows; ++row) {
        for (uint8_t element = 0u; element < nElements; ++element) {
            ifs >> colID;
            
            if (colID < 0) {
                ++matrix.nMOnes;
                matrix.mOnes[row].emplace_back(static_cast<Index_t>(-colID));
            }
            else {
                ++matrix.nOnes;
                matrix.ones[row].emplace_back(static_cast<Index_t>(colID));
            }
            
            // Skip comma
            if (element != nElements - 1u)
                ifs.ignore(1);
        }
        
        // Gain some space
        matrix.ones[row].shrink_to_fit();
        matrix.mOnes[row].shrink_to_fit();

        matrix.nnz += nElements;
    }
}

//! Four elements per line
// FIXME -0 is not correctly interpreted (because it is read as 0)
inline
void readCADOMatrix(Modulus_t p, std::istream& ifs, SMatrix_t& matrix)
{
    // Reinit
    reinitMatrix(matrix);
    
    // indexes are stored in hexidecimal
    ifs >> std::hex;

    // Content with 4 per row
    const size_t nElements(4);
    uint64_t colID[nElements];
    Element_t a,b;
    char c;
    
//     ifs >> a >> c >> b >> c >> colID[0] >> c >> colID[1] >> c >> colID[2] >> c >> colID[3];
    ifs >> a >> c >> b;
    while( !ifs.eof() ) {
//         std::cerr << std::dec << a << ',' << b << ' ' 
//                   << colID[0] << ',' 
//                   << colID[1] << ','
//                   << colID[2] << ',' 
//                   << colID[3] << std::endl;
       matrix.ones.resize(matrix.ones.size()+1);
        matrix.mOnes.resize(matrix.mOnes.size()+1);
        for (uint8_t element = 0u; element < nElements; ++element) {
            ifs >> c >> c;
            bool negative(false);
            if (c == '-') negative=true;
            else ifs.putback(c);
            ifs >> colID[element];
            
            if (negative){
                ++matrix.nMOnes;
                matrix.mOnes[matrix.nRows].emplace_back(static_cast<Index_t>(colID[element]));
                if (colID[element]>matrix.nCols) matrix.nCols=colID[element];
            }
            else {
                ++matrix.nOnes;
                matrix.ones[matrix.nRows].emplace_back(static_cast<Index_t>(colID[element]));
                if (colID[element]>matrix.nCols) matrix.nCols=colID[element];
            }
        }       
        // Gain some space
        matrix.ones[matrix.nRows].shrink_to_fit();
        matrix.mOnes[matrix.nRows].shrink_to_fit();

        matrix.nnz += nElements;

        ++matrix.nRows;
        ifs >> a >> c >> b;        
    }

    ++matrix.nCols;
    
}

// Important: We need to be in a parallel region before calling this
inline
std::istream& readSRDMatrix(Modulus_t p, std::istream& ifs, SMatrix_t& matrix)
{
    // Reinit
    reinitMatrix(matrix);
    
    // Header
    ifs >> matrix.nRows >> matrix.nCols;
    
    matrix.ones.resize(matrix.nRows);
    matrix.mOnes.resize(matrix.nRows);

    // Content
    Index_t col;
    Index_t row = 0u;
    Element_t data;
    Element_t one(1u), mOne(-1);
    Size_t nElements;

    while (ifs >> nElements) {
        // For each element on line, add it to the correct list
        for (Index_t i = 0u; i < nElements; ++i) {
            ifs >> col;
            ifs.ignore(1);   // Separator character ':'
            ifs >> data;
            
            // Data is 1
            if (data == one) {
                ++matrix.nOnes;
                matrix.ones[row].emplace_back(col);
            }
            
            // Data is -1
            else if (data == mOne) {
                ++matrix.nMOnes;
                matrix.mOnes[row].emplace_back(col);
            }
            
            // Data is something else
            else {
                ++matrix.nOthers;
                
                // Add the data
                SData_t sData;
                sData.rowIndex = row;
                sData.colIndex = col;
                sData.data = data;
                matrix.others.emplace_back(std::move(sData));
            }
        }
        
        // Gain some space
        matrix.ones[row].shrink_to_fit();
        matrix.mOnes[row].shrink_to_fit();
        
        // Sort indexes
        // Well, in fact we're slower during the computation if we do that...
        /*std::sort(std::begin(matrix.ones[row]),  std::end(matrix.ones[row]));
        std::sort(std::begin(matrix.mOnes[row]), std::end(matrix.mOnes[row]));*/
        
        matrix.nnz += nElements;
        ++row;
    }
    
    matrix.others.shrink_to_fit();
    
    if (row != matrix.nRows) {
        std::cerr << "/!\\ Effective number of matrix rows does not match size." << std::endl;
        std::cerr << "    The matrix has really " << row << " rows but indicated " << matrix.nRows << "." << std::endl;
    }
    
    
    /* This has not a great impact in the code, just using numactl --interleave=all is better.
    // Mapping memory to different CPU
    // We force data locality to each core
    // Program have to be launched with GOMP_CPU_AFFINITY envvar set.
    const Size_t refLength = matrix.nRows / DLP_TASKS;
    
    for (Index_t l = 0u; l < DLP_TASKS; ++l)
    {
    #pragma omp task shared(matrix)
    {
        const Size_t offset = l * refLength;
        const Size_t length = (l == DLP_TASKS - 1u)? matrix.nRows - offset : refLength;
        
        for (Index_t k = 0; k < length; ++k) {
            std::vector<Index_t> tmp1, tmp2;
            tmp1 = matrix.ones[offset + k];
            tmp2 = matrix.mOnes[offset + k];
            
            matrix.ones[offset + k].swap(tmp1);
            matrix.mOnes[offset + k].swap(tmp2);
        }
    }
    }
    
    #pragma omp taskwait*/
    return ifs;
}

// Note: The SMS have to be row-sorted.
inline
std::istream& readSMSMatrix(Modulus_t p, std::istream& ifs, SMatrix_t& matrix)
{
    // Reinit
    reinitMatrix(matrix);
    
    // Header
    std::string tmp;
    ifs >> matrix.nRows >> matrix.nCols >> tmp;
    
    matrix.ones.resize(matrix.nRows);
    matrix.mOnes.resize(matrix.nRows);

    // Content
    Index_t row, cRow;
    Index_t col;
    Element_t data;
    Element_t one(1u), mOne(-1);

    ifs >> row;

    while (true) {
        // Stop on 0 0 0
        if (row == 0u)
            break;
            
        // Current row
        cRow = row;
    
        // Read until new row index
        while (row == cRow) {
            ifs >> col;
            ifs >> data;
            
            // Data is 1
            if (data == one) {
                ++matrix.nOnes;
                matrix.ones[row - 1u].emplace_back(col - 1u);
            }
            
            // Data is -1
            else if (data == mOne) {
                ++matrix.nMOnes;
                matrix.mOnes[row - 1u].emplace_back(col - 1u);
            }
            
            // Data is something else
            else {
                ++matrix.nOthers;
                
                // Add the data
                SData_t sData;
                sData.rowIndex = row - 1u;
                sData.colIndex = col - 1u;
                sData.data = data;
                matrix.others.emplace_back(std::move(sData));
            }
            
            ++matrix.nnz;
            ifs >> row;
        }
        
        // Gain some space
        matrix.ones[cRow - 1u].shrink_to_fit();
        matrix.mOnes[cRow - 1u].shrink_to_fit();
    }
    
    matrix.others.shrink_to_fit();
    return ifs;
}

//-------------------//
//----- Getters -----//

int getMatrix(Modulus_t p, const std::string& matrixFile, SMatrix_t& matrix, bool checkSquared=true)
{
    // Select stream
    std::ifstream matrixStream(matrixFile);
    if (!matrixStream.is_open()) {
        std::cerr << "/!\\ No valid matrix file specified with option -f." << std::endl;
        exit(EXIT_FAILURE);
    }
    
    std::cerr << "[GMAT] Reading matrix from file " << matrixFile << "..." << std::endl;

    TTimer timer;
    timer.clear(); timer.start();
    
    auto extension = matrixFile.substr(matrixFile.find_last_of(".") + 1);
    if (extension == "srd")
        readSRDMatrix(p, matrixStream, matrix);
    else if (extension == "s4o")
        readS4OMatrix(p, matrixStream, matrix);
    else if (extension == "cado")
        readCADOMatrix(p, matrixStream, matrix);
    else if (extension == "sms")
        readSMSMatrix(p, matrixStream, matrix);
    else {
        std::cerr << "/!\\ Unknown extension: " << extension << std::endl;
        std::cerr << "    Please use lowercase and sms, srd or s4o matrix file format." << std::endl;
        exit(EXIT_FAILURE);
    }
    
    timer.stop();
    
    std::cerr << "[GMAT] Read matrix " << matrix.nRows << "x" << matrix.nCols << " in ";
    std::cerr << timer.usertime() << "s." << std::endl;
    std::cerr << "[GMAT] Matrix has " << matrix.nnz << " nnz: ";
    std::cerr << matrix.nOnes << " ones " << matrix.nMOnes << " mOnes and " << matrix.nOthers << " others" << std::endl;

    if (checkSquared && (matrix.nRows != matrix.nCols)) {
        std::cerr << "/!\\ Matrix is not squared." << std::endl;
        return (EXIT_FAILURE);
    }
    
    return 0;
}


/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@gmail.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#pragma once
#include <recint/rint.h>
#include <recint/ruint.h>

#include <cstdint>
#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip> // std::setw()

//--------------------//
//----- Typedefs -----//

#ifdef __USE_128bits
using Element_t = RecInt::rint128;
using Modulus_t = RecInt::rint128;
#define __DLP_ABS(a) ((a)<0?-(a):(a))
#else
using Element_t = int64_t;      // Sufficiently large to store a (field) element.
using Modulus_t = int64_t;     // Sufficiently large to store the modulus.
#define __DLP_ABS(a) std::abs(a)
#endif
using Index_t = uint64_t;       // Sufficiently large to store an array index.
using Size_t = uint64_t;        // Sufficiently large to store an array size.

//-----------------//
//----- Utils -----//

inline std::string formatNumber(const Index_t n)
{
    std::stringstream str;
    // Note: 20 is base ten max number of digits for a 64-bits number.
    str << "_" << std::setw(20) << std::setfill('0') << n;
    return str.str();
}

//----------------//
//----- Info -----//

inline
void saveInfo(const std::string& filename, const std::string& matrixFile, Modulus_t q,
              Size_t m, Size_t n, Size_t s1, Size_t s2, Size_t d, Size_t K)
{
    std::ofstream infoFile(filename);
    if (!infoFile.is_open())
        throw std::runtime_error("Cannot open info file.");

    infoFile << "matrixFile: " << matrixFile << std::endl;
    infoFile << "modulus: " << q << std::endl;
    infoFile << "nRows: " << m << std::endl;
    infoFile << "nCols: " << n << std::endl;
    infoFile << "s1: " << s1 << std::endl;
    infoFile << "s2: " << s2 << std::endl;
    infoFile << "sequenceLength: " << d << std::endl;
    infoFile << "checkpointsStep: " << K << std::endl;
}

inline
void loadInfo(const std::string& filename, std::string& matrixFile, Modulus_t& q,
              Size_t& m, Size_t& n, Size_t& s1, Size_t& s2, Size_t& d, Size_t& K)
{
    std::ifstream infoFile(filename);
    if (!infoFile.is_open())
        throw std::runtime_error("Cannot open info file.");

    std::string tmp;
    infoFile >> tmp >> matrixFile;
    infoFile >> tmp >> q;
    infoFile >> tmp >> m;
    infoFile >> tmp >> n;
    infoFile >> tmp >> s1;
    infoFile >> tmp >> s2;
    infoFile >> tmp >> d;
    infoFile >> tmp >> K;
}


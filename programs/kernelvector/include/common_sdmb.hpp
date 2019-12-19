/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@imag.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#pragma once

#include <vector>
#include <iostream>

//----------------//
//----- SDMB -----//

void writeSDMBHeader(std::ostream& os, uint64_t d, uint64_t m, uint64_t n)
{
    os.write(reinterpret_cast<const char*>(&m), sizeof(m));
    os.write(reinterpret_cast<const char*>(&n), sizeof(n));
    os.write(reinterpret_cast<const char*>(&d), sizeof(d));
}

void readSDMBHeader(std::istream& is, uint64_t& d, uint64_t& m, uint64_t& n)
{
    is.read(reinterpret_cast<char*>(&m), sizeof(m));
    is.read(reinterpret_cast<char*>(&n), sizeof(n));
    is.read(reinterpret_cast<char*>(&d), sizeof(d));
}

template<class T>
inline void writeSDMBMatrix(std::ostream& os, const std::vector<T>& M)
{
    os.write(reinterpret_cast<const char*>(M.data()), M.size() * sizeof(T));
}

template<class T>
inline void readSDMBMatrix(std::istream& is, std::vector<T>& M)
{
    is.read(reinterpret_cast<char*>(M.data()), M.size() * sizeof(T));
}


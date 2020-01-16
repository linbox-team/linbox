/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@imag.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#pragma once

#include <vector>
#include <iostream>
#include "common_common.hpp"
// Left projection: a list of line indexes
using Size_t = uint64_t;
using Index_t = uint64_t;
using LeftProj_t = std::set<Index_t>;

//----------------//
//----- SDMP -----//

std::ostream& writeSDMPHeader(std::ostream& os, uint64_t d, uint64_t m, uint64_t n)
{
    return os << m << ' ' << n << ' ' << d << std::endl;
}

std::istream& readSDMPHeader(std::istream& is, uint64_t& d, uint64_t& m, uint64_t& n)
{
    return is >> m >> n >> d;
}

std::istream& readSDMPHeader(std::istream& is, uint64_t& d)
{
    return is >> d;
}

void getProj(std::istream& projStream, LeftProj_t& U, const Size_t s1)
{
    Index_t r;
    // Header
    projStream >> r;
    if (r != s1)
        throw std::logic_error("Size conflict between info and left projection files.");

    // Content
    for (Index_t i = 0u; i < s1; ++i) {
        projStream >> r;
        U.insert(r);
    }
}


template<class T>
inline std::ostream& writeSDMPMatrix(std::ostream& os, const std::vector<T>& M)
{
    for (const auto& a : M)
        os << a << ' ';
    return os << std::endl;
}

template<class T>
inline std::istream& readSDMPMatrix(std::istream& is, std::vector<T>& M)
{
    for (auto& a : M)
        is >> a;
    return is;
}

template<class T>
inline std::istream& readSDMPVector(std::istream& is, std::vector<T>& V, const size_t stride)
{
    T ignore; const size_t skip=stride-1;
    for (auto& a : V) {
        is >> a;
        for(size_t ii=0;ii<skip;++ii) is >> ignore;
    }
    return is;
}

template<class T, class Field>
inline std::ostream& writeSDMPMatrix(const Field& field, std::ostream& os, const std::vector<T>& M)
{
    for (const auto& a : M)
        field.write(os, a) << ' ';
    return os << std::endl;
}

template<class T, class Field>
inline std::ostream& writeSDMPMatrixMaple(const Field& field, std::ostream& os, const std::vector<T>& M, const size_t m)
{
    os << '[';
    for(size_t i=0; i<m; ++i) {
        os << '[';
        for(size_t j=0; j<m; ++j) {
            field.write(os, M[i*m+j]);
            if (j != (m-1)) os << ',';
        }
        os << ']';
        if (i != (m-1)) os << ',';
        os << std::endl;
    }
    return os;
}

template<class T, class Field>
inline std::istream& readSDMPMatrix(const Field& field, std::istream& is, std::vector<T>& M)
{
    for (auto& a : M)
        field.read(is, a);
    return is;
}


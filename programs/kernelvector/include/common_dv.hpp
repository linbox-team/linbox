/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@gmail.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#pragma once

#include "common_define.hpp"

inline
void readDV(std::istream& stream, DVector_t& V)
{ 
    Size_t size;
    stream >> size;
    V.resize(size);
    
    for (auto& v : V)
        stream >> v;
}

inline
void writeDV(std::ostream& stream, const DVector_t& V)
{
    stream << V.size() << std::endl;
    for (const auto& v : V)
        stream << v << std::endl;
}

template <class Field, class Element> inline
void writeDV(const Field& field, std::ostream& stream, const std::vector<Element>& V)
{
    stream << V.size() << std::endl;
    for (const auto& v : V)
        field.write(stream, v) << std::endl;
}

template <class Vector_t> inline
bool isNullVector(const Modulus_t& p, const Vector_t& v) {
//     for(auto vals:v) {
//         if (! (vals%p)) {
//             std::cerr << "value: " << vals << " mod " << p << " is " << (vals%p) << " test: " << (! (vals%p)) << ", t2: " << (! (vals==Modulus_t(0))) << std::endl;
//             return false;
//         }
//     }
//     for(auto vals:v) {
//         if (! ( (vals%p)==Modulus_t(0))) {
//             std::cerr << "value: " << vals << " mod " << p << " is " << (vals%p) << " test: " << (! (vals%p)) << ", t2: " << (! (vals==Modulus_t(0))) << std::endl;
//             return false;
//         }
//     }
//     for(auto vals:v) if (! (vals%p)) return false;
    for(auto vals:v) if (! ( (vals%p)==Modulus_t(0))) return false;
    return true;
}
    
#include <givaro/udl.h>
template<typename Element>
bool isFirstCanonical(const std::vector<Element>& v) {
    const Element elementOne(1u), elementZero(0u);
    auto iter(v.begin());
    if ( static_cast<Element>(*iter) != elementOne) return false;
    for(++iter ; iter != v.end(); ++iter) {
        if (static_cast<Element>(*iter) != elementZero) return false;
    }
    return true;
}
    

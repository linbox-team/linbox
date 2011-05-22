/* Copyright (C) 2010 LinBox
 * Written by JG Dumas <Jean-Guillaume.Dumas@imag.fr>
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

// Pair of I and T : struct { column index, value }
// ======================================================================= 
#ifndef __LINBOX_pair_H
#define __LINBOX_pair_H
#include <iostream>

// ---------------------------------------------------
//
/// Pair of I and T : struct { column index, value }
template<class I, class T> class Pair {
public:
    typedef Pair<I, T> Self_t;
    typedef T                      Type_t;
    typedef I                      first_type;
    typedef T                      second_type;
    
//    ~Pair() {};
    Pair() {};

    Pair(const I jj, const T& val) :first(jj),second(val){};
    Pair(const Self_t& p) :first(p.first),second(p.second){};


    T getvalue() const { return second; };

    I getindex() const { return first; };
    I j() const { return first; };
    
    T affect(const T& val) { return second = val; };
    T change_value(const T& val) { return second = val; };
   
    I change_j(const I jj) { return first = jj; };      
    I change_index(const I jj) { return first = jj; };      
            
            
    Self_t assign(const T& val) {
        second = val;
        return *this;
    };      
            
    Self_t assign(const I jj, const T& val) {
        second = val;
        first = jj;
        return *this;
    };      
            
    I decr() { return --first; };      
    I operator--() { return --first; };      
    I operator--(int) { return first--; };      
    I incr() { return ++first; };      
    I operator++() { return ++first; };      
    I operator++(int) { return first++; };      
            
    friend inline std::istream& operator>> (std::istream& is, Pair<I, T>& a) {
        is >> a.first >> a.second;
        return is;
    }
    
    friend inline std::ostream& operator<< (std::ostream& o, const Pair<I, T> a){
        return o << a.first << " " << a.second ;
    }
    

public:
    I first;
    T second;
};



#endif // __LINBOX_pair_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

/* linbox/solutions/minpoly.h
 * Copyright (C) 2000, 2010 Jean-Guillaume Dumas
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
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

// ======================================================================= 
// Sparse Vector      : vector< Pair<T> > and an additional actual size
// ======================================================================= 
#ifndef __LINBOX_sparse_vector_H
#define __LINBOX_sparse_vector_H
#include <iostream>

#include <linbox/vector/vector-traits.h>

// ---------------------------------------------------
//
/// Default container for sparse vectors is LightContainer
#ifndef _IBB_VECTOR_
// #include <vector>
// #define _IBB_VECTOR_ std::vector
#include "linbox/vector/light_container.h"
#define _IBB_VECTOR_ LightContainer
#endif // _IBB_VECTOR_
#ifndef _IBB_PAIR_
#include <utility>
#define _IBB_PAIR_ std::pair
// #include "linbox/vector/pair.h"
// #define _IBB_PAIR_ Pair
#endif // _IBB_PAIR_



namespace LinBox
{
// ---------------------------------------------------
//
/** \brief vector< Pair<T,I> > and actualsize
\ingroup vector
*/
template<class T, class I = unsigned int>
class Sparse_Vector : public _IBB_VECTOR_< _IBB_PAIR_<I, T> > {
public:
    typedef _IBB_PAIR_<I, T>       Element;
    typedef T                      Type_t;
    typedef Sparse_Vector<T, I>    Self_t;

    // Dan Roche 6-30-04
    typedef VectorCategories::SparseSequenceVectorTag VectorCategory;


    Sparse_Vector() {};
    Sparse_Vector(size_t n) : _IBB_VECTOR_< _IBB_PAIR_<I, T> >(n), _rsize(0) {};
    Sparse_Vector(size_t n, size_t rn) : _IBB_VECTOR_< _IBB_PAIR_<I, T> >(n), _rsize(rn) {};
    ~Sparse_Vector() {};
    
            

        /// Actual dimension of the vector, 0 for infinite or unknown
    inline size_t actualsize() const { return _rsize; };            
    inline size_t reactualsize( const size_t s ) { return _rsize = s; };
    template<class XX> inline size_t reactualsize( const XX s ) { return _rsize = (size_t)s; };

    friend inline std::ostream& operator<< (std::ostream& o, const Sparse_Vector<T, I> v) {
        if (v.size())
            for(long i=0;i<v.size();i++)
                o << v[i] << std::endl;
        return o;
    }


private:
    size_t _rsize;
};    

} //end of namespace LinBox

#endif // __LINBOX_sparse_vector_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

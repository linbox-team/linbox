/* Copyright (C) 2010 LinBox
 * Written by JG Dumas <Jean-Guillaume.Dumas@imag.fr>
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


// =================================================================== //
// LightContainer : std::vector like container
// =================================================================== //

#ifndef __LINBOX_light_vector_container_H 
#define __LINBOX_light_vector_container_H 

#include <iostream>
#include <cstdlib>
#include "linbox/util/contracts.h"
#include "linbox/vector/vector-traits.h"


namespace LinBox 
{

template<typename Elem> struct LightContainer {
private:
    typedef LightContainer<Elem> Self_t;
    size_t allocated;
    Elem * _container;
    Elem * _finish;
public:
    typedef Elem value_type;
    typedef Elem* iterator;
    typedef const Elem* const_iterator;

    LightContainer() : allocated(2), _container(new Elem[allocated]), _finish(_container) { 
        ENSURE( (allocated == 2) && (size() == 0) );

    }
    LightContainer(size_t s) : allocated(s), _container(new Elem[s]), _finish(_container+s) {
        
        ENSURE( (allocated == s) && (size() == s) );
        ENSURE( allocated >= size() );
     
    }

    Self_t& operator=(const Self_t& v) {
        this->resize(v.size());
        const_iterator vit = v.begin();
        for(iterator tit = begin(); tit != end(); ++tit, ++vit)
            *tit = *vit;
        return *this;
    }

    LightContainer(const Self_t& v) : allocated(v.allocated) {
        _container = new Elem[allocated];
        _finish = _container +v.size();
        const_iterator vit = v.begin();
        for(iterator tit = begin(); tit != end(); ++tit, ++vit)
            *tit = *vit;
    }

    void reserve(size_t s) { 
        STATE( size_t oldsize = size() );
        reallocate(s, size() );
        ENSURE( (allocated >= s) && (size() == oldsize) && ( allocated >= size() ) );
    }

    size_t size() const { return size_t(_finish - _container); }

    Elem& operator[](size_t i) { 
        REQUIRE( (i >= 0) && (i < allocated) && (i < size()) );
        return _container[i]; }

    const Elem& operator[](size_t i) const { 
        REQUIRE( (i >= 0) && (i < allocated) && (i < size()) );
        return _container[i]; }

    void clear() { _finish = _container; 
            ENSURE( (size() == 0) );

    }
    void resize(size_t s) {
        if (s>allocated) reallocate( s+(s>>1), s );
        else _finish = _container + s;
        ENSURE( allocated >= size() );
    }
    iterator begin() { return _container; }
    iterator end() { return _finish; }
    const_iterator begin() const { return const_iterator(_container); }
    const_iterator end() const { return const_iterator(_finish); }
    Elem& front() { return *_container; }
    const Elem& front() const { return *_container; }
    Elem& back() { return *(_finish-1); }
    const Elem& back() const { return *(_finish-1); }

    void push_back(const Elem& c) {
        STATE( size_t  oldsize = size() );
        if (size() == allocated) reserve(allocated+(allocated>>1));
        *(_finish) = c; ++_finish; 

        ENSURE( size() == (oldsize+1) );
        ENSURE( allocated >= size() );
    }
    void pop_back() {
        STATE( size_t  oldsize = size() );
        REQUIRE( oldsize >= 1 );
        --_finish;
        ENSURE( size() == (oldsize-1) );
        ENSURE( allocated >= size() );
    }

    ~LightContainer() { delete[] _container; }


    iterator insert(iterator pos, const Elem& c) {
        REQUIRE( (pos-begin()) <= (end()-begin()) );
        REQUIRE( (pos-begin()) >= 0 );
        STATE( size_t oldsize = size() );
        iterator newpos;
        if (pos == _finish) {
            push_back(c);
            newpos = _finish-1;
        } else {
            if (allocated > size())
                newpos = insertwithspace(pos,c);
            else
                newpos = insertwithrealloc(pos,c);
        }
        ENSURE( size() == oldsize+1 );
        ENSURE( allocated >= size() );
        return newpos;
    }
    
    iterator insert(iterator pos, const_iterator beg, const_iterator end) {
        REQUIRE( (pos-begin()) <= (end()-begin()) );
        REQUIRE( (pos-begin()) >= 0 );
        if (pos == _finish) {
            for(const_iterator iter=beg; iter != end; ++iter)
                push_back(*iter);
        } else {
            for(const_iterator iter=beg; iter != end; ++iter)
                insert(pos, *iter);
        }
        ENSURE( allocated >= size() );
        return pos;
    }


    iterator erase(iterator pos) {
        REQUIRE( (pos-begin()) < (end()-begin()) );
        REQUIRE( (pos-begin()) >= 0 );
        STATE( size_t oldsize = size() );
        iterator ppos=pos+1;
        if ( ppos == _finish) {
            pop_back();
        } else {
            *(pos)=*(ppos);
            erase(ppos);
        }
        ENSURE( size() == oldsize-1 );
        ENSURE( allocated >= size() );
        ENSURE( _finish >= _container );
        return pos;
    }
    
    iterator erase(iterator first, iterator last) {
        REQUIRE( (first-begin()) < (end()-begin()) );
        REQUIRE( (last-begin()) <= (end()-begin()) );
        REQUIRE( (last-first) > 0 );
        REQUIRE( (first-begin()) >= 0 );
        STATE( size_t oldsize = size() );
        const size_t lmf = last-first;
        if ( last == _finish) {
            for(size_t i=0; i<lmf; ++i) pop_back();
        } else {
            iterator nf(first), nl(last);
            for( ; (nf != last) && (nl != _finish); ++nf, ++nl)
                *(nf) = *(nl);
            erase(nf,nl);
        }
        ENSURE( size() == oldsize-lmf );
        ENSURE( allocated >= size() );
        ENSURE( _finish >= _container );
        return first;
    }
    
/*
    friend std::ostream& operator<< (std::ostream& o, const Self_t& C) {
        o << '[';
        const_iterator refs =  C.begin();
        for( ; refs != (C.end()-1) ; ++refs )
            o << (*refs) << ',';
        return o << (*refs) << ']';
    }
    friend std::ostream& operator<< (std::ostream& o, const Self_t& C) {
        o << '[';
        for(size_t i=0; i<(C.size()-1); ++i) 
            o << C[i] << ',';
        return o << C[C.size()-1] << ']';
    }
*/

    friend std::ostream& operator<< (std::ostream& o, const Self_t& C) {
        o << '[';
        for(const_iterator refs =  C.begin(); refs != C.end() ; ++refs )
            o << (*refs) << ',';
        return o << ']';
    }



protected:
    void reallocate(size_t s, size_t endc) {
        REQUIRE( (s >= endc) );
        if (allocated < s) {
            Elem * futur = new Elem[s];
            for(size_t i=0; (i<s) && (i < allocated); ++i)
                futur[i] = _container[i];
            size_t olds = size();
            delete [] _container;
            _container = futur;
            _finish = _container + olds;
            allocated = s;
        }
        _finish = _container + endc;
        ENSURE( allocated >= size() );
    }

    iterator insertwithspace(iterator pos, const Elem& c) {
        STATE( size_t oldsize = size() );
        STATE( size_t oldalloc = allocated );
        REQUIRE( (size()+1) <= allocated );
        REQUIRE( (pos-begin()) <= (end()-begin()) );
        REQUIRE( (pos-begin()) >= 0 );
        if (pos == _finish) {
            push_back(c);
        } else {
            insertwithspace(pos+1, *(pos));
            *(pos) = c;
        }
        ENSURE( size() == oldsize+1 );
        ENSURE( allocated >= size() );
        ENSURE( allocated == oldalloc );
        return pos;
    }
    
    iterator insertwithrealloc(iterator pos, const Elem& c) {
        REQUIRE( (pos-begin()) <= (end()-begin()) );
        REQUIRE( (pos-begin()) >= 0 );
        allocated += (allocated>>1);
        Elem * futur = new Elem[allocated];
        iterator newcont = futur;
        iterator oldcont=_container;
        for( ; oldcont != pos; ++oldcont,++newcont)
            *newcont = *oldcont;
        *newcont = c; 
        iterator newpos = newcont;
        for(++newcont ; oldcont != _finish; ++oldcont,++newcont)
            *newcont = *oldcont;
        size_t olds = size();
        delete [] _container;
        _container = futur;
        _finish = _container + (++olds);
        ENSURE( allocated >= size() );
        return newpos;
    }
   


};  


    // Specialization for LightContainer
template <class Element>
struct VectorTraits< LightContainer<Element> >
{ 
    typedef LightContainer<Element> VectorType;
    typedef typename VectorCategories::DenseVectorTag VectorCategory; 
};

    // Specialization for LightContainer of pairs of size_t and elements
template <class Element> 
struct VectorTraits< LightContainer< std::pair<size_t, Element> > >
{ 
    typedef LightContainer< std::pair<size_t, Element> > VectorType;
    typedef typename VectorCategories::SparseSequenceVectorTag VectorCategory; 
    
    static void sort (VectorType& v) { std::stable_sort(v.begin(), v.end(), SparseSequenceVectorPairLessThan<Element>()); }
};

} // namespace LinBox

#endif //__LINBOX_light_vector_container_H 

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

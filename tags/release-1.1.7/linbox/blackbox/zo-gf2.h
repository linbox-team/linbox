/* linbox/blackbox/zo-gf2.h
 * Copyright (C) 2009,2010 The LinBox group
 *
 * Time-stamp: <23 Jun 10 16:02:02 Jean-Guillaume.Dumas@imag.fr> 
 *
 * See COPYING for license information.
 *
 */
#ifndef __LINBOX_zo_gf2_H
#define __LINBOX_zo_gf2_H

#include <algorithm>
#include "linbox/blackbox/zero-one.h"
#include "linbox/field/gf2.h"
#include "linbox/field/unparametric.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/vector/stream.h"
#include "linbox/matrix/sparse.h"
#include "linbox/vector/light_container.h"

namespace LinBox 
{
    
 /** \brief Time and space efficient representation of sparse matrices over GF2.
   * Representation if a full row array containing vector of non-zero locations
   * A 0-1 matrix is a matrix with all 0's and 1's as entries.  
   \ingroup blackbox
  */ 
    template<>
    struct ZeroOne<GF2> : public LightContainer< LightContainer< size_t > > 
    {
        typedef LightContainer< LightContainer< size_t > > Father_t;
        typedef LightContainer< size_t > Row_t;
        typedef GF2::Element Element;
        typedef size_t Index;
        typedef ZeroOne<GF2> Self_t;
        typedef GF2 Field;

        const GF2 _F;

        ZeroOne(const GF2& ) : _nnz(0) {}
        ZeroOne(const GF2& , const size_t m) : Father_t(m), _rowdim(m), _coldim(m),_nnz(0) {}
        ZeroOne(const GF2& , const size_t m, const size_t n) : Father_t(m), _rowdim(m), _coldim(n),_nnz(0) {}

        ZeroOne(): _nnz(0)  {}
        ZeroOne(const size_t m) : Father_t(m), _rowdim(m), _coldim(m),_nnz(0) {}
        ZeroOne(const size_t m, const size_t n) : Father_t(m), _rowdim(m), _coldim(n),_nnz(0) {}

        ZeroOne(const GF2& , VectorStream<Row_t>& stream) : Father_t(stream.m()), _rowdim(stream.m()), _coldim(stream.n()), _nnz(0) {
            for (Father_t::iterator row=begin(); row != end(); ++row) {
                stream >> *row;
		_nnz += row->size();
            }
        }        

        ZeroOne(const Self_t& A) : Father_t(static_cast<const Father_t&>(A)), _rowdim(A._rowdim), _coldim(A._coldim), _nnz(A._nnz) {
        }

        ZeroOne(const GF2& , size_t* rowP, size_t* colP, const size_t m, const size_t n, const size_t nnz, const bool ,const bool) : Father_t(m), _rowdim(m), _coldim(n), _nnz(nnz) {
            for(size_t k=0; k<nnz; ++k) 
                this->operator[](rowP[k]).push_back(colP[k]);
        }

            
        size_t rowdim() const { return _rowdim; }
        size_t coldim() const { return _coldim; }


	void setEntry(size_t i, size_t j, const Element& v) ;
	const Element& getEntry(size_t i, size_t j) const ;
	Element& getEntry(Element&, size_t i, size_t j) const ;

        template<class OutVector, class InVector>
        OutVector& apply(OutVector& y, const InVector& x) const; // y = A x
    
        template<class OutVector, class InVector>
        OutVector& applyTranspose(OutVector& y, const InVector& x) const; // y = A^T x    

            /** Read the matrix from a stream in ANY format
             *  entries are read as "long int" and set to 1 if they are odd,
             *  0 otherwise
             *  @param is Input stream from which to read the matrix
             *  @return Reference to input stream 
             */
        std::istream &read (std::istream &is) ;
	std::ostream& write (std::ostream& out, FileFormatTag format=FORMAT_GUILLAUME) const ;
   
        const Field& field() const { return _F; }

    template<typename _Tp1>
    struct rebind
    {
      typedef ZeroOne<_Tp1> other;
      void operator() (other & Ap,
                       const Self_t& A,
                       const _Tp1& F) {
              // ZeroOne does not store any Field element
      }
    };

        
        template<typename _Tp1>
        ZeroOne(ZeroOne<_Tp1>& A, const GF2 F2) : Father_t(A.rowdim()), _rowdim(A.rowdim()), _coldim(A.coldim()), _nnz(0) {
            for(typename ZeroOne<_Tp1>::RawIndexIterator it = A.indexBegin();
                it != A.indexEnd(); ++it,++_nnz) {
                this->operator[]( it->first ).push_back( it->second );
            }
        }

        size_t nnz() const { return _nnz; }
        bool isRowSorted() const { return true; }
        bool isColSorted() const { return true; }

    /** RawIterator class.  Iterates straight through the values of the matrix
     */
    class RawIterator;

    RawIterator rawBegin();
    RawIterator rawEnd();
    const RawIterator rawBegin() const;
    const RawIterator rawEnd() const;

     /** RawIndexIterator - Iterates through the i and j of the current element
   * and when accessed returns an STL pair containing the coordinates
   */
    class RawIndexIterator;
    RawIndexIterator indexBegin();
    const RawIndexIterator indexBegin() const;
    RawIndexIterator indexEnd();
    const RawIndexIterator indexEnd() const;


    private:
        size_t _rowdim, _coldim, _nnz;
  };
    
}

#include "linbox/blackbox/zo-gf2.inl"

#endif //__LINBOX_zo_gf2_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

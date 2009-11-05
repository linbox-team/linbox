/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/blackbox/zo-gf2.h
 * Copyright (C) 2009 The LinBox group
 *
 * Time-stamp: <05 Nov 09 10:34:40 Jean-Guillaume.Dumas@imag.fr> 
 *
 * See COPYING for license information.
 *
 */
#ifndef __ZO_GF2_H
#define __ZO_GF2_H

#include <algorithm>
#include "linbox/field/gf2.h"
#include "linbox/field/unparametric.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/vector/stream.h"
#include "linbox/matrix/sparse.h"
#include "linbox/blackbox/zero-one.h"
#include "linbox/vector/light_container.h"

namespace LinBox 
{
    
 /** \brief Time and space efficient representation of sparse matrices over GF2.
   * Representation if a full row array containing vector of non-zero locations
   * A 0-1 matrix is a matrix with all 0's and 1's as entries.  
   \ingroup blackbox
  */ 
    template<>
    struct ZeroOne<GF2> : LightContainer< LightContainer< size_t > > 
    {
        typedef LightContainer< LightContainer< size_t > > Father_t;
        typedef LightContainer< size_t > Row_t;
        typedef GF2::Element Element;
        typedef size_t Index;
        typedef ZeroOne<GF2> Self_t;
        typedef GF2 Field;

        const GF2 _F;

        ZeroOne(const GF2& F) {}
        ZeroOne(const GF2& F, const size_t m) : Father_t(m), _rowdim(m), _coldim(m) {}
        ZeroOne(const GF2& F, const size_t m, const size_t n) : Father_t(m), _rowdim(m), _coldim(n) {}

        ZeroOne() {}
        ZeroOne(const size_t m) : Father_t(m), _rowdim(m), _coldim(m) {}
        ZeroOne(const size_t m, const size_t n) : Father_t(m), _rowdim(m), _coldim(n) {}

        ZeroOne(const GF2& F, VectorStream<Row_t>& stream) : Father_t(stream.m()), _rowdim(stream.m()), _coldim(stream.n()) {
            for (Father_t::iterator row=begin(); row != end(); ++row)
                stream >> *row;
        }        

        ZeroOne(const Self_t& A) : Father_t(static_cast<const Father_t&>(A)), _rowdim(A._rowdim), _coldim(A._coldim) {}

        ZeroOne(const GF2& F, size_t* rowP, size_t* colP, const size_t m, const size_t n, const size_t nnz, const bool notused=true,const bool notused2=true) : Father_t(m), _rowdim(m), _coldim(n) {
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
      void operator() (other *& Ap,
                       const Self_t& A,
                       const _Tp1& F) {
	size_t nnz(0);
	for(Father_t::const_iterator ro=A.begin(); ro!= A.end(); ++ro)
		nnz+=ro->size();
	size_t * rowP = new size_t[nnz], * colP = new size_t[nnz];
	size_t cur=0;
	for(size_t i=0; i<A.rowdim(); ++i) {
		for( Row_t::const_iterator it=A[i].begin(); it != A[i].end() ; ++it, ++cur) {
			rowP[cur] = i;
			colP[cur] = *it;
		}
	}

        Ap = new other(F, rowP, colP, A.rowdim(), A.coldim(), nnz, true, true);
      }
    };


    private:
        size_t _rowdim, _coldim;
  };
    
}

#include "linbox/blackbox/zo-gf2.inl"
#endif

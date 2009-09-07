/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/blackbox/zo-gf2.h
 * Copyright (C) 2009 The LinBox group
 *
 * Time-stamp: <04 Sep 09 16:00:15 Jean-Guillaume.Dumas@imag.fr> 
 *
 * See COPYING for license information.
 *
 */
#ifndef __ZO_GF2_H
#define __ZO_GF2_H

#include "linbox/field/gf2.h"
#include "linbox/field/unparametric.h"
#include "linbox/util/matrix-stream.h"
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

        ZeroOne() {}
        ZeroOne(const GF2& F) {}
        ZeroOne(const size_t m) : Father_t(m), _rowdim(m), _coldim(m) {}
        ZeroOne(const size_t m, const size_t n) : Father_t(m), _rowdim(m), _coldim(n) {}
            
        size_t rowdim() const { return _rowdim; }
        size_t coldim() const { return _coldim; }

        template<class OutVector, class InVector>
        OutVector& apply(OutVector& y, const InVector& x) const; // y = Ax;
    
        template<class OutVector, class InVector>
        OutVector& applyTranspose(OutVector& y, const InVector& x) const; // y = ATx
    

            /** Read the matrix from a stream in ANY format
             *  entries are read as "long int" and set to 1 if they are odd,
             *  0 otherwise
             *  @param is Input stream from which to read the matrix
             *  @return Reference to input stream 
             */
        std::istream &read (std::istream &is) {
		// Reads a long int and take it mod 2 afterwards (v&1)
            UnparametricField<long> Ints;
            MatrixStream<UnparametricField<long> > S(Ints, is);
            S.getDimensions( _rowdim, _coldim );
            this->resize(_rowdim);
            Index r, c; 
            long v;
            while( S.nextTriple(r, c, v) ) {
                if (v&1) this->operator[](r).push_back(c);
            }
            return is;
        }
   
        const Field& field() const { return _F; }

    private:
        size_t _rowdim, _coldim;
  };
    
}

#include "linbox/blackbox/zo-gf2.inl"
#endif

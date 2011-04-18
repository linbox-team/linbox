/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 *    ntl-sylvester.h
 *    Copyright (C) 2003 Austin Lobo, B. David Saunders
 *
 *    Template for sylvester matrix specification for ntl Arithmetic,
 *    for polynomials in one variable.
 *    Linbox version 2003
 *    See COPYING for licence information
 *-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/

#ifndef __LINBOX_ntl_sylvester_H
#define __LINBOX_ntl_sylvester_H

#include <iostream>
#include <fstream>
#include <vector>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_p.h>

#include <linbox/blackbox/blackbox-interface.h>
#include "linbox/vector/vector-traits.h"

namespace LinBox 
{
    /// This is a representation of the Sylvester matrix of two polynomials.
/// \ingroup blackbox
    template <class _Field>
    class Sylvester : public  BlackboxInterface
    {
    public:
	  typedef _Field Field;
	  typedef typename Field::Element element;      // Currently restricted to ZZ_p
	  ~Sylvester();                                 // Destructor
	  Sylvester();                                  // Default Constructor
	  Sylvester( const Field F, 
			 const std::vector<element> &vpx,
			 const std::vector<element> &vpy ); // Constructor given 2 polys and Field

      
	  void   print( std::ostream& os = std::cout ) const; // Print to stream or stdout      
	  void   print( char *outFileName ) const;            // Print to file or stdout
	  void   printcp( char *outFileName) const;           // Print to file in sparse format

	  inline size_t rowdim() const { return rowDim; }
	  inline size_t coldim() const { return colDim; }
	  inline size_t sysdim() const { return sysDim; }
	  const Field& field() const { return K; }

        template<typename _Tp1>
        struct rebind
        { typedef Sylvester<_Tp1> other; };
        
	  template <class OutVector, class InVector>
	  OutVector& apply( OutVector &v_out, const InVector& v_in ) const;

	  template <class OutVector, class InVector>
	  OutVector& applyTranspose( OutVector &v_out, const InVector& v_in ) const;
      
	  //Sylvester(char *dataFileName ); // read from a file -- not implemented yet

    protected:
	  Field         K;

	  size_t        rowDim,
		colDim,
		sysDim;

	  NTL::ZZ_pX    pxdata,                // "Upper" Polynomial 
		qxdata;                // "Lower" Polynomial in Sylvester matrix

	  std::vector<NTL::ZZ_p>    pdata,     // Input vector of polynomial coeff
		qdata;     // Input vector of coeffs for second poly
      
	  size_t pxdeg() const { return pdata.size() - 1; }
	  size_t qxdeg() const { return qdata.size() - 1; }


      
    };// End,   template <class Field, class Vector>
}

#include <linbox/blackbox/ntl-sylvester.inl>     
    
#endif //__LINBOX_ntl_sylvester_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

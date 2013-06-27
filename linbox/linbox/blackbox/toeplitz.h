/* *******************************************************************
 *    ntl-toeplitz.h
 *    Copyright (C) 2002 Austin Lobo, B. David Saunders
 *
 *    Template for Toeplitz specification for ntl Arithmetic
 *    Linbox version 2001 and 2002 from a version
 *    Designed by A.Lobo and B.D. Saunders in 4/98
 *  ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/

#ifndef __LINBOX_toeplitz_H
#define __LINBOX_toeplitz_H

#include <iostream>
#include <vector>
#include "linbox/vector/vector-traits.h"
#include "linbox/solutions/methods.h"  // for shape
#include "linbox/solutions/solution-tags.h"  // to offer trace, det
#include "linbox-config.h"
#include "linbox/blackbox/blackbox-interface.h"

#ifdef __LINBOX_HAVE_NTL
#include "linbox/field/ntl.h"
#endif

//#define DBGMSGS 1

/*! @file blackbox/toeplitz.h
 *   @brief Partial Specialization of Toeplitz for Dense vectors from  an FFT based on Shoup's NTL library.
 *
 *   This file is included at the end of the template specification
 *   for Toeplitz.h
 */

namespace LinBox
{
	template <class _Field, class _PRing>
	class ToeplitzBase : public  BlackboxInterface {
	protected: // Constructors etc. are protected because stand-alone instances of this
		// class should not be constructed; use Toeplitz instead.
		typedef _PRing PRing;
		typedef typename _PRing::Element Poly;

		typedef typename _PRing::CoeffField Field;

		typedef typename _PRing::Coeff Element;

		//------- CONSTRUCTORS AND DESTRUCTORS

		virtual ~ToeplitzBase();                // Destructor
		// ToeplitzBase();                 // Zero Param Constructor
		ToeplitzBase( const Field& F);// Field only cstr. JGD 30.09.2003
		ToeplitzBase( const PRing& PF );
		// Constructor using a polynomial field, a polynomial in that
		// field representing the matrix, and the matrix dimensions
		// If no n is given, it will default to the same as m.
		ToeplitzBase( const PRing& PF, const Poly& p,
			      size_t m, size_t n=0 );


	public:
		//------- READ-ONLY ACCESSOR, and OBSERVER METHODS

		template<typename _Tp1, typename _Tp2>
		struct rebind
		{ typedef ToeplitzBase<_Tp1,_Tp2> other; };



		inline size_t rowdim() const { return rowDim; }// Number of Rows
		inline size_t coldim() const { return colDim; }// Number of Cols
		inline size_t sysdim() const { return sysDim; } // Max of rows & columns;
		const Field& field() const {return *field_;}

		//------- MUTATOR METHODS

		void setToUniModUT() ;      // Convert to UTriang matrix with det 1
		void setToUniModLT() ;      // Convert to LTriang matrix with det 1

		//------ SERVICE METHODS

		//      void convert(NTL::ZZ_pX &pout, const std::vector<Element> &vin);
	protected:
		PRing P; 	  // Polynomial field
		const Field *field_;           // Field parameter (coefficient field)

		size_t rowDim;             // row dimension
		size_t colDim;             // col dimension
		size_t sysDim;             // max{row,col} dimension

		Poly  pdata;         // Poly rep of the toeplitz matrix
		Poly  rpdata;        // reverse polynomial

		static const int UnimodUT=1;
		static const int UnimodLT=2;
		BlackboxSpecifier shape; // Helps us deduce what our shape is
		//std::vector<NTL::ZZ_p> data;    // The vector of coeffs of the polynomial

		/* These were only used by the XML stuff and are more or less
		 * obsolete with the presence of polynomial rings.

		 void convert(Poly &pout, const std::vector<Element> &vin);
		// CONVERTS the input vector of field Elements to a ZZ_pX

		void convert( std::vector<Element> &vout, const Poly &pin);
		// Converts from polynomial rep to a vector rep
		// inverse of the convert from Element to ZZ_pX
		*/

	}; //  class ToeplitzBase

	/** \brief This is the blackbox representation of a Toeplitz matrix.

	  \ingroup blackbox
	 * It stores the 2n-1 values of the first row and column.
	 * The _PRing template parameter should be a polynomial ring.
	 * Computations on the matrix will be performed using this polynomial ring.
	 * The apply is a call to polynomial multiplication and if the poly ring is
	 * FFT based for large n, apply will run in O(n lg(n)) time.
	 */
#ifdef __LINBOX_HAVE_NTL
	template< class _CField, class _PRing = NTL_ZZ_pX >
	#else
	template< class _CField, class _PRing >
	#endif
	class Toeplitz : public ToeplitzBase<_CField,_PRing> {
	protected:
		typedef ToeplitzBase<_CField,_PRing> TBase;
	public:
		typedef TBase Father_t;
		typedef typename TBase::PRing PRing;
		typedef typename TBase::Poly Poly;

		typedef typename TBase::Field Field;

		typedef typename TBase::Element Element;

		// Toeplitz() :
			// TBase()
		// {}
		Toeplitz(const Field& F) :
			TBase(F)
		{}
		Toeplitz( const PRing& PF ) :
			TBase(PF)
		{}
		Toeplitz( const PRing& PF, const Poly& p,
			  size_t m, size_t n=0 ) :
			TBase(PF,p,m,n)
		{}

		using TBase::field_;
		using TBase::field;
	}; // end of class Toeplitz

	/** Specialization for when the field of matrix elements is the same
	 * as the coefficient field of the polynomial field.
	 */
	template< class _PRing >
	class Toeplitz< typename _PRing::CoeffField, _PRing > : public ToeplitzBase< typename _PRing::CoeffField, _PRing > {
	protected:
		typedef ToeplitzBase< typename _PRing::CoeffField, _PRing > TBase;
		typedef TBase Father_t;
		using TBase::P;
		using TBase::rowDim;
		using TBase::colDim;
		using TBase::sysDim;
		using TBase::shape;
		using TBase::pdata;
		using TBase::rpdata;
		using TBase::field_;
	public:
		using TBase::rowdim;
		using TBase::coldim;
		typedef typename TBase::PRing PRing;
		typedef typename TBase::Poly Poly;

		typedef typename TBase::Field Field;

		typedef typename TBase::Element Element;

		//------- CONSTRUCTORS AND DESTRUCTORS
		// Toeplitz() :
			// TBase()
		// {}
		Toeplitz(const Field& F) :
			TBase(F)
		{}
		Toeplitz( const PRing& PF ) :
			TBase(PF)
		{}
		Toeplitz( const PRing& PF, const Poly& p,
			  size_t m, size_t n=0 ) :
			TBase(PF,p,m,n)
		{}

		Toeplitz( const Field& F,    // Cnstr. with Field and STL vec. of elems
			  const std::vector<Element>& v) :
			TBase(F)
		{ init_vector(v); }

		Toeplitz( const Field& F,    // Cnstr. with Field and STL vec. of elems
			  const BlasVector<Field>& v) :
			TBase(F)
		{ init_vector(v.getRep()); }

		void   write( std::ostream& os = std::cout) const;        // Print the contents to the screen

		void   write( char *outFileName) const;
		// Print the contents to a file

		template<class OutVector, class InVector>
		OutVector& apply( OutVector &v_out, const InVector& v_in) const;

		template<class OutVector, class InVector>
		OutVector& applyTranspose( OutVector &v_out, const InVector& v_in) const;

		// Get the determinant of the matrix
		Element& det( Element& res ) const;

		Element& trace( Element& res ) const;
	protected:
		// initialization via a vector. Usually called by a constructor
		// Moved in a separate protected function to enable easier
		// inherited constructor calls. JGD 30.09.2003
		void init_vector( const std::vector<Element>& v );
	}; //  Toeplitz specialization

	template<class Field, class PD>
	struct DetCategory<Toeplitz<Field,PD> >	{ typedef typename SolutionTags::Local Tag; };

	template<class Field, class PD>
	struct TraceCategory<Toeplitz<Field,PD> >	{ typedef typename SolutionTags::Local Tag; };


} // namespace LinBox

#include "linbox/blackbox/toeplitz.inl"

#endif //__LINBOX_toeplitz_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

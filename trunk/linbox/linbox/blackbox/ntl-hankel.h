/* *******************************************************************
 *    ntl-hankel.h
 * Copyright (C) 2003 Austin Lobo, B. David Saunders
 *    Template for Hankel specification for ntl Arithmetic
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

#ifndef __LINBOX_ntl_hankel_H
#define __LINBOX_ntl_hankel_H

#include "linbox/blackbox/blackbox-interface.h"
#include "toeplitz.h" // we inherit everything from ntl-toeplitz
// #include "linbox/vector/blas-vector.h"

//#define DBGMSGS 1

/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 *   Partial Specialization of Hankel for Dense vectors from
 *   an FFT based on Shoup's NTL library. Extends toeplitz matrix
 *-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/

namespace LinBox
{
	/// \ingroup blackbox
	template <class _Field>
	class Hankel: public Toeplitz<_Field> {
	public:
		typedef _Field Field;
		typedef typename Field::Element Element;

		typedef Toeplitz<_Field> TBase;
		typedef TBase Father_t;
                using TBase::P;
                using TBase::rowDim;
                using TBase::colDim;
                using TBase::sysDim;
                using TBase::shape;
                using TBase::pdata;
                using TBase::rpdata;
                using TBase::field_;
                using TBase::field;


		template<typename _Tp1>
		struct rebind
		{ typedef Hankel<_Tp1> other; };

		//------- CONSTRUCTORS AND DESTRUCTORS

		~Hankel();
		// Hankel();// : Toeplitz<_Field>(){}
		Hankel( const Field F) : Father_t (F) {};
		// Cnstr. with Field and STL vec. of elems
		Hankel( const Field F,    const std::vector<Element>&v);// : Toeplitz<_Field>(F, v){}
		Hankel( const BlasVector<Field>&v);// : Toeplitz<_Field>(F, v){}

		Hankel( const Field F,    size_t n);// : Toeplitz<_Field>(F, n) {}
		//	  Hankel(char *dataFileName ); // read from a file
		//void init( const Field F,    size_t n) { Toeplitz<_Field>::init(F, n) }

		//------- INHERITED READ-ONLY ACCESSOR, and OBSERVER METHODS

		void   print( std::ostream& os = std::cout) const; // Print to stdout
		void   print( char *outFileName) const;            // Print to file
		/*      inline size_t this->rowdim() const;// Number of Rows
		 *      inline size_t this->coldim() const;// Number of Cols
		 *      inline size_t sysdim() const;// Max of rows & columns;
		 */

		//------ SERVICE METHODS
		template<class OutVector, class InVector>
		OutVector& apply( OutVector &v_out, const InVector& v_in) const;

		template<class OutVector, class InVector>
		OutVector& applyTranspose( OutVector &v_out, const InVector& v_in) const;

	}; //  class Hankel

} // namespace Linbox

#include "linbox/blackbox/ntl-hankel.inl"

#endif //__LINBOX_ntl_hankel_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

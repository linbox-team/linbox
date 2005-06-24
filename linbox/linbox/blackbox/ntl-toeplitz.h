/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */

/* *******************************************************************
 *    ntl-toeplitz.h 
 *    Copyright (C) 2002 Austin Lobo, B. David Saunders
 *
 *    Template for Toeplitz specification for ntl Arithmetic
 *    Linbox version 2001 and 2002 from a version 
 *    Designed by A.Lobo and B.D. Saunders in 4/98
 *-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/

#ifndef NTL_TOEPLITZ_H
#define NTL_TOEPLITZ_H
#include <iostream>
#include <vector>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_p.h>
#include "linbox/vector/vector-traits.h"
#include "linbox-config.h"
#include <linbox/blackbox/blackbox-interface.h>

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

#include <algorithm>
#include <string>

#endif

#define SPARSE    0X1    // 0001
#define DENSE     0X2    // 0010
#define DIAGONAL  0X4    // 0100
#define TOEPLITZ  0X8    // 1000
#define HANKEL    0XC    // 1100
#define UNIMOD_UT 0XA    // 1010 -- unimodular upper triang. Toeplitz
#define UNIMOD_LT 0X9    // 1001 -- unimodular lower triang. Toeplitz
#define UNIMOD_UH 0XE    // 1110 -- unimodular upper triang. Hankel
#define UNIMOD_LH 0XD    // 1101 -- unimodular lower triang. Hankel
#define BLKVECTOR 0X10


//#define DBGMSGS 1

/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 *   Partial Specialization of Toeplitz for Dense vectors from 
 *   an FFT based on Shoup's NTL library.
 *   this file is included at the end of the template specification
 *   for Toeplitz.h
 *-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/

namespace LinBox
{
	/** \brief This is the blackbox representation of a Toeplitz matrix.

\ingroup blackbox
	 * It stores the 2n-1 values of the first row and column.
	 * The apply is a call to polynomial multiplication and for large n
	 * will be FFT based, running in O(lg(n)) time.
	 */
	template <class _Field>
	class Toeplitz : public  BlackboxInterface
	{
	public:
		typedef _Field Field;

		typedef typename Field::Element Element;
		
		//------- CONSTRUCTORS AND DESTRUCTORS
		
		~Toeplitz();                // Destructor
		Toeplitz();                 // Zero Param Constructor
                Toeplitz( const Field& F) : K(F) {} // Field only cstr. JGD 30.09.2003
            
		Toeplitz( const Field& F,    // Cnstr. with Field and STL vec. of elems
				  const std::vector<Element>&v);

		
		//------- READ-ONLY ACCESSOR, and OBSERVER METHODS 

#ifndef __LINBOX_XMLENABLED		
		void   print( std::ostream& os = std::cout) const;        // Print the contents to the screen
		void   print( char *outFileName) const; 
#else
		Toeplitz(LinBox::Reader &);
		Toeplitz(const Toeplitz<Field, Vector>&);

		std::ostream &write(std::ostream &) const;
		bool toTag(LinBox::Writer &) const;
#endif
            template<typename _Tp1>
            struct rebind
            { typedef Toeplitz<_Tp1> other; };


		
		inline size_t rowdim() const;// Number of Rows
		inline size_t coldim() const;// Number of Cols
		inline size_t sysdim() const;// Max of rows & columns; 
		const Field& field() const {return K;}
		// Print the contents to a file
		//------- MUTATOR METHODS
		
		void setToUniModUT() ;      // Convert to UTriang matrix with det 1
		void setToUniModLT() ;      // Convert to LTriang matrix with det 1
		
		//------ SERVICE METHODS

		template<class OutVector, class InVector>
		OutVector& apply( OutVector &v_out, const InVector& v_in) const;

		template<class OutVector, class InVector>
		OutVector& applyTranspose( OutVector &v_out, const InVector& v_in) const;
		//      void convert(NTL::ZZ_pX &pout, const std::vector<Element> &vin);
    protected:
		Field K;                   // Field parameter
		
		size_t rowDim;             // row dimension 
		size_t colDim;             // col dimension 
		size_t sysDim;             // max{row,col} dimension
		
		NTL::ZZ_pX  pdata;         // Poly rep of the toeplitz matrix
		NTL::ZZ_pX  rpdata;        // reverse polynomial
		
		static const int UnimodUT=1;
		static const int UnimodLT=2;
		int shape;                // Helps us deduce what our shape is
		std::vector<NTL::ZZ_p> data;    // The vector of coeffs of the polynomial


                // initialization via a vector. Usually called by a constructor
                // Moved in a separate protected function to enable easier
                // inherited constructor calls. JGD 30.09.2003
                void init_vector( const std::vector<typename Field::Element>&v) ;


		void convert(NTL::ZZ_pX &pout, const std::vector<Element> &vin);
		// CONVERTS the input vector of field Elements to a ZZ_pX
		// use the convert for the field Element to integer and use
		// this integer to initialize a ZZ_p Element which will be the coeff
		// of a ZZ_pX poly
		
		void convert( const std::vector<Element> &vout, 
					  class NTL::ZZ_pX &pin);
		// Converts from polynomial rep to a vector rep
		// inverse of the convert from Element to ZZ_pX
		
	}; //  class Toeplitz
	
} // namespace Linbox

#include <linbox/blackbox/ntl-toeplitz.inl>     
// Hide the implementation; include it here because
// older compilers want everything in one template file
    
#endif

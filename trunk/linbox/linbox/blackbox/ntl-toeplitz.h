/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */

/* *******************************************************************
 *    ntl-toeplitz.h 
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
#include "linbox/blackbox/archetype.h"
#include "linbox/vector/vector-traits.h"

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
	template <class Field, class Vector>
	class Toeplitz: public BlackboxArchetype<Vector>
	{
	public:
		typedef typename Field::Element element;
		
		//------- CONSTRUCTORS AND DESTRUCTORS
		
		~Toeplitz();                // Destructor
		Toeplitz();                 // Zero Param Constructor
		Toeplitz( const Field F,    // Cnstr. with Field and STL vec. of elems
				  const std::vector<element>&v);
		//	  Toeplitz(char *dataFileName ); // read from a file
		BlackboxArchetype<Vector>* clone() const;
		
		//------- READ-ONLY ACCESSOR, and OBSERVER METHODS 
		
		void   print( std::ostream& os = std::cout) const;        // Print the contents to the screen
		void   print( char *outFileName) const; 
		
		inline size_t rowdim() const;// Number of Rows
		inline size_t coldim() const;// Number of Cols
		inline size_t sysdim() const;// Max of rows & columns; 
		
		// Print the contents to a file
		//------- MUTATOR METHODS
		
		void setToUniModUT() ;      // Convert to UTriang matrix with det 1
		void setToUniModLT() ;      // Convert to LTriang matrix with det 1
		
		//------ SERVICE METHODS
		Vector& apply( Vector &v_out, const Vector& v_in) const;
		Vector& applyTranspose( Vector &v_out, const Vector& v_in) const;
		//      void convert(NTL::ZZ_pX &pout, const std::vector<element> &vin);
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
		void convert(NTL::ZZ_pX &pout, const std::vector<element> &vin);
		// CONVERTS the input vector of field elements to a ZZ_pX
		// use the convert for the field element to integer and use
		// this integer to initialize a ZZ_p element which will be the coeff
		// of a ZZ_pX poly
		
		void convert( const std::vector<element> &vout, 
					  class NTL::ZZ_pX &pin);
		// Converts from polynomial rep to a vector rep
		// inverse of the convert from element to ZZ_pX
		
	}; //  class Toeplitz
	
} // namespace Linbox

#include <linbox/blackbox/ntl-toeplitz.inl>     
// Hide the implementation; include it here because
// older compilers want everything in one template file
    
#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */

/* *******************************************************************
 *    ntl-hankel.h 
 * Copyright (C) 2003 Austin Lobo, B. David Saunders
 *    Template for Hankel specification for ntl Arithmetic
 *    Linbox version 2001 and 2002 from a version 
 *    Designed by A.Lobo and B.D. Saunders in 4/98
 *-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/

#ifndef NTL_HANKEL_H
#define NTL_HANKEL_H

#include "ntl-toeplitz.h" // we inherit everything from ntl-toeplitz

//#define DBGMSGS 1

/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 *   Partial Specialization of Hankel for Dense vectors from 
 *   an FFT based on Shoup's NTL library. Extends toeplitz matrix
 *-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/

namespace LinBox
{
  template <class Field, class Vector>
    class Hankel: public Toeplitz<Field,Vector>
    {
    public:
      typedef typename Field::Element element;
      
      //------- CONSTRUCTORS AND DESTRUCTORS
      
      ~Hankel();                // Destructor
      Hankel();                 // Zero Param Constructor
      Hankel( const Field F,    // Cnstr. with Field and STL vec. of elems
	      const std::vector<element>&v);
      //	  Hankel(char *dataFileName ); // read from a file
      BlackboxArchetype<Vector>* clone() const;
      
      //------- INHERITED READ-ONLY ACCESSOR, and OBSERVER METHODS 
      
		void   print( std::ostream& os = std::cout) const; // Print to stdout
		void   print( char *outFileName) const;            // Print to file
       /*      inline size_t rowdim() const;// Number of Rows
		*      inline size_t coldim() const;// Number of Cols
		*      inline size_t sysdim() const;// Max of rows & columns; 
		*/

      //------- MUTATOR METHODS
      
      void setToUniModUT() ;      // Convert to UTriang matrix with det 1
      void setToUniModLT() ;      // Convert to LTriang matrix with det 1
      
      //------ SERVICE METHODS
      Vector& apply( Vector &v_out, const Vector& v_in) const;
      Vector& applyTranspose( Vector &v_out, const Vector& v_in) const;
      
    }; //  class Hankel
  
} // namespace Linbox

#include <linbox/blackbox/ntl-hankel.inl>     

#endif

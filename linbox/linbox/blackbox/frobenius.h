/* linobx/blackbox/frobenuis.h
 *
 * Written by Austin Lobo <alobo@cis.udel.edu> and
 *            B.D. Saunders <saunders@cis.udel.edu>
 * See COPYING
 */
#ifndef __FROBENIUS_H
#define __FROBENIUS_H

#include "linbox/blackbox/companion.h"
#include "linbox/blackbox/direct-sum.h"
#include <vector>

namespace LinBox {
  template <class Field, class Vector>
  class Frobenius: public DirectSum<Field,Vector>
  {
  public:
    Frobenius() { ; }   // default constructor

    /**
     *    Build a matrix in Frobenius form whose block sizes are
     *    specified by vlist, generated from random polynomials 
     *    @param vlist diagonal-block sizes, positive ints in non-increasing order
     */
    template <class VDegList>
      Frobenius( const Field &F, const VDegList &vlist)
      {
	; 
      }
    
    /**
     *    Build a square, block-diagonal matrix as a direct sum of the companion
     *    matrices of the polynomials. The dimension is the sum of the degrees.
     *    @param pbegin iterator pointing to the start of a list of polynomials
     *    @param pend   iterator pointing after end   of a list of polynomials
     */
    template <class PolyIterator>
      Frobenius( const Field &F, PolyIterator pbegin, PolyIterator pend)
      { 
	
	if ( pbegin != pend) {
	  _Ap = new Companion<Field, Vector>( F, *pbegin);
	  ++pbegin;       // go to next polynomial in list
	  _Bp = new Frobenius( F, pbegin, pend);
	} 
      }
    

    BlackboxArchetype<Vector>* clone () const
      { 
	return new Frobenius(*this); 
      }
    
  private:
  }; // class Frobenius
  
}// Namespace LinBox
#endif

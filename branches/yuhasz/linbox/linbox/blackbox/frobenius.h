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
  template <class _Field>
  class Frobenius: public DirectSum<Companion<_Field> >
  {
  public:
    Frobenius() { }   // default constructor

    /**
     *    Build a matrix in Frobenius form whose block sizes are
     *    specified by vlist, generated from random polynomials 
     *    @param vlist diagonal-block sizes, positive ints in non-increasing order
     */
    template <class VDegList>
      Frobenius( const _Field &F, const VDegList &vlist)
      {
      }
    
    /**
     *    Build a square, block-diagonal matrix as a direct sum of the companion
     *    matrices of the polynomials. The dimension is the sum of the degrees.
     *    @param pbegin iterator pointing to the start of a list of polynomials
     *    @param pend   iterator pointing after end   of a list of polynomials
     */
    template <class PolyIterator>
      Frobenius( const _Field &F, PolyIterator pbegin, PolyIterator pend) {
		_VB.resize(pend - pbegin);
		PolyIterator pp = pbegin;
		typename std::vector<const Companion<_Field>* >::iterator vp;
		m = 0;
		n = 0;
		for(vp = _VB.begin(); vp != _VB.end(); ++vp,++pp)  {
			*vp = new  Companion<_Field>(F,*pp);
			m += (*vp) -> rowdim();
			n += (*vp) -> coldim();
		}
	}


	~Frobenius() {
		typename std::vector< const Companion<_Field>* >::iterator vp;
		for(vp = _VB.begin(); vp != _VB.end(); ++vp)
			delete (*vp);
	}

  }; // class Frobenius
  
}// Namespace LinBox
#endif

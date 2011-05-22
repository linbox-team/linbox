/* linobx/blackbox/frobenuis.h
 * Copyright(c) 2010 LinBox
 * Written by Austin Lobo <alobo@cis.udel.edu> and
 *            B.D. Saunders <saunders@cis.udel.edu>
 * See COPYING for licence information
 */
#ifndef __LINBOX_frobenius_H
#define __LINBOX_frobenius_H

#include <linbox/blackbox/blackbox-interface.h>
#include "linbox/blackbox/companion.h"
#include "linbox/blackbox/direct-sum.h"
#include <vector>

namespace LinBox 
{
  /// \ingroup blackbox
  template <class _Field>
  class Frobenius: public BlackboxInterface, public DirectSum<Companion<_Field> >
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
		this->_VB.resize(pend - pbegin);
		PolyIterator pp = pbegin;
		typename std::vector<const Companion<_Field>* >::iterator vp;
		this->m = 0;
		this->n = 0;
		for(vp = this->_VB.begin(); vp != this->_VB.end(); ++vp,++pp)  {
			*vp = new  Companion<_Field>(F,*pp);
			this->m += (*vp) -> rowdim();
			this->n += (*vp) -> coldim();
		}
	}


	~Frobenius() {
		typename std::vector< const Companion<_Field>* >::iterator vp;
		for(vp = this->_VB.begin(); vp != this->_VB.end(); ++vp)
			delete (*vp);
	}


      template<typename _Tp1>
      struct rebind
      { typedef Frobenius<_Tp1> other; };



  }; // class Frobenius
  
}// Namespace LinBox

#endif //__LINBOX_frobenius_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

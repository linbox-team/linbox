/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/dense.inl
 * Copyright (C) 2001 B. David Saunders, 
 *               2001-2002 Bradford Hovinen, 
 *               2002 Zhendong Wan
 *
 * Written by B. David Saunders <saunders@cis.udel.edu>, 
 *            Bradford Hovinen <hovinen@cis.udel.edu>, 
 *            Zhendong Wan <wan@mail.eecis.udel.edu>
 *
 * evolved from dense-matrix.h by -bds, Zhendong Wan
 *
 * --------------------------------------------------------
 * 2002-10-27  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Split out container/iterator functionality into DenseMatrixBase
 * --------------------------------------------------------
 * 2002-08-09  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Renamed file from dense-matrix1.C to dense.inl
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#ifndef __BLACKBOX_DENSE_INL
#define __BLACKBOX_DENSE_INL

#include <iostream>
#include <vector>
#include <fstream>

#include "linbox/blackbox/dense.h"
#include "linbox/util/debug.h"

namespace LinBox
{

template <class Field>
template<class Vect1, class Vect2>
Vect1& DenseMatrix<Field>::apply (Vect1& y, const Vect2& x) const
{
	typename DenseMatrixBase<Element>::ConstRowIterator p;
	typename Vect1::iterator p_y = y.begin ();  

	for (p = rowBegin (); p != rowEnd (); ++p, ++p_y)
		_VD.dot (*p_y, *p, x);
    
	return y;
}
 
template <class Field>
template<class Iterator1, class Iterator2>
Iterator1& DenseMatrix<Field>::apply (Iterator1        in,
					      const Iterator2 &outbegin,
					      const Iterator2 &outend) const
{
	linbox_check (coldim () == (outend - outbegin));

	typename DenseMatrixBase<Element>::ConstRowIterator rowp;
	Iterator2 p_out;
	typename DenseMatrixBase<Element>::Row::const_iterator pe;

	for (rowp = rowBegin (); rowp != rowEnd (); ++rowp, ++in) {
		_F.init (*in, 0);
		for (pe = rowp->begin (), p_out = outbegin; pe != rowp->end (); ++pe, ++p_out)
			_F.axpyin (*in, *pe, *p_out);
	}
    
	return in;
}
 
template <class Field>
template<class Vect1, class Vect2>
Vect1& DenseMatrix<Field>::applyTranspose (Vect1& y, const Vect2& x) const
{
	typename DenseMatrixBase<Element>::ConstColIterator colp;
	typename Vect1::iterator p_y = y.begin ();  

	for (colp = colBegin (); colp != colEnd (); ++colp, ++p_y)
		_VD.dot (*p_y, *colp, x);
    
	return y;
}
  
template <class Field>
template<class Iterator1, class Iterator2>
Iterator1& DenseMatrix<Field>::applyTranspose (Iterator1        in,
						       const Iterator2 &outbegin,
						       const Iterator2 &outend) const
{
	linbox_check (rowdim () == (outend - outbegin));

	typename DenseMatrixBase<Element>::ConstColIterator colp;
	Iterator2 p_out;
	typename DenseMatrixBase<Element>::Col::const_iterator pe;

	for (colp = colBegin (); colp != colEnd (); ++colp, ++in) {
		_F.init (*in, 0);
		for (pe = colp->begin (), p_out = outbegin; pe != colp->end (); ++pe, ++p_out)
			_F.axpyin (*in, *pe, *p_out);
	}
    
	return in;
}

} // namespace LinBox

#endif // __BLACKBOX_DENSE_INL

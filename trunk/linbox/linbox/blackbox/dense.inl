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
Vect1& DenseMatrix<Field>::apply (Vect1& y, const Vect2& x) const {

#ifdef __LINBOX_PARALLEL

	return BlackboxParallel (y, *this, x, BBBase::Apply);
#else

	_MD. vectorMul (y, *this, x);

#endif
	return y;
}

 
template <class Field>
template<class Vect1, class Vect2>
Vect1& DenseMatrix<Field>::applyTranspose (Vect1& y, const Vect2& x) const {

#ifdef __LINBOX_PARALLEL

        return BlackboxParallel (y, *this, x, BBBase::ApplyTranspose);
#else

	return _MD.vectorMul (y, TransposeMatrix<DenseMatrix<Field> > (*this), x);
#endif
	    
}
  
} // namespace LinBox

#endif // __BLACKBOX_DENSE_INL

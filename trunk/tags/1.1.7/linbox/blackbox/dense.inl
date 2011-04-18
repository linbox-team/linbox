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

#ifndef __LINBOX_blackbox_dense_INL
#define __LINBOX_blackbox_dense_INL

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

	return _MD.vectorMul (y, _AT, x);
#endif
	    
}

template< class Field, class BElement >
DenseMatrix<Field>*
	DenseMatrixFactory<Field,BElement>::makeBlackbox( const Field& F )
{
	DenseMatrixBase<typename Field::Element> newBase ( rowdim(), coldim() );
	
	typename DenseMatrixBase<BElement>::ConstRawIterator i;
	typename DenseMatrixBase<typename Field::Element>::RawIterator j;

	for( i = _A.rawBegin(), j = newBase.rawBegin();
	     i != _A.rawEnd(), j != newBase.rawEnd();
	     ++i, ++j )
		F.init( *j, *i );

	return new DenseMatrix<Field>( F, newBase );
}

template< class Field, class BElement >
integer& DenseMatrixFactory<Field,BElement>::maxNorm( integer& res ) {
	typename DenseMatrixBase<BElement>::ConstRawIterator i;
	res = 0L;
	integer tmp;
	
	for( i = _A.rawBegin(); i != _A.rawEnd(); ++i ) {
		tmp = abs( *i );
		if( res < tmp ) res = tmp;
	}

	return res;
}

template< class Field, class BElement >
integer& DenseMatrixFactory<Field,BElement>::hadamardBound(integer& res) const {
	typename DenseMatrixBase<BElement>::ConstRowIterator r;
	typename DenseMatrixBase<BElement>::ConstRow::const_iterator c;

	res = 1L;
	integer temp;
	
	for( r = _A.rowBegin(); r != _A.rowEnd(); ++r ) {
		temp = 0;
		for( c = r->begin(); c != r->end(); ++c )
			temp += static_cast<integer>((*c)) * (*c);
		res *= temp;
	}

	res = sqrt(res);
	return res;
}
  
} // namespace LinBox

#endif // __LINBOX_blackbox_dense_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

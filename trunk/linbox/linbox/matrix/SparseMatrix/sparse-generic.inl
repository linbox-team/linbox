/* linbox/matrix/sparse.inl
 * Copyright (C) 2001-2002 Bradford Hovinen
 *               1999-2001 William J Turner,
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 * Based on sparse-base.h by William J Turner <wjturner@math.ncsu.edu>
 *
 * --------------------------------------------------------
 * 2003-01-11  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Move from blackbox/sparse-base.inl to matrix/sparse.inl
 * ------------------------------------
 * 2002-11-28  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 *   - Renamed ColOfRowsIterator to RowIterator
 *   - Named template argument _Row rather than Row; add a typedef to Row
 * ------------------------------------
 *
 *
 * ========LICENCE========
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
 *.
 */

#ifndef __LINBOX_matrix_sparse_INL
#define __LINBOX_matrix_sparse_INL

#include "linbox/linbox-config.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstring>

#include "linbox/matrix/sparse.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/util/debug.h"





namespace LinBox { namespace Protected {
	template <class Field, class Row, class Tag>
	SparseMatrixGeneric<Field,Row,Tag> ::SparseMatrixGeneric( MatrixStream<Field>& ms ) :
		_matA(0), _m(0), _n(0)
	{
		Element val;
		size_t i, j;
		while( ms.nextTriple(i,j,val) ) {
			if( i >= _m ) {
				_m = i + 1;
				_matA.resize( _m );
			}
			if( j >= _n ) _n = j + 1;
			setEntry(i,j,val);
		}
		if( ms.getError() > END_OF_MATRIX )
			throw ms.reportError(__func__,__LINE__);
		if( !ms.getDimensions( i, _n ) )
			throw ms.reportError(__func__,__LINE__);
		if( i > _m ) {
			_m = i;
			_matA.resize(_m);
		}
	}



#if 0 /*  weird */
	template <class Field,  class _Row, class BRow>
	SparseMatrixGeneric<Field,_Row> *SparseMatrixFactory<Field, _Row, BRow>::makeBlackbox (const Field &F)
	{
		SparseMatrixGeneric<Field, _Row> *A = new SparseMatrixGeneric<Field, _Row> (F, rowdim (), coldim ());

		typename SparseMatrixGeneric<Field, BRow>::ConstIterator i;
		typename SparseMatrixGeneric<Field, BRow>::ConstIndexedIterator j;

		for (i = _matA.Begin (), j = _matA.IndexedBegin (); i != _matA.End (); ++i, ++j)
			A.field().init (A->refEntry (j.rowIndex (), j.colIndex ()), *i);

		return A;
	}
#endif

} // namespace LinBox
} // namespace Protected

#endif // __LINBOX_matrix_sparse_INL


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/blackbox/sparse.inl
 * Copyright (C) 1999-2001 William J Turner,
 *               2001-2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * Modified by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Refactoring:
 *   - Eliminated SparseMatrixAux and moved that functionality into Sparse0
 *   - Made SparseMatrixBase parameterized only on the element type
 *   - New read/write implementations for SparseMatrixBase, supporting multiple
 *     formats
 *   - Eliminated Gaussian elimination code
 *   - Added iterators, including ColOfRowsIterator, Iterator, and
 *     IndexIterator
 *   - Eliminated operator []; added getEntry; changed put_value to setEntry
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

#ifndef __LINBOX_blackbox_sparse_INL
#define __LINBOX_blackbox_sparse_INL

#include "linbox/blackbox/sparse.h"

namespace LinBox
{

	template <class Field, class BElement,  class _Row, class BRow>
	SparseMatrix<Field,_Row> *SparseMatrixFactory<Field, BElement, _Row, BRow>::makeBlackbox (const Field &F)
	{
		SparseMatrix<Field, _Row> *A = new SparseMatrix<Field, _Row> (F, rowdim (), coldim ());

		typename SparseMatrixBase<BElement, BRow>::ConstIterator i;
		typename SparseMatrixBase<BElement, BRow>::ConstIndexedIterator j;

		for (i = _matA.Begin (), j = _matA.IndexedBegin (); i != _matA.End (); ++i, ++j)
			F.init (A->refEntry (j.rowIndex (), j.colIndex ()), *i);

		return A;
	}

}

#endif // __LINBOX_blackbox_sparse_INL


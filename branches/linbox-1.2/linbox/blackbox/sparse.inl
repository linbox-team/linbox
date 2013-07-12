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
 * See COPYING for license information.
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

		for (i = _A.Begin (), j = _A.IndexedBegin (); i != _A.End (); ++i, ++j)
			F.init (A->refEntry (j.rowIndex (), j.colIndex ()), *i);

		return A;
	}

}

#endif // __LINBOX_blackbox_sparse_INL

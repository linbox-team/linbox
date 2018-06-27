/* linbox/algorithms/blockbb.h
 * Copyright (C) 2018 Gavin Harrison
 *
 * Written by Gavin Harrison <gavin.har@gmail.com>
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
 */


#ifndef __LINBOX_blockbb_H
#define __LINBOX_blockbb_H

#include <iostream>
#include "linbox/util/error.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
namespace LinBox {
	
template<class _BB> 
struct is_blockbb { 
	static const bool value = false;
};

/// converts a black box into a block black box
template<class _BB>
class BlockBB 
{
public:
	typedef _BB BB;
	typedef typename BB::Field Field;
	
protected:
	BB _bb;
	
public:
	BlockBB(BB &bb) : _bb(bb) {}
	BlockBB(BlockBB &bb) : _bb(bb._bb) {}

	size_t rowdim() const {
		return _bb.rowdim();
	}
	
	size_t coldim() const {
		return _bb.coldim();
	}
	
	const Field& field() const {
		return _bb.field();
	}
	
	// Y = A*X
	template<class Matrix>
	Matrix& applyLeft(Matrix& Y, const Matrix& X) const {
		typename Matrix::ColIterator p1 = Y.colBegin();
		typename Matrix::ConstColIterator p2 = X.colBegin();
		
		for (; p2 != X.colEnd(); ++p1, ++p2) {
			_bb.apply(*p1, *p2);
		}
		
		return Y;
	}
	
	// Y = X*A
	template<class Matrix>
	Matrix& applyRight(Matrix& Y, const Matrix& X) const {
		typename Matrix::ColIterator p1 = Y.rowBegin();
		typename Matrix::ConstColIterator p2 = X.rowBegin();
		
		for (; p2 != X.rowEnd(); ++p1, ++p2) {
			_bb.applyTranspose(*p1, *p2);
		}
		
		return Y;
	}
	
	template<class OutVector, class InVector>
	OutVector& apply(OutVector& y, const InVector& x) const {
		return _bb.apply(y, x);
	}
	
	template<class OutVector, class InVector>
	OutVector& applyTranspose(OutVector& y, const InVector& x) const {
		return _bb.applyTranspose(y, x);
	}
	
}; // class BlockBB

template<class _BB>
struct is_blockbb<BlockBB<_BB>> {
	static const bool value = true;
};

} // LinBox
#endif // __LINBOX_blockbb_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

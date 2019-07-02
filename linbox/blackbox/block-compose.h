/* linbox/blackbox/compose.h
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_block_compose_H
#define __LINBOX_block_compose_H

#include "linbox/util/debug.h"
#include "linbox/linbox-config.h"
#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/vector/blas-vector.h"

#include "linbox/blackbox/blockbb.h"

namespace LinBox
{

/**
 * Blackbox of a product: \f$C = AB\f$, i.e \f$Cx \gets A(Bx)\f$.

 * This is a class that multiplies two matrices by implementing an
 * apply method that calls the apply methods of both of the consituent
 * matrices, one after the other.
 *
 * This class, like the Black Box archetype from which it is derived,
 * is templatized by the vector type to which the matrix is applied.
 * Both constituent matrices must also use this same vector type.
 * For specification of the blackbox members see \ref BlackboxArchetype.
 *
 * <b> Template parameter:</b> must meet the \ref Vector requirement.
 \ingroup blackbox
 */
template <class _Blackbox1, class _Blackbox2>
class BlockCompose : public BlackboxInterface {
	typedef BlockCompose<_Blackbox1, _Blackbox2> Self_t;
public:

	typedef _Blackbox1 Blackbox1;
	typedef _Blackbox2 Blackbox2;

	typedef typename Blackbox2::Field Field;
	typedef typename Field::Element Element;
	
protected:
	Blackbox1 _A;
	Blackbox2 _B;
	
public:

	/** Constructor of C := A*B from blackbox matrices A and B.
	 * Build the product A*B of any two black box matrices of compatible dimensions.
	 * @pre <code>A.coldim() == B.rowdim()</code>.
	 * @param A blackbox
	 * @param B blackbox
	 */
	BlockCompose(const Blackbox1 &A, const Blackbox2 &B) :
		_A(A), _B(B)
	{
	}
	
	size_t rowdim() const {
		return _A.rowdim();
	}
	
	size_t coldim() const {
		return _B.coldim();
	}
	
	const Field& field() const {
		return _A.field();
	}

	// Y = A*X
	template<class Matrix>
	Matrix& applyLeft(Matrix& Y, const Matrix& X) const {
		BlasMatrix<Field> Z(_A.field(), _A.rowdim(), X.coldim());
		
		_B.applyLeft(Z, X);
		_A.applyLeft(Y, Z);
		
		return Y;
	}
	
	// Y = X*A
	template<class Matrix>
	Matrix& applyRight(Matrix& Y, const Matrix& X) const {
		BlasMatrix<Field> Z(_A.field(), X.rowdim(), _B.coldim());
		
		_A.applyRight(Z, X);
		_B.applyRight(Y, Z);
		
		return Y;
	}
};

template<class A, class B>
struct is_blockbb<BlockCompose<A, B>> {
	static const bool value = true;
};
} // LinBox

#endif // __LINBOX_block_compose_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

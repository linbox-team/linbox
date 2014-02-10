/*
 * Copyright (C) 2014 the LinBox
 *
 * Written by :
 * BB <bbboyer@ncsu.edu>
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

/*! @file matrix/matrix-domain.h
 * @ingroup matrixdomain
 * @brief NODOC
 */

#ifndef __LINBOX_matrix_matrix_domain_H
#define __LINBOX_matrix_matrix_domain_H

#include "linbox/vector/vector-domain.h"
// #include "linbox/matrix/DenseMatrix/blas-matrix.h"

namespace LinBox {
#if 1
	template<class Field, class Rep> class BlasMatrix;
	template<class Matrix> class BlasSubmatrix;
	template<class _Field, class _Rep> class BlasVector ;
	template<class Field, class _Rep> class TriangularBlasMatrix ;
	template<class Matrix> class TransposedBlasMatrix ;
#endif


	/** \brief Helper class to allow specializations of certain matrix-vector products
	 *
	 * This class implements a method mulColSPD that multiplies a
	 * column-represented matrix by a dense vector
	 */
	template <class Field>
	class MVProductDomain {
	public:
		typedef typename Field::Element Element;

		MVProductDomain () {}

	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense (const VectorDomain<Field> &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const;
	};
}

#include "MatrixDomain/matrix-domain.h" // needed by BlasMatrix
#include "MatrixDomain/blas-matrix-domain.h"
#include "MatrixDomain/matrix-domain-gf2.h"
// #include "MatrixDomain/matrix-domain-m4ri.h"
#include "MatrixDomain/plain-domain.h"

#include "MatrixDomain/opencl-domain.h"


#endif // __LINBOX_matrix_matrix_domain_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

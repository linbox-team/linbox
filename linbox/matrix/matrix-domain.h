/*
 * Copyright (C) 2014 the LinBox
 *
 * Written by :
 * Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

#include <linbox/linbox-config.h>
#include "linbox/vector/vector-domain.h"
// #include "linbox/matrix/densematrix/blas-matrix.h"

namespace LinBox {
#if 1
	//! @bug necessitates the inclusion of those things after of before. We should not include directly this but just the matrices and they automatically need their domain (for now)
	template<class Field, class Rep> class BlasMatrix;
	template<class Matrix> class BlasSubmatrix;
	template<class _Field, class _Rep> class BlasVector ;
	//template<class Field, class _Rep> class TriangularBlasMatrix ;
    template<class Matrix> class TriangularBlasMatrix ;
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

#include "linbox/matrix/matrixdomain/matrix-domain.h" // needed by BlasMatrix
#include "linbox/matrix/matrixdomain/blas-matrix-domain.h"
#include "linbox/matrix/matrixdomain/matrix-domain-gf2.h"
// #include "matrixdomain/matrix-domain-m4ri.h"
#include "linbox/matrix/matrixdomain/plain-domain.h"

#include "linbox/matrix/matrixdomain/opencl-domain.h"


#endif // __LINBOX_matrix_matrix_domain_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

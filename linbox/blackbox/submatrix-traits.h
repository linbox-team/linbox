/* -*- mode:C++ -*- */

/* File: submatrix-traits.h
 *  Author: Zhendong Wan
 */
#ifndef __SUBMATRIX_TRAITS_H__
#define __SUBMATRIX_TRAITS_H__
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/submatrix.h>

namespace LinBox {

	template<class Matrix>
	class SubMatrixTraits;

	template<class Field>
	class SubMatrixTraits<DenseMatrix<Field> > {
	
	public:

		typedef  Submatrix<DenseMatrix<Field> > value_type;
	};


	template<class Field>
	class SubMatrixTraits<Submatrix<DenseMatrix<Field> > > {
	
	public:
	
		typedef Submatrix<DenseMatrix<Field> > value_type;
	};

}

#endif

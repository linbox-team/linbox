#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixAddSub_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixAddSub_H

#include <vector>
#include "linbox/matrix/matrixdomain/blas-matrix-domain.h"
#include "linbox/matrix/slicedpolynomialmatrix/SlicedPolynomialMatrix.h"

namespace LinBox
{
	/* internal
	 * Adding two matrices
	 */
	template< class Field, class Operand1, class Operand2, class Operand3>
	class SlicedPolynomialMatrixAdd
	{
	public:
		Operand1 &operator() (const Field &F, Operand1 &C, const Operand2 &A, const Operand3 &B) const;

	};

	/* internal
	 * Substracting two matrices
	 */
	template< class Field, class Operand1, class Operand2, class Operand3>
	class SlicedPolynomialMatrixSub
	{
	public:
		Operand1 &operator() (const Field &F, Operand1 &C, const Operand2 &A, const Operand3 &B) const;

	};

	//! C += A
	template< class Field, class Operand1, class Operand2>
	class SlicedPolynomialMatrixAddin
	{
	public:
		Operand1 &operator() (const Field &F, Operand1 &C, const Operand2 &A) const;

	};

	//! C -= A
	template< class Field, class Operand1, class Operand2>
	class SlicedPolynomialMatrixSubin
	{
	public:
		Operand1 &operator() (const Field &F, Operand1 &C, const Operand2 &A) const;

	};
} /* end of namespace LinBox */

#include "SlicedPolynomialMatrixAddSub.inl"

#endif

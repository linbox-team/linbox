#ifndef __LINBOX_matrix_SlicedPolynomialVector_SlicedPolynomialVectorAddSub_H
#define __LINBOX_matrix_SlicedPolynomialVector_SlicedPolynomialVectorAddSub_H

#include <vector>
#include "linbox/matrix/matrixdomain/vector-domain.h"
#include "linbox/matrix/slicedpolynomialvector/SlicedPolynomialVector.h"

namespace LinBox
{
	/* internal
	 * Adding two vectors
	 */
	template< class Field, class Operand1, class Operand2, class Operand3>
	class SlicedPolynomialVectorAdd
	{
	public:
		Operand1 &operator() (const Field &F, Operand1 &C, const Operand2 &A, const Operand3 &B) const;

	};

	/* internal
	 * Substracting two vectors
	 */
	template< class Field, class Operand1, class Operand2, class Operand3>
	class SlicedPolynomialVectorSub
	{
	public:
		Operand1 &operator() (const Field &F, Operand1 &C, const Operand2 &A, const Operand3 &B) const;

	};

	//! C += A
	template< class Field, class Operand1, class Operand2>
	class SlicedPolynomialVectorAddin
	{
	public:
		Operand1 &operator() (const Field &F, Operand1 &C, const Operand2 &A) const;

	};

	//! C -= A
	template< class Field, class Operand1, class Operand2>
	class SlicedPolynomialVectorSubin
	{
	public:
		Operand1 &operator() (const Field &F, Operand1 &C, const Operand2 &A) const;

	};
} /* end of namespace LinBox */

#include "SlicedPolynomialVectorAddSub.inl"

#endif

#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulToomCook_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulToomCook_H

#include "SlicedPolynomialMatrix.h"

namespace LinBox
{ 
	template< class Field, class Vector3, class Vector1, class Vector2>
	class SlicedPolynomialMatrixMulToomCook
	{
	private:
		typedef SlicedPolynomialMatrixtoBlasMatrix<Field, Vector1::Rep, Vector1::MatrixElement, Vector1::IntField, Vector1::Rep> conversionAtoAm;
		typedef SlicedPolynomialMatrixtoBlasMatrix<Field, Vector2::Rep, Vector2::MatrixElement, Vector2::IntField, Vector2::Rep> conversionBtoBm;
		typedef Vector3::IntField IntField;
		typedef BlasMatrix<IntField> Matrix;
		typedef BlasMatrixtoSlicedPolynomialMatrix<Vector3::IntField, Vector3::Rep, Field, Vector3::Rep, Vector3::MatrixElement> conversionCmtoC;
		typedef Vector3::polynomial polynomial;
		Matrix& EvaluationInterpolationMatrices (Matrix& TC, Matrix& iTC);
		template<class Field, class Vector3, class Vector1, class Vector2>
		Matrix& mul (IntField& F, Matrix& CMatBloc, const Matrix& AMatBloc, const Matrix& BMatBloc,
								   const size_t m, const size_t k, const size_t n, const size_t e, polynomial irreducible);
	public:
		Operand1 &operator() (const Field &GF, Operand1 &C, const Operand2 &A, const Operand3 &B) const;
	}; 
} /* end of namespace LinBox */

#include "SlicedPolynomialMatrixMulToomCook.inl"

#endif
#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixAddSub_INL
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixAddSub_INL

namespace LinBox
{
	template<class Field, class Vector1, class Vector2, class Vector3>
	Vector1& SlicedPolynomialMatrixAdd<Field, Vector1, Vector2, Vector3 >::operator()(const Field& F,
									   Vector1& C,
									   const Vector2& A,
									   const Vector3& B) const
	{
		//check dimensions
		BlasMatrix<Vector1::IntField> M(C.fieldF(), A.rowdim(), B.coldim());
		for (int m = 0; m < C.length(); m++)
		{
			BlasMatrixDomainAdd<Vector1::IntField, BlasMatrix<Vector1::IntField>, BlasMatrix<Vector2::IntField>, BlasMatrix<Vector3::IntField>>()(C.fieldF(),
				M, A.getMatrixCoefficient(m), B.getMatrixCoefficient(m));
			C.setMatrixCoefficient(m, M);
		}
		return C;
	}

	template<class Field, class Vector1, class Vector2, class Vector3>
	Vecto1r& SlicedPolynomialMatrixSub<Field, Vector1, Vector2, Vector3 >::operator()(const Field& F,
									   Vector1& C,
									   const Vector2& A,
									   const Vector3& B) const
	{
		//check dimensions
		BlasMatrix<Field> M(F, A.rowdim(), B.coldim());
		for (int m = 0; m < C.length(); m++)
		{
			BlasMatrixDomainSub<Vector1::IntField, BlasMatrix<Vector1::IntField>, BlasMatrix<Vector2::IntField>, BlasMatrix<Vector3::IntField>>()(C.fieldF(),
				M, A.getMatrixCoefficient(m), B.getMatrixCoefficient(m));
			C.setMatrixCoefficient(m, M);
		}
		return C;
	}

	template<class Field, class Vector1, class Vector3>
	Vector1& SlicedPolynomialMatrixAddin<Field, Vector1, Vector3 >::operator()(const Field& F,
									   Vector1& C,
									   const Vector3& B) const
	{
		//check dimensions
		BlasMatrix<Field> M(F, B.rowdim(), B.coldim());
		for (int m = 0; m < C.length(); m++)
		{
			M = C.getMatrixCoefficient(m);
			BlasMatrixDomainAddin<Vector1::IntField, BlasMatrix<Vector1::IntField>, BlasMatrix<Vector3::IntField>>()(C.fieldF(),
				M, A.getMatrixCoefficient(m), B.getMatrixCoefficient(m));
			C.setMatrixCoefficient(m, M);
		}
		return C;
	}

	template<class Field, class Vector1, class Vector3>
	Vector1& SlicedPolynomialMatrixSubin<Field, Vector1, Vector3 >::operator()(const Field& F,
									   Vector1& C,
									   const Vector3& B) const
	{
		//check dimensions
		BlasMatrix<Field> M(F, B.rowdim(), B.coldim());
		for (int m = 0; m < C.length(); m++)
		{
			M = C.getMatrixCoefficient(m);
			BlasMatrixDomainAddin<Vector1::IntField, BlasMatrix<Vector1::IntField>, BlasMatrix<Vector3::IntField>>()(C.fieldF(),
				M, A.getMatrixCoefficient(m), B.getMatrixCoefficient(m));
			C.setMatrixCoefficient(m, M);
		}
		return C;
	}
} // LinBox

#endif
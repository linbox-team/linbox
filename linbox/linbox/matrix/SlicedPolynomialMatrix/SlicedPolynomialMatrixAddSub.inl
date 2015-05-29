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
		BlasMatrix<Field> M(F, A.rowdim(), B.coldim());
		for (int m = 0; m < C.length(); m++)
		{
			BlasMatrixDomainAdd<Field, Modular<int>, Modular<int>, Modular<int>>()(C.fieldF(),
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
			BlasMatrixDomainSub<Field, Modular<int>, Modular<int>, Modular<int>>()(C.fieldF(),
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
			BlasMatrixDomainAddin<Field, Modular<int>, Modular<int>>()(C.fieldF(),
				M, B.getMatrixCoefficient(m));
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
			BlasMatrixDomainSubin<Field, Modular<int>, Modular<int>>()(C.fieldF(),
				M, B.getMatrixCoefficient(m));
			C.setMatrixCoefficient(m, M);
		}
		return C;
	}
} // LinBox

#endif
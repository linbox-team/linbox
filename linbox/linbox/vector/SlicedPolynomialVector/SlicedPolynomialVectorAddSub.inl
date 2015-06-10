#ifndef __LINBOX_matrix_SlicedPolynomialVector_SlicedPolynomialVectorAddSub_INL
#define __LINBOX_matrix_SlicedPolynomialVector_SlicedPolynomialVectorAddSub_INL

namespace LinBox
{
	template<class Field, class Vector1, class Vector2, class Vector3>
	Vector1& SlicedPolynomialVectorAdd<Field, Vector1, Vector2, Vector3 >::operator()(const Field& F,
									   Vector1& C,
									   const Vector2& A,
									   const Vector3& B) const
	{
		//check dimensions
		BlasVector<Vector1::IntField> M(C.fieldF(), C.rowdim());
		for (int m = 0; m < C.length(); m++)
		{
			VectorDomain<Vector1::IntField>.add(C.fieldF(),
				M, A.getVectorCoefficient(m), B.getVectorCoefficient(m));
			C.setVectorCoefficient(m, M);
		}
		return C;
	}

	template<class Field, class Vector1, class Vector2, class Vector3>
	Vecto1r& SlicedPolynomialVectorSub<Field, Vector1, Vector2, Vector3 >::operator()(const Field& F,
									   Vector1& C,
									   const Vector2& A,
									   const Vector3& B) const
	{
		//check dimensions
		BlasVector<Vector1::IntField> M(C.fieldF(), C.rowdim());
		for (int m = 0; m < C.length(); m++)
		{
			VectorDomain<Vector1::IntField>.sub(C.fieldF(),
				M, A.getVectorCoefficient(m), B.getVectorCoefficient(m));
			C.setVectorCoefficient(m, M);
		}
		return C;
	}

	template<class Field, class Vector1, class Vector3>
	Vector1& SlicedPolynomialVectorAddin<Field, Vector1, Vector3 >::operator()(const Field& F,
									   Vector1& C,
									   const Vector3& B) const
	{
		//check dimensions
		BlasVector<Vector1::IntField> M(C.fieldF(), C.rowdim());
		for (int m = 0; m < C.length(); m++)
		{
			M = C.getVectorCoefficient(m);
			VectorDomain<Vector1::IntField>.addin(C.fieldF(),
				M, A.getVectorCoefficient(m), B.getVectorCoefficient(m));
			C.setVectorCoefficient(m, M);
		}
		return C;
	}

	template<class Field, class Vector1, class Vector3>
	Vector1& SlicedPolynomialVectorSubin<Field, Vector1, Vector3 >::operator()(const Field& F,
									   Vector1& C,
									   const Vector3& B) const
	{
		//check dimensions
		BlasVector<Vector1::IntField> M(C.fieldF(), C.rowdim());
		for (int m = 0; m < C.length(); m++)
		{
			M = C.getVectorCoefficient(m);
			VectorDomain<Vector1::IntField>.add(C.fieldF(),
				M, A.getVectorCoefficient(m), B.getVectorCoefficient(m));
			C.setVectorCoefficient(m, M);
		}
		return C;
	}
} // LinBox

#endif
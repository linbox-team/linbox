#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulKaratsuba_INL
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulKaratsuba_INL

#include <givaro/givpoly1dense.h>
#include <givaro/givpoly1denseops.inl>
#include "linbox/matrix/matrixdomain/blas-matrix-domain.h"
#include "linbox/vector/vector-domain.h"
#include "fflas-ffpack/fflas/fflas.h"

namespace LinBox
{
	template<class Field, class Vector3, class Vector1, class Vector2>
	SlicedPolynomialMatrixMulKaratsuba<Field, Vector3, Vector1, Vector2 >::vec&
		SlicedPolynomialMatrixMulKaratsuba<Field, Vector3, Vector1, Vector2 >::karatsuba(IntField& F, VectorDomain<Vector3::IntField>& VD, 
		                                                                     vecV& C, vecM& A, vecV& B, size_t rd, size_t cd)
	{
		if (A.size() == 1)
		{
			for (int i = 0; i <B.size(); i++)
			{
				FFLAS::fgemv(F, FflasNoTrans, rd, cd, F.one, A[0].getPointer(), cd, B[i].getPointer(), 1, F.zero, C[i].getPointer(), 1);
			}
			return C;
		}
		if (B.size() == 1)
		{
			for (int i = 0; i < A.size(); i++)
			{
				FFLAS::fgemv(F, FflasNoTrans, rd, cd, F.one, A[i].getPointer(), cd, B[0].getPointer(), 1, F.zero, C[i].getPointer(), 1);
			}
			return C;
		}
		int m = (A.size() < B.size()) ? (B.size() / 2) : (A.size() / 2);
		if ((m < A.size()) && (m < B.size()))
		{
			vecM A1(A.begin(), A.begin() + m);
			vecM A2(A.begin() + m, A.end());
			vecV B1(B.begin(), B.begin() + m);
			vecV B2(B.begin() + m, B.end());
			vecM A3;
			int minlength_a = (A1.size() < A2.size()) ? A1.size() : A2.size();
			int minlength_a = (A1.size() < A2.size()) ? A2.size() : A1.size();
			for (int i = 0; i < minlength_a; i++)
			{
				Matrix AA(F, rd, cd);
				A3.push_back(BlasMatrixDomainAdd<IntField, Matrix, Matrix, Matrix>()(F, AA, A1[i], A2[i]));
			}
			if (maxlength_a == A1.size())
			{
				for (int i = minlength_a; i < maxlength_a; i++)
				{
					A3.push_back(A1[i]);
				}
			}
			else
			{
				for (int i = minlength_a; i < maxlength_a; i++)
				{
					A3.push_back(A2[i]);
				}
			}
			vecV B3;
			int minlength_b = (B1.size() < B2.size()) ? B1.size() : B2.size();
			int maxlength_b = (B1.size() < B2.size()) ? B2.size() : B1.size();
			for (int i = 0; i < minlength_b; i++)
			{
				Vector BB(F, cd);
				B3.push_back(VD.add<Vector, Vector, Vector>(BB, B1[i], B2[i]));
			}
			if (maxlength_b == B1.size())
			{
				for (int i = minlength_b; i < maxlength_b; i++)
				{
					B3.push_back(B1[i]);
				}
			}
			else
			{
				for (int i = minlength_b; i < maxlength_b; i++)
				{
					B3.push_back(B1[i]);
				}
			}
			vecV C1;
			vecV C2;
			vecV C3;
			Matrix CC(F, rd);
			int xx;
			xx = A1.size() + B1.size() - 1;
			for (int i = 0; i < xx; i++)
			{
				C1.push_back(CC);
			}
			xx = A2.size() + B2.size() - 1;
			for (int i = 0; i < xx; i++)
			{
				C2.push_back(CC);
			}
			xx = A3.size() + B3.size() - 1;
			for (int i = 0; i < xx; i++)
			{
				C3.push_back(CC);
			}
			karatsuba(F, VD, C1, A1, B1);
			karatsuba(F, VD, C2, A2, B2);
			karatsuba(F, VD, C3, A3, B3);
			for (int i = 0; i < C1.size(); i++)
			{
				VD.addin<Vector, Vector>(C[i], C1[i]);
			}
			int mm = 2 * m;
			for (int i = 0; i < C2.size(); i++)
			{
				VD.addin<Vector, Vector>(C[mm + i], C2[i]);
			}
			for (int i = 0; i < C3.size(); i++)
			{
				VD.addin<Vector, Vector>(C[m + i], C3[i]);
			}
			for (int i = 0; i < C1.size(); i++)
			{
				VD.subin<Vector, Vector>(C[m + i], C1[i]);
			}
			for (int i = 0; i < C2.size(); i++)
			{
				VD.subin<Vector, Vector>(C[m + i], C2[i]);
			}
			return C;
		}
		if (A.size() <= m)
		{
			vecV B1(B.begin(), B.begin() + m);
			vecV B2(B.begin() + m, B.end());
			vecV C1;
			vecV C2;
			Vector CC(F, rd);
			int xx;
			xx = A.size() + B1.size() - 1;
			for (int i = 0; i < xx; i++)
			{
				C1.push_back(CC);
			}
			xx = A.size() + B2.size() - 1;
			for (int i = 0; i < xx; i++)
			{
				C2.push_back(CC);
			}
			karatsuba(F, VD, C1, A, B1);
			karatsuba(F, VD, C2, A, B2);
			for (int i = 0; i < C1.size(); i++)
			{
				VD.addin<Vector, Vector>(C[i], C1[i]);
			}
			for (int i = 0; i < C2.size(); i++)
			{
				VD.addin<Vector, Vector>(C[m + i], C2[i]);
			}
			return C;
		}
		if (B.size() <= m)
		{
			vecM A1(A.begin(), A.begin() + m);
			vecM A2(A.begin() + m, A.end());
			vecV C1;
			vecV C2;
			Matrix CC(F, rd);
			int xx;
			xx = A1.size() + B.size() - 1;
			for (int i = 0; i < xx; i++)
			{
				C1.push_back(CC);
			}
			xx = A2.size() + B.size() - 1;
			for (int i = 0; i < xx; i++)
			{
				C2.push_back(CC);
			}
			karatsuba(F, VD, C1, A1, B);
			karatsuba(F, VD, C2, A2, B);
			for (int i = 0; i < C1.size(); i++)
			{
				VD.addin<Vector, Vector>(C[i], C1[i]);
			}
			for (int i = 0; i < C2.size(); i++)
			{
				VD.addin<Vector, Vector>(C[m + i], C2[i]);
			}
			return C;
		}
		return C;
	}

	template<class Field, class Vector3, class Vector1, class Vector2>
	Vector1& SlicedPolynomialMatrixMulKaratsuba<Field, Vector3, Vector1, Vector2 >::operator()(const Field& GF,
									   Vector3& C,
									   const Vector1& A,
									   const Vector2& B) const
	{
		//check dimensions
		int xx;
		vecV A1;
		vecM B1;
		xx = A.length();
		for (int m = 0; m < xx; m++)
		{
			A1.push_back(A.getMatrixCoefficient(m));
		}
		xx = B.length();
		for (int m = 0; m < xx; m++)
		{
			B1.push_back(B.getVectorCoefficient(m));
		}
		vecV C1;
		xx = A1.size() + B1.size() - 1;
		Vector CC(C.fieldF(), A.rowdim());
		for (int i = 0; i < xx; i++)
		{
			C1.push_back(CC);
		}
		VectorDomain<IntField> VD(F);
		karatsuba(C.fieldF(), VD, C1, A1, B1);
		int mk = C1.size();
		int mi = C1[0].rowdim();
		int n = C.size();
		Poly1Dom<IntField, Dense> PD(F);
		for (int i = 0; i < mi; i++)
		{
			
				polynomial entry;
				for (int k = 0; k < mk; k++)
				{
					entry.push_back(C1[k].getEntry(i));
				}
				PD.modin(entry, C.irreducible);
				for (int k = 0; k < n; k++)
				{
					C[k].setEntry(i, w1[k]);
				}
			
		}
		return C;
	}
} // LinBox

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

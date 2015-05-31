#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulKaratsuba_INL
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulKaratsuba_INL

#include <givaro/givpoly1dense.h>
#include <givaro/givpoly1denseops.inl>

namespace LinBox
{
	template<class Field, class Vector3, class Vector1, class Vector2>
	SlicedPolynomialMatrixMulKaratsuba<Field, Vector3, Vector1, Vector2 >::vec& 
		SlicedPolynomialMatrixMulKaratsuba<Field, Vector3, Vector1, Vector2 >::modulo(vec& C, int n, polynomial irreducible)
	{
		int m = C.size();
		int mi = C[0].rowdim();
		int mj = C[0].coldim();
		vec result(n, mi, mj);
		for (int i = 0; i < mi; i++)
		{
			for (int j = 0; j < mj; j++)
			{
				polynomial entry;
				for (int k = 0; k < mk; k++)
				{
					entry.push_back(C[k].getEntry(i, j));
				}
				polynomial w1;
				Poly1Dom<IntField,Dense>::div(w1, entry, irreducible);
				polynomial w2;
				Poly1Dom<IntField,Dense>::mul(w2, w1, irreducible);
				Poly1Dom<IntField,Dense>::sub(w1, entry, w2);
				for (int k = 0; k < n; k++)
				{
					result[k].setEntry(i, j, w1[k]);
				}
			}
		}
		C = result;
		return C;
	}

	template<class Field, class Vector3, class Vector1, class Vector2>
	SlicedPolynomialMatrixMulKaratsuba<Field, Vector3, Vector1, Vector2 >::vec&
		SlicedPolynomialMatrixMulKaratsuba<Field, Vector3, Vector1, Vector2 >::karatsuba(IntField& F, vec& C, vec& A, vec& B)
	{
		if (A.size() == 1)
		{
			for (int i = 0; i <B.size(); i++)
			{
				BlasMatrixDomainMul<IntField, Matrix, Matrix, Matrix>()(F, C[i], A[0], B[i]);
			}
			return C;
		}
		if (B.size() == 1)
		{
			for (int i = 0; i < A.size(); i++)
			{
				BlasMatrixDomainMul<IntField, Matrix, Matrix, Matrix>()(F, C[i], A[i], B[0]);
			}
			return C;
		}
		int m = (A.size() < B.size()) ? (B.size() / 2) : (A.size() / 2);
		if ((m < A.size()) && (m < B.size()))
		{
			vec A1(A.begin(), A.begin() + m);
			vec A2(A.begin() + m, A.end());
			vec B1(B.begin(), B.begin() + m);
			vec B2(B.begin() + m, B.end());
			vec A3;
			int minlength_a = (A1.size() < A2.size()) ? A1.size() : A2.size();
			int minlength_a = (A1.size() < A2.size()) ? A2.size() : A1.size();
			for (int i = 0; i < minlength_a; i++)
			{
				Matrix AA(F, A.rowdim(), A.coldim());
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
			vec B3;
			int minlength_b = (B1.size() < B2.size()) ? B1.size() : B2.size();
			int maxlength_b = (B1.size() < B2.size()) ? B2.size() : B1.size();
			for (int i = 0; i < minlength_b; i++)
			{
				Matrix BB(F, B.rowdim(), B.coldim());
				B3.push_back(BlasMatrixDomainAdd<IntField, Matrix, Matrix, Matrix>()(F, BB, B1[i], B2[i]));
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
			vec C1;
			vec C2;
			vec C3;
			Matrix CC(F, A.rowdim(), B.coldim());
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
			karatsuba(F, C1, A1, B1);
			karatsuba(F, C2, A2, B2);
			karatsuba(F, C3, A3, B3);
			for (int i = 0; i < C1.size(); i++)
			{
				BlasMatrixDomainAddin<IntField, Matrix, Matrix>()(F, C[i], C1[i]);
			}
			int mm = 2 * m;
			for (int i = 0; i < C2.size(); i++)
			{
				BlasMatrixDomainAddin<IntField, Matrix, Matrix>()(F, C[mm + i], C2[i]);
			}
			for (int i = 0; i < C3.size(); i++)
			{
				BlasMatrixDomainAddin()<IntField, Matrix, Matrix>(F, C[m + i], C3[i]);
			}
			for (int i = 0; i < C1.size(); i++)
			{
				BlasMatrixDomainSubin()<IntField, Matrix, Matrix>(F, C[m + i], C1[i]);
			}
			for (int i = 0; i < C2.size(); i++)
			{
				BlasMatrixDomainSubin()<IntField, Matrix, Matrix>(F, C[m + i], C2[i]);
			}
			return C;
		}
		if (A.size() <= m)
		{
			vec B1(B.begin(), B.begin() + m);
			vec B2(B.begin() + m, B.end());
			vec C1;
			vec C2;
			Matrix CC(F, A.rowdim(), B.coldim());
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
			karatsuba(F, C1, A, B1);
			karatsuba(F, C2, A, B2);
			for (int i = 0; i < C1.size(); i++)
			{
				BlasMatrixDomainAddin<IntField, Matrix, Matrix>()(F, C[i], C1[i]);
			}
			for (int i = 0; i < C2.size(); i++)
			{
				BlasMatrixDomainAddin<IntField, Matrix, Matrix>()(F, C[m + i], C2[i]);
			}
			return C;
		}
		if (B.size() <= m)
		{
			vec A1(A.begin(), A.begin() + m);
			vec A2(A.begin() + m, A.end());
			vec C1;
			vec C2;
			Matrix CC(F, A.rowdim(), B.coldim());
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
			karatsuba(F, C1, A1, B);
			karatsuba(F, C2, A2, B);
			for (int i = 0; i < C1.size(); i++)
			{
				BlasMatrixDomainAddin<IntField, Matrix, Matrix>()(F, C[i], C1[i]);
			}
			for (int i = 0; i < C2.size(); i++)
			{
				BlasMatrixDomainAddin<IntField, Matrix, Matrix>()(F, C[m + i], C2[i]);
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
		vec A1;
		vec B1;
		xx = A.length();
		for (int m = 0; m < xx; m++)
		{
			A1.push_back(A.getMatrixCoefficient(m));
		}
		xx = B.length();
		for (int m = 0; m < xx; m++)
		{
			B1.push_back(B.getMatrixCoefficient(m));
		}
		vec C1;
		xx = A1.size() + B1.size() - 1;
		Matrix CC(C.fieldF(), A.rowdim(), B.coldim());
		for (int i = 0; i < xx; i++)
		{
			C1.push_back(CC);
		}
		karatsuba(C.fieldF(), C1, A1, B1);
		modulo(C1, C.length(), C.irreducible);
		for (int m = 0; m < C.length(); m++)
		{
			C.setMatrixCoefficient(m, C1[m]);
		}
		return C;
	}
} // LinBox

#endif
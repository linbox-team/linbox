#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulToomCook_INL
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulToomCook_INL

#include <fflas-ffpack/fflas/fflas.h>
#include <fflas-ffpack/ffpack/ffpack.h>
#include "linbox/algorithms/matrix-hom.h"
#include <givaro/givpoly1denseops.inl>
#include "linbox/matrix/SlicedPolynomialMatrix/conversion.h"

namespace LinBox
{
	template<class Field, class Vector3, class Vector1, class Vector2>
	BlasMatrix<Vector3::IntField>&
	SlicedPolynomialMatrixMulToomCook<Field, Vector3, Vector1, Vector2 >::EvaluationInterpolationMatrices
				(BlasMatrix<Vector3::IntField>& TC, BlasMatrix<Vector3::IntField>& iTC)
	{
		size_t E = TC.rowdim();
		for (size_t i = 0 ; i < E ; ++i)
		{
			for (size_t j = 0 ; j < E ; ++j)
			{
				TC.field().init(TC.refEntry(i,j), pow((Integer)i,j));
			}
		}
		int null;
		FFPACK::Invert(TC.field(),E,TC.getPointer(),E,iTC.getWritePointer(),E,null);
		return TC;
	}

	template<class Field, class Vector3, class Vector1, class Vector2>
		BlasMatrix<Vector3::IntField>&
		SlicedPolynomialMatrixMulToomCook<Field, Vector3, Vector1, Vector2 >::mul
								   (Vector3::IntField& F,
							       BlasMatrix<Vector3::IntField>& CMatBloc,
								   const BlasMatrix<Vector3::IntField>& AMatBloc,
								   const BlasMatrix<Vector3::IntField>& BMatBloc,
								   const size_t m,
								   const size_t k,
								   const size_t n, const size_t e,
								   polynomial irreducible)
	{
		#if (__LINBOX_FFLAS_FFPACK_VERSION < 10501)
				#warning "Invert is buggy in your fflas-ffpack version. please consider upgrading to >=1.5.1."
		#endif
		size_t E = 2*e - 1 ;

		Matrix TC    (F, E, E);
		Matrix iTC   (F, E, E);
		Matrix iEval (F, E, E);
		EvaluationInterpolationMatrices(TC, iTC);

		// each row is a result matrix
		Matrix TMatBloc( F, E, m*n);
		Matrix AEval( F , E, m*k);
		Matrix BEval( F , E, k*n);

		FFLAS::fgemm(F,
					 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
					 E, m*k, e,
					 F.one,
					 TC.getPointer(),E,
					 AMatBloc.getPointer(), m*k,
					 F.zero,
					 AEval.getWritePointer(), m*k);

		FFLAS::fgemm(F,
					 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
					 E ,n*k, e,
					 F.one,
					 TC.getPointer(),E,
					 BMatBloc.getPointer(), n*k,
					 F.zero,
					 BEval.getWritePointer(), n*k);

		for (size_t i = 0 ; i < E ; ++i)
		{
			FFLAS::fgemm(F,
						 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
						 m, n, k,
						 F.one,
						 AEval.getPointer()+i*m*k, k,
						 BEval.getPointer()+i*n*k, n,
						 F.zero,
						 TMatBloc.getWritePointer()+i*m*n, n);
		}

		FFLAS::fgemm(F,
					 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				     E, m * n, E,
					 F.one,
					 iTC.getPointer(),E, //lda
					 TMatBloc.getWritePointer()+i*m*n, m*n, //ldb
					 F.zero,
					 CMatBloc.getWritePointer()+i*m*n, m*n);

		return CMatBloc;
	}

	// all matrix classes should be SlicedPolynomialMatrices
	template<class Field, class Vector3, class Vector1, class Vector2>
	Vector3& SlicedPolynomialMatrixMulToomCook<Field, Vector3, Vector1, Vector2 >::operator()
									   (const Field& GF,
									   Vector3& C,
									   const Vector1& A,
									   const Vector2& B)
	{
		size_t e = C.length();
		size_t m = C.rowdim();
		size_t k = B.rowdim();
		size_t n = C.coldim();

		IntField F = A.filedF();

		if (e == 1) {
			BlasMatrix<Vector1::IntField, Vector1::Rep> Am(A.fieldF(), A.rowdim(), A.coldim());
			conversionAtoAm()(Am, A);
			BlasMatrix<Vector2::IntField, Vector2::Rep> Bm(B.fieldF(), B.rowdim(), B.coldim());
			conversionBtoBm()(Bm, B);
			BlasMatrix<Vector3::IntField, Vector3::Rep> Cm(C.fieldF(), C.rowdim(), C.coldim());
			Matrix Af(Am,F);
			Matrix Bf(Bm,F);
			Matrix Cf(Cm,F);
			FFLAS::fgemm((typename IntField::Father_t)F,
						 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
						 // m, k , n,
						 m,n,k,
						 F.one,
						 Af.getPointer(), Am.getStride(), //lda
						 Bf.getPointer(), Bm.getStride(), //ldb
						 F.zero,
						 Cf.getWritePointer(), Cm.getStride());
			MatrixHom::map(Cm,Cf);
			conversionCmtoC()(C, Cm); //change when convertation is done
			return C;
		}

		Matrix Cbloc(F,e,m*n);
		Matrix Abloc(F,e,m*k);
		Matrix Bbloc(F,e,k*n);


		for (size_t l = 0 ; l < e ; ++l)
		{
			for (size_t i = 0 ; i < m ; ++i)
			{
				for (size_t j = 0 ; j < k ; ++j)
				{
					Abloc.setEntry(l, i*k+j, A.getEntry(l, i, j));
				}
			}
		}

		for (size_t l = 0 ; l < e ; ++l)
		{
			for (size_t i = 0 ; i < k ; ++i)
			{
				for (size_t j = 0 ; j < n ; ++j)
				{
					Bbloc.setEntry(l, i*n+j, B.getEntry(l, i, j));
				}
			}
		}

		Protected::mul(C.fieldF(), Cbloc, Abloc, Bbloc, m, k, n, e, C.irreducible);
		
		size_t E = 2 * E - 1;
		polynomial x(E);
		for (size_t i = 0 ; i < m ; ++i)
		{
			for (size_t j = 0 ; j < n ; ++j)
			{
				for (size_t l = 0 ; l < E ; ++l)
				{
					x[l] = Cbloc.getEntry(l, i*n+j);
				}
				Poly1Dom<IntField>.modin(x, C.irreducible);
				for (size_t l = 0 ; l < e ; ++l)
				{
					C.setEntry(l, i, j, x[l]);
				}
			}
		}

		return C ;
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

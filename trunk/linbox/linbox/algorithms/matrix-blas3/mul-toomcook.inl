/*  Copyright (C) 2012 the members of the LinBox group
 * Written by B. Boyer < bboyer@imag.fr >
 *
 * This file is part of the LinBox library.
 *
 * ========LICENCE========
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * LinBox is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
 */

/*! @internal
 * @file matrix-blas3/mul-toomcook.inl
 * @ingroup algorithm
 * @ingroup blas
 * @brief Implementation of Toom-Cook.
 */

#ifndef __LINBOX_matrix_blas3_mul_toomcook_INL
#define __LINBOX_matrix_blas3_mul_toomcook_INL



#include <fflas-ffpack/fflas/fflas.h>
#include <fflas-ffpack/ffpack/ffpack.h>

namespace LinBox {
	namespace BLAS3 {

		template<class Zpz>
		BlasMatrix<Zpz>& ToomCook(BlasMatrix<Zpz>& TC, BlasMatrix<Zpz>& iTC)
		{
			linbox_check(TC.rowdim() == TC.coldim());
			size_t l = TC.rowdim();

			for (size_t i = 0 ; i < l ; ++i)
				for (size_t j = 0 ; j < l ; ++j) {
					TC.field().init(TC.refEntry(i,j), pow((Integer)i,j));
				}
			BlasMatrixDomain<Zpz> BMD(TC.field()) ;
			// TC.write(std::cout << "TC ") << std::endl;
			// BMD.invert(iTC, TC);
			int null;
			FFPACK::Invert(TC.field(),l,TC.getPointer(),l,iTC.getWritePointer(),l,null);
			// iTC.write(std::cout << "TC^(-1) ") << std::endl;
			// TC.write(std::cout << "TC ") << std::endl;
			return TC;
		}

		namespace Protected {


			template<class Zpz, class GFpe>
			BlasMatrix<Zpz >& mul (BlasMatrix<Zpz>      & CMatBloc,
								   const BlasMatrix<Zpz>& AMatBloc,
								   const BlasMatrix<Zpz>& BMatBloc,
								   const size_t m,
								   const size_t k,
								   const size_t n,
								   const mulMethod::ToomCook<GFpe> & T)
			{
				const Zpz & F  = CMatBloc.field();
				const GFpe& GF = T._myF ;
				// linbox_check(T._myF.characacteristic() == F.characteristic());
				// TODO si e = 1 on matmul !
				size_t e = (size_t) GF.exponent() ; // extension degree
				size_t l = 2*e - 1 ; // sure ?

				BlasMatrix<Zpz> TC    (F,l,l);
				BlasMatrix<Zpz> iTC   (F,l,l);
				BlasMatrix<Zpz> iEval (F,l,l);
				// FWD = Matrix(K, l, l, [K(i**j) for i in range(l) for j in range(l)])
				ToomCook(TC,iTC);
				// each row is a result matrix
				BlasMatrix<Zpz> TMatBloc( F, l, m*n);


				// AY = [sum(FWD[i,j]*A[j] for j in range(len(A))) for i in range(l)]
				// BY = [sum(FWD[i,j]*B[j] for j in range(len(B))) for i in range(l)]
				if (!T.memory_unlimited)	{			/* space efficient */
					BlasMatrix<Zpz> AEval( F , m, k);
					BlasMatrix<Zpz> BEval( F , k, n);


					for (size_t i = 0 ; i < l ; ++i) {
						FFLAS::fgemv(F, FFLAS::FflasTrans,
									 e, m*k,
									 F.one,
									 AMatBloc.getPointer(), m*k,
									 TC.getPointer()+ i*l, 1,
									 F.zero,
									 AEval.getWritePointer(), 1);

						FFLAS::fgemv(F, FFLAS::FflasTrans,
									 e, k*n,
									 F.one,
									 BMatBloc.getPointer(), k*n,
									 TC.getPointer()+ i*l, 1,
									 F.zero,
									 BEval.getWritePointer(), 1);

						FFLAS::fgemm(F,
									 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
									 // m, k , n,
									 m,n,k,
									 F.one,
									 AEval.getPointer(), AEval.coldim(), //lda
									 BEval.getPointer(), BEval.coldim(), //ldb
									 F.zero,
									 TMatBloc.getWritePointer()+i*m*n, n);
					}
				}
				else { /* time efficient (matmul) */
					BlasMatrix<Zpz> AEval( F , l, m*k);
					BlasMatrix<Zpz> BEval( F , l, k*n);


					FFLAS::fgemm(F,
								 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
								 // m, k , n,
								 l,m*k,e,
								 F.one,
								 // AEval.getPointer(), AEval.coldim(), //lda
								 TC.getPointer(),l,
								 AMatBloc.getPointer(), m*k,
								 F.zero,
								 AEval.getWritePointer(), m*k);
					// TMatBloc.getWritePointer()+i*m*n, n);

					FFLAS::fgemm(F,
								 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
								 // m, k , n,
								 l,n*k,e,
								 F.one,
								 // AEval.getPointer(), AEval.coldim(), //lda
								 TC.getPointer(),l,
								 BMatBloc.getPointer(), n*k,
								 F.zero,
								 BEval.getWritePointer(), n*k);

					for (size_t i = 0 ; i < l ; ++i) {

						FFLAS::fgemm(F,
									 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
									 // m, k , n,
									 m,n,k,
									 F.one,
									 AEval.getPointer()+i*m*k, k, //lda
									 BEval.getPointer()+i*n*k, n, //ldb
									 F.zero,
									 TMatBloc.getWritePointer()+i*m*n, n);
					}
				}

				std::vector<double>  minpoly_vector  =  GF.irreducible() ;
				for (size_t i = 0 ; i < minpoly_vector.size(); ++i)
					F.negin(minpoly_vector[i]);

				// BlasSubmatrix<Zpz> CompMat(CMatBloc,0,0,l,l);
				BlasMatrix<Zpz> CompMat(F,l,l);
				for (size_t i = 0 ; i < e ; ++i) { // degree == l ?
					CompMat.setEntry(i,i,F.one);
				}
				for (size_t j = 0 ; j < e-1 ; ++j) {
					CompMat.setEntry(j,e,minpoly_vector[j]);
				}
				typename Zpz::Element coeff ;
				typename Zpz::Element tmp_coeff, tmp_coeff2 ;
				F.init(tmp_coeff);
				for (size_t i = 1 ; i < e-1 ; ++i){
					for (size_t j = 1 ; j < e+1 ; ++j) {
						CompMat.setEntry(j, i + e, CompMat.getEntry((j-1),i - 1 + e) ) ;
					}
					if (!F.isZero(CompMat.getEntry(e,i+e))){
						F.init(coeff, CompMat.getEntry( e, i + e) );
						for(size_t j = 0 ; j < e-1 ; ++j){
							F.mul(tmp_coeff, coeff,minpoly_vector[j]);
							F.init(tmp_coeff2,CompMat.getEntry(j,i+e));
							F.addin(tmp_coeff2, tmp_coeff);
							CompMat.setEntry(j,i+e, tmp_coeff2 );
						}
					}
				}

				// BCK = ~FWD
				// XXX some stuff here
				FFLAS::fgemm(F,
							 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
							 l, l, l, F.one,
							 CompMat.getPointer(), l,
							 iTC.getPointer(), l,
							 F.zero,
							 iEval.getWritePointer(), l);


				// Y = [sum(BCK[i,j]*Y[j] for j in range(l)) for i in range(l)]
				FFLAS::fgemm(F,
							 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
							 e,m*n,l,
							 F.one,
							 iEval.getPointer(), l,
							 TMatBloc.getPointer(), m*n,
							 F.zero,
							 CMatBloc.getWritePointer(), m*n);
				return CMatBloc;

			}
		}

		template<class Zpz, class GF>
		std::vector<BlasMatrix<Zpz> >& mul (std::vector<BlasMatrix<Zpz> >& C,
											const std::vector<BlasMatrix<Zpz> >& A,
											const std::vector<BlasMatrix<Zpz> >& B,
											const mulMethod::ToomCook<GF> & T)
		{
			size_t m = C[0].rowdim();
			size_t k = B[0].rowdim();
			size_t n = C[0].coldim();
			const Zpz & F = C[0].field();
			BlasMatrix<Zpz> Cbloc(F,C.size(),m*n);
			BlasMatrix<Zpz> Abloc(F,C.size(),m*k);
			BlasMatrix<Zpz> Bbloc(F,C.size(),k*n);
			// convert
			for (size_t l = 0 ; l < C.size() ; ++l){
				for (size_t i = 0 ; i < m ; ++i) {
					for (size_t j = 0 ; j < k ; ++j) {
						Abloc.setEntry(l,i*k+j,A[l].getEntry(i,j));
					}
				}
			}

			for (size_t l = 0 ; l < C.size() ; ++l){
				for (size_t i = 0 ; i < k ; ++i) {
					for (size_t j = 0 ; j < n ; ++j) {
						Bbloc.setEntry(l,i*n+j,B[l].getEntry(i,j));
					}
				}
			}


			Protected::mul(Cbloc,Abloc,Bbloc,m,k,n,T);

			for (size_t l = 0 ; l < C.size() ; ++l){
				for (size_t i = 0 ; i < m ; ++i) {
					for (size_t j = 0 ; j < n ; ++j) {
						C[l].setEntry(i,j,Cbloc.getEntry(l,i*n+j));
					}
				}
			}
			// convert back
		}

		template<class Zpz>
		BlasMatrix<GivaroExtension<Zpz> >&
		mul (BlasMatrix<GivaroExtension<Zpz> >& C,
			 const BlasMatrix<GivaroExtension<Zpz> >& A,
			 const BlasMatrix<GivaroExtension<Zpz> >& B,
			 const mulMethod::ToomCook<GivaroExtension<Zpz> > & T)
		{
			size_t m = C.rowdim();
			size_t k = B.rowdim();
			size_t n = C.coldim();
			Zpz F ( A.field().characteristic() ); // BaseField ?
			size_t e = (size_t) A.field().exponent();
			BlasMatrix<Zpz> Cbloc(F,e,m*n);
			BlasMatrix<Zpz> Abloc(F,e,m*k);
			BlasMatrix<Zpz> Bbloc(F,e,k*n);


			for (size_t l = 0 ; l < e ; ++l){
				for (size_t i = 0 ; i < m ; ++i) {
					for (size_t j = 0 ; j < k ; ++j) {
						if (l< A.getEntry(i,j).size())
							Abloc.setEntry(l,i*k+j,
										   A.getEntry(i,j)[l]);
					}
				}
			}
			// for (size_t l = 0 ; l < e ; ++l){
				// std::cout << Abloc.getEntry(l,0) << ';';
			// }
			// std::cout << std::endl;

			for (size_t l = 0 ; l < e ; ++l){
				for (size_t i = 0 ; i < k ; ++i) {
					for (size_t j = 0 ; j < n ; ++j) {
						if (l< B.getEntry(i,j).size())
							Bbloc.setEntry(l,i*n+j,B.getEntry(i,j)[l]);
					}
				}
			}

			// C.field().init(C.refEntry(0,0));

			Protected::mul(Cbloc,Abloc,Bbloc,m,k,n,T);
			// convert back

			typedef typename GivaroExtension<Zpz>::Element Element ;
			for (size_t i = 0 ; i < m ; ++i) {
				for (size_t j = 0 ; j < n ; ++j) {
					Element x(e) ;
					for (size_t l = 0 ; l < e ; ++l){
						x[l] = Cbloc.getEntry(l,i*n+j);
					}
					C.field().polynomial_domain().modin(x,C.field().irreducible());
					// A.field().convert((Element&)C.refEntry(i,j),x);
					// A.field().init(x);
					C.setEntry(i,j,x);
				}
			}
			// for (size_t l = 0 ; l < e ; ++l){
				// std::cout << Cbloc.getEntry(l,0) << ';';
			// }
			// std::cout << std::endl;

		}

	}
}

#endif // __LINBOX_matrix_blas3_mul_toomcook_INL

//Local Variables:
//mode: C++
//tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


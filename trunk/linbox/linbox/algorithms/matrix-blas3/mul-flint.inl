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

#ifndef __LINBOX_matrix_blas3_mul_flint_INL
#define __LINBOX_matrix_blas3_mul_flint_INL


namespace LinBox {
	namespace BLAS3 {
		namespace Protected {
			template<class _ZZ>
			void getFlintMatrix(BlasMatrix<_ZZ> &C , const FLINT::fmpz_mat_t Cf)
			{
				Integer toto(0) ;
				for (size_t i = 0 ; i< C.rowdim() ; ++i)
					for (size_t j = 0 ; j< C.coldim() ; ++j) {
						FLINT::fmpz_get_mpz(toto.get_mpz(), fmpz_mat_entry(Cf,i,j));
						C.setEntry(i,j, toto );
					}
			}

			template<class _ZZ>
			void setFlintMatrix(FLINT::fmpz_mat_t Cf, const BlasMatrix<_ZZ> &C )
			{
				Integer toto(0) ;
				for (size_t i = 0 ; i< C.rowdim() ; ++i)
					for (size_t j = 0 ; j< C.coldim() ; ++j) {
						FLINT::fmpz_set_mpz(fmpz_mat_entry(Cf,i,j), C.getEntry(i,j).get_mpz_const());
					}
			}

#if 0
			template<class _Mat>
			void getFlintMatrix(BlasSubmatrix<_Mat> &C , const FLINT::fmpz_mat_t Cf)
			{
				Integer toto(0) ;
				for (size_t i = 0 ; i< C.rowdim() ; ++i)
					for (size_t j = 0 ; j< C.coldim() ; ++j) {
						FLINT::fmpz_get_mpz(toto.get_mpz(), fmpz_mat_entry(Cf,i,j));
						C.setEntry(i,j, toto );
					}
			}


			template<class _ZZ>
			void setFlintMatrix(FLINT::fmpz_mat_t Cf, const BlasSubmatrix<_ZZ> &C )
			{
				Integer toto(0) ;
				for (size_t i = 0 ; i< C.rowdim() ; ++i)
					for (size_t j = 0 ; j< C.coldim() ; ++j) {
						FLINT::fmpz_set_mpz(fmpz_mat_entry(Cf,i,j), C.getEntry(i,j).get_mpz_const());
					}
			}
#endif

		}

		template<class DenseIntMat>
		DenseIntMat &
		mul (DenseIntMat& C,
		     const DenseIntMat& A,
		     const DenseIntMat& B,
		     const mulMethod::FLINT & )
		{
			FLINT::fmpz_mat_t Af ;
			FLINT::fmpz_mat_t Bf ;
			FLINT::fmpz_mat_t Cf ;
			FLINT::fmpz_mat_init(Af, A.rowdim(), A.coldim());
			FLINT::fmpz_mat_init(Bf, B.rowdim(), B.coldim());
			FLINT::fmpz_mat_init(Cf, C.rowdim(), C.coldim());
			Protected::setFlintMatrix<typename DenseIntMat::Field>(Af,A);
			Protected::setFlintMatrix<typename DenseIntMat::Field>(Bf,B);
			// Timer Tim ;Tim.clear(); Tim.start();
			FLINT::fmpz_mat_mul(Cf,Af,Bf);
			// Tim.stop() ; std::cout << "inside : " << Tim << std::endl;
			Protected::getFlintMatrix<typename DenseIntMat::Field>(C,Cf);
			FLINT::fmpz_mat_clear(Af);
			FLINT::fmpz_mat_clear(Bf);
			FLINT::fmpz_mat_clear(Cf);

			return C;
		}
	} // BLAS3
}


#endif //  __LINBOX_matrix_blas3_mul_flint_INL

//Local Variables:
//mode: C++
//tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


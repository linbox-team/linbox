/* Copyright (C) 2010 LinBox
 * Written by <brice.boyer@imag.fr>
 *
 *
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_dense_nullspace_INL
#define __LINBOX_dense_nullspace_INL

#include "linbox/matrix/blas-matrix.h"

#include "fflas-ffpack/ffpack/ffpack.h" // LU
#include "fflas-ffpack/fflas/fflas.h" // trsm

#include <fflas-ffpack/utils/Matio.h> // write_field ;
#include <iostream>
#include <cassert>

namespace LinBox
{

	namespace Protected {

		/**
		 * @brief
		 *
		 * @param F
		 * @param Z
		 * @param ldZ
		 * @param lig1
		 * @param col1
		 * @param lig2
		 * @param col2
		 */
		template <class Field >
		void Zero(const Field & F,
			  typename Field::Element * Z, const size_t ldZ,
			  const size_t lig1, const size_t col1,
			  const size_t lig2, const size_t col2)
		{
			assert(lig1<lig2);
			assert(col1<col2);
			assert(col2<=ldZ);

			//! @todo use fzero
			for (size_t i = lig1 ; i < lig2 ; ++i)
				for (size_t j = col1; j < col2 ; ++j) // F.assign(*(Id+i*ldI+j),zero)
					*(Z+i*ldZ+j) = F.zero ;
			return;
		}


		/*! Creates identity matrix in \p F of size \p dim1 \p x \p dim2.
		 * @warning diag_num peut être < 0 !
		 * @bug     long et size_t ne cohabitent pas bien.
		 */
		template <class Field >
		void Identity(const Field & F,
			      typename Field::Element * Id, const size_t ldI,
			      const size_t lig1, const size_t col1,
			      const size_t lig2, const size_t col2)
		{
			assert(lig1<lig2);
			assert(col1<col2);
			assert(col2<=ldI);
			for (size_t i = lig1 ; i < lig2 ; ++i)
				for (size_t j = col1; j < col2 ; ++j) // F.assign(*(Id+i*ldI+j),zero)
					*(Id+i*ldI+j) = F.zero ;
			// Zero(F,Id,ldI,lig1,col1,lig2,col2);
			size_t nb_un = std::min(col2-col1,lig2-lig1)-1;

			typename Field::Element * Un_ici = Id+lig1*ldI+col1 ;
			for (size_t i = 0 ; i < nb_un ; ++i){
				*(Un_ici) = F.one ;
				Un_ici += ldI ;
				++Un_ici;
			}
			*(Un_ici) = F.one ;
			return;
		}


		/*!
		 * @brief The right or left nullspace (kernel or cokernel) of a matrix A
		 * We use the LU decomposition
		 * @param F the field in which \p A lives
		 * @param A is a matrix whose nullspace we look for.
		 * @param M number of lines in \p A
		 * @param N number of column of \p A
		 * @param lda the leading dimension of matrix \p A
		 * @param ker_dim the dimension of the kernel
		 * @tparam Field -
		 * @return a matrix of leading dimension ker_dim whose column vectors span the nullspace of A. Returns \p NULL (and not \f$\mathbf{0}\f$) if <code> ker_dim == 0 </code>.
		 */
		template<class Field>
		typename Field::Element *
		RightNullspaceDirect ( const Field & F,		// in place
				       typename Field::Element * A,
				       const size_t & M, // rows
				       const size_t & N, // colm
				       const size_t & lda, // leading dimension
				       //typename  Field::Element &* V ,
				       size_t & ker_dim )
		{
			typename Field::Element * NS;
			FFPACK::NullSpaceBasis((typename Field::Father_t)F, FFLAS::FflasRight, M, N, A, lda, NS, ker_dim);
			return NS ;

		}

		template<class Field>
		typename Field::Element *
		RightNullspaceIndirect ( const Field & F,		// in place
					 typename Field::Element * A,
					 const size_t & M, // rows
					 const size_t & N, // colm
					 const size_t & lda, // leading dimension
					 //typename  Field::Element &* V ,
					 size_t & ker_dim )
		{
			size_t *P = new size_t[M];
			size_t *Qt = new size_t[N];
			size_t R = FFPACK::LUdivine ((typename Field::Father_t)F, FFLAS::FflasNonUnit, FFLAS::FflasTrans,
						     M, N, A, lda, P, Qt, FFPACK::FfpackLQUP);


			// write_field (F, std::cout<<"ALU :="<<std::endl, A, M, N, N, true);
			// PrintLQUP (F,FFLAS::FflasNonUnit,FFLAS::FflasTrans,M,N,A,R,std::cout<<"L.Q.U.P"<<std::endl,Qt,P,true);

			// std::cout << "ker_dim = " << ker_dim << std::endl;

			delete[] P;

			ker_dim = N-R ;								// dimension of kernel
			if (ker_dim == 0) {
				delete[] Qt ;
				return NULL ;							// only 0 in kernel
			}
			size_t ldV = ker_dim ;
			typename Field::Element * V = new typename Field::Element[ker_dim*N];	// Result here.
			if (R == 0) {
				delete[] Qt ;
				Identity(F,V,ker_dim,0,0,N,ker_dim);
				return V ;
			}

			//////////////////////////

			Zero    (F,V,ldV,0,0,R,ker_dim);
			Identity(F,V,ldV,R,0,N,ker_dim);
			// write_field (F, std::cout<<"V init   ="<<std::endl, V, N, ker_dim, ker_dim,true);
			FFPACK::applyP((typename Field::Father_t)F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
				       ker_dim, 0, (int)R, V, ldV, Qt);			//  V = Q V

			// write_field (F, std::cout<<"V reordered   ="<<std::endl, V, N,ker_dim, ker_dim,true);

			// actually we just select a line in the inverse of L.
			//size_t wda = M ;


			if (N <= M) { // on a de la place...
				for ( size_t i=0; i< M; ++i ) {
					size_t j = 0 ;
					for ( ; j<std::min(i,N); ++j )
						*(A+i*N+j) = F.zero;
					if (i==j)
						*(A+i*N+j) = F.zero;
				}

				//FFPACK::applyP((typename Field::Father_t) F, FFLAS::FflasRight, FFLAS::FflasNoTrans, M,0,M, A, N, Q );
				FFPACK::applyP((typename Field::Father_t) F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
					       N,0,(int)N, A, N, Qt );
				for ( size_t i=0; i< N; ++i )
					*(A+N*i+i) = F.one ;

				//write_field (F, std::cout<<"A avant trsm   ="<<std::endl, A, M, N, N,true);
				FFLAS::ftrsm((typename Field::Father_t)F, FFLAS::FflasLeft, FFLAS::FflasUpper,
					     FFLAS::FflasNoTrans, FFLAS::FflasUnit,
					     N, ker_dim, F.one,  A, lda , V, ldV) ; // V = inv(Lower) V ;
				//write_field (F, std::cout<<"V if after trsm  ="<<std::endl, V, N,ker_dim, ker_dim,true);
			}
			else { // N > M we can't ftrsm because we can't add 0's to the lower part...
				typename Field::Element * L = new typename Field::Element[N*N];					// L_inf
				// début de L
				size_t i = 0 ;
				for ( ; i< M; ++i ){
					size_t j=0;
					for (; j< std::min(i,N) ; ++j )
						*(L+i*N+j) = F.zero ;
					if (i==j) {
						*(L+i*N+j) = F.zero ;
						j++;
						for (; j<N; ++j )
							*(L+i*N+j) = *(A+N*i+j);
					}
				}
				for ( ; i< N; ++i )
					for (size_t j = 0 ; j<N; ++j )
						*(L+i*N+j) = F.zero ;

				FFPACK::applyP((typename Field::Father_t) F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
					       N,0,(int)R, L, N, Qt );
				for ( size_t ii=0; ii< N; ++ii )
					*(L+N*ii+ii) = F.one ;
				// fin de L.
				//write_field (F, std::cout<<"U_1="<<std::endl, L, M, M, M,true);
				FFLAS::ftrsm((typename Field::Father_t)F, FFLAS::FflasLeft, FFLAS::FflasUpper,
					     FFLAS::FflasNoTrans, FFLAS::FflasUnit,
					     N, ker_dim, F.one,  L, N , V, ldV); 	// V = inv(Lower) V ;

				// write_field (F, std::cout<<"V else after trsm   ="<<std::endl, V, N,ker_dim, ker_dim,true);

				delete[] L;
			}
			//write_field (F, std::cout<<"V final   ="<<std::endl, V, ker_dim, M, M,true);
			delete[] Qt ;
			return V;

		}

		template<class Field>
		typename Field::Element *
		LeftNullspaceIndirect ( const Field & F,		// in place
					typename Field::Element * A,
					const size_t & M, // rows
					const size_t & N, // colm
					const size_t & lda, // leading dimension
					//typename  Field::Element &* V ,
					size_t & coker_dim )
		{
			typename Field::Element * NS;
			FFPACK::NullSpaceBasis((typename Field::Father_t)F, FFLAS::FflasLeft, M, N, A, lda, NS, coker_dim);
			return NS ;
		}

		// directement à gauche (X LQUP noyau de U puis invQ invL...)
		template<class Field>
		typename Field::Element *
		LeftNullspaceDirect ( const Field & F,		// in place
				      typename Field::Element * A,
				      const size_t & M, // rows
				      const size_t & N, // colm
				      const size_t & lda, // leading dimension
				      size_t & coker_dim)
		{

			size_t *P = new size_t[N];
			size_t *Q = new size_t[M];

			//write_field (F, std::cout<<"A avant LU   ="<<std::endl, A, M, N, N, true);
			size_t R = FFPACK::LUdivine ((typename Field::Father_t)F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans,
						     M, N, A, lda, P, Q, FFPACK::FfpackLQUP);
			assert(R<=std::min(M,N));

			//write_field (F, std::cout<<"ALU :="<<std::endl, A, M, N, N, true);
			//PrintLapackPermutation(P,N,std::cout<<"Permutation P := ");
			//PrintPermutation(F , FFLAS::FflasNoTrans, P, N, 0, N, std::cout, true);
			//PrintLapackPermutation(Q,M,std::cout<<"Permutation Q := ");
			//PrintPermutation(F , FFLAS::FflasTrans, Q, M, 0, M, std::cout, true);
			//PrintLQUP (F,FFLAS::FflasUnit,FFLAS::FflasNoTrans,M,N,A,R,std::cout<<"L.Q.U.P"<<std::endl,Q,P,true);
			coker_dim = M -R ; // dimension of co-kernel.
			//std::cout << "coker_dim = " << coker_dim << std::endl;
			delete[] P; // on s'en fout de P !
			if (coker_dim == 0){
				delete[] Q ;
				return NULL;	// CoKernel is \f$\{\mathbf{0}_n\}\f$
			}


			size_t ldV = M ;
			typename Field::Element * V = new typename Field::Element[coker_dim * M]; // le résultat sera ici.
			if (R == 0) {
				delete[] Q ;
				Identity(F,V,ldV,0,0,coker_dim,M);
				return V ;
			}
			Zero    (F,V,ldV,0,0,coker_dim,R);
			Identity(F,V,ldV,0,R,coker_dim,M);
			// write_field (F, std::cout<<"V init   ="<<std::endl, V, coker_dim, M, M,true);
			FFPACK::applyP((typename Field::Father_t)F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				       coker_dim, 0, (int)R, V, ldV, Q); // V = V  tQ

			//write_field (F, std::cout<<"V reordered   ="<<std::endl, V, coker_dim, M, M,true);

			//size_t wda = M ;
			if (M <= N) {
				for ( size_t i=0; i< M; ++i )
					for (size_t j = i ; j<N; ++j )
						*(A+i*N+j) = F.zero;

				FFPACK::applyP((typename Field::Father_t) F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
					       M,0,(int)M, A, N, Q );
				for ( size_t i=0; i< M; ++i )
					*(A+N*i+i) = F.one ;

				//write_field (F, std::cout<<"A avant trsm   ="<<std::endl, A, M, N, N,true);
				FFLAS::ftrsm((typename Field::Father_t)F, FFLAS::FflasRight, FFLAS::FflasLower,
					     FFLAS::FflasNoTrans, FFLAS::FflasUnit,
					     coker_dim , M , F.one,  A, lda , V, ldV) ; // V = V inv(Lower) ;
			}
			else { // M > N we can't ftrsm because we can't add 0's to the lower part...
				typename Field::Element * L = new typename Field::Element[M*M]; // L_inf

				for ( size_t i=0; i< M; ++i ){ // copying A_inf to L
					size_t j=0;
					for (; j< std::min(i,N) ; ++j )
						*(L+i*M+j) = *(A+N*i+j);
					for (; j<M; ++j )
						*(L+i*M+j) = F.zero;
				}
				FFPACK::applyP((typename Field::Father_t) F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
					       M,0,(int)R, L, M, Q );
				for ( size_t i=0; i< M; ++i )
					*(L+M*i+i) = F.one ;

				//write_field (F, std::cout<<"U_1="<<std::endl, L, M, M, M,true);
				FFLAS::ftrsm((typename Field::Father_t)F, FFLAS::FflasRight, FFLAS::FflasLower,
					     FFLAS::FflasNoTrans, FFLAS::FflasUnit,
					     coker_dim , M , F.one,  L, M , V, ldV) ; // V = V inv(Lower) ;

				delete[] L;
			}
			//write_field (F, std::cout<<"V final   ="<<std::endl, V, coker_dim, M, M,true);
			delete[] Q ;
			return V;
		}
	} // Protected

		/** Computes the kernel of a dense matrix using \c LQUP.
		 *
		 * Acccording to the dimensions of the input matrix, we chose different methods.
		 * @warning timings may vary and these choices were made on an experimental basis.
		 *
		 * @param F  Field
		 * @param Side  left or right from \c LinBox::SideTag
		 * @param m rows
		 * @param n cols
		 * @param A input matrix
		 * @param lda leading dimension of A
		 * @param Ker Kernel. \c NULL if \c kerdim==0
		 * @param ldk leading dimension of the kernel.
		 * @param kerdim dimension of the kernel.
		 * @return dimension of the kernel.
		 *
		 * @warning A is modified.
		 */
	template<class Field>
	size_t
	NullSpaceBasisIn ( const Field & F, const Tag::Side Side,
			const size_t & m, const size_t & n,
			typename Field::Element * A, const size_t & lda,
			typename Field::Element *& Ker, size_t& ldk,
			size_t & kerdim)
	{
#if 0
		if (Side == Tag::Side::Right){
			if (m < n)
				Ker = RightNullspaceDirect(F,A,m,n,lda,kerdim) ;
			else
				Ker = RightNullspaceIndirect(F,A,m,n,lda,kerdim) ;
			ldk = kerdim;

		}
		else {
			if (m < n)
				Ker = LeftNullspaceDirect(F,A,m,n,lda,kerdim) ;
			else
				Ker = LeftNullspaceIndirect(F,A,m,n,lda,kerdim) ;
			ldk = m;
		}
#else

		FFPACK::NullSpaceBasis ((typename Field::Father_t) F, (FFLAS::FFLAS_SIDE) Side,
						 m,n, A, lda, Ker, ldk, kerdim);
#endif
		return kerdim;
	}


	//!@todo uses too much memory
	template<class DenseMat>
	size_t&
	NullSpaceBasisIn (const Tag::Side Side,
			BlasSubmatrix<DenseMat> & A,
			BlasMatrix<typename DenseMat::Field> & Ker,
			size_t & kerdim)
	{

		typedef typename DenseMat::Field Field;

		typename Field::Element * Ker_ptr;
		size_t ldk;
		NullSpaceBasisIn(A.field(),Side,A.rowdim(),A.coldim(), A.getWritePointer(),A.getStride(), Ker_ptr,ldk,kerdim);
		if (Side == Tag::Side::Right){
			Ker.resize(A.coldim(),kerdim);
		}
		else {
			assert(Side == Tag::Side::Left);
			Ker.resize(kerdim,A.rowdim());
		}
		//! @todo use copy
		const typename Field::Element * Ker_ptri = Ker_ptr;
		for(typename BlasMatrix<Field>::Iterator it=Ker.Begin(); it!= Ker.End(); ++it,++Ker_ptri)
			A.field().assign(*it,*Ker_ptri);

		delete[] Ker_ptr ;

		return kerdim;
	}

	template<class Field>
	size_t&
	NullSpaceBasisIn (const Tag::Side Side,
			BlasMatrix<Field> & A,
			BlasMatrix<Field> & Ker,
			size_t & kerdim)
	{
		BlasSubmatrix< BlasMatrix<Field>  > Asub(A);
		return NullSpaceBasisIn(Side,Asub,Ker,kerdim);
	}

	template<class Field>
	size_t&
	NullSpaceBasis (const Tag::Side Side,
			const BlasMatrix<Field> & A,
			BlasMatrix<Field> & Ker,
			size_t & kerdim)
	{
		const BlasSubmatrix<BlasMatrix<Field> > Asub (A);
		return NullSpaceBasis<Field>(Side,Asub,Ker,kerdim);
	}

	template<class DenseMat>
	size_t&
	NullSpaceBasis (const Tag::Side Side,
			const BlasSubmatrix<DenseMat> & A,
			BlasMatrix<typename DenseMat::Field> & Ker,
			size_t & kerdim)
	{
		BlasMatrix<typename DenseMat::Field> B (A);
		return NullSpaceBasisIn<typename DenseMat::Field>(Side,B,Ker,kerdim);
	}


} // LinBox

#endif // __LINBOX_dense_nullspace_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

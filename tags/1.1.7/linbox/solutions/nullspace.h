/* Copyright (C) 2009 LinBox
 * Written by <brice.boyer@imag.fr>
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_modulardense_nullspace_H
#define __LINBOX_modulardense_nullspace_H

/** \file linbox/solutions/nullspace.h
 * @brief The right or left nullspace (kernel or cokernel) of a matrix A
 * @details Provides :
 *	- the nullspace of a matrix \p A
 *	- (soon) a random vector within the nullspace of \p A
 */

//#include "vector-traits.h"
//#include "linbox/blackbox/blas-blackbox.h"
#include "linbox/matrix/blas-matrix.h"

#include "linbox/ffpack/ffpack.h" // LU
#include "linbox/fflas/fflas.h" // trsm

//#include "linbox/blackbox/dense.h"
//#include "linbox/blackbox/submatrix.h"

//#include "Matio.h" // write_field ;
#include <iostream>
#include <cassert>
namespace LinBox
{

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
		typename Field::Element zero;
		F.init(zero,0UL);

		for (size_t i = lig1 ; i < lig2 ; ++i)
			for (size_t j = col1; j < col2 ; ++j) // F.assign(*(Id+i*ldI+j),zero) 
				*(Z+i*ldZ+j) = zero ;
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
		typename Field::Element one,zero;
		F.init(one,1UL);
		F.init(zero,0UL);
		for (size_t i = lig1 ; i < lig2 ; ++i)
			for (size_t j = col1; j < col2 ; ++j) // F.assign(*(Id+i*ldI+j),zero) 
				*(Id+i*ldI+j) = zero ;
		//        Zero(F,Id,ldI,lig1,col1,lig2,col2);
		size_t nb_un = std::min(col2-col1,lig2-lig1)-1;

		typename Field::Element * Un_ici = Id+lig1*ldI+col1 ; 
		for (size_t i = 0 ; i < nb_un ; ++i){
			*(Un_ici) = one ;
			Un_ici += ldI ;
			++Un_ici;
		}
		*(Un_ici) = one ;
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
		size_t *P = new size_t[N];
		size_t *Qt = new size_t[M];
		size_t R = FFPACK::LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans,
					     M, N, A, lda, P, Qt, FFPACK::FfpackLQUP);
		delete[] Qt;

		ker_dim = N-R ;								// dimension of kernel
		if (ker_dim == 0) {
			delete[] P ;
			return NULL ;							// only 0 in kernel
		}
		typename Field::Element * V = new typename Field::Element[ker_dim*N];	// Result here.
		size_t ldV = ker_dim ;
		if (R == 0) {
			delete[] P ;
			Identity(F,V,ldV,0,0,N,ker_dim);
			return V ;
		}

		for (size_t i = 0 ; i < R ; ++i)					// V <- U2 
			FFLAS::fcopy (F, ker_dim, V + i * ldV, 1, A + R + i*lda, 1);
		typename Field::Element one ;
		F.init(one,1UL);
		FFLAS::ftrsm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
			     R, ker_dim, one,  A, lda , V, ldV) ;			// (V = U2, A = U1) : U2 <- inv(U1)*U2 ;

		typename Field::Element minus_one, zero;
		F.init(zero,0UL);
		F.neg(minus_one, one);
		for ( size_t i = R ; i < N; ++i){					// filling the rest of V with -Id
			for (size_t j = 0 ; j < ker_dim ; ++j)
				*(V+ker_dim*i+j) = zero ;
			*(V+ker_dim*i+i-R) = minus_one ;
			/* // FIXME : faster ?
			   for (size_t j = 0 ; j < ker_dim ; ++j) {

			   if (i-R == j)
			 *(V+ker_dim*i+j) = minus_one ;
			 else
			 *(V+ker_dim*i+j) = zero ;
			 }
			 */

		}
		FFPACK::applyP(F, FFLAS::FflasLeft, FFLAS::FflasTrans,
			       ker_dim, 0, R, V, ldV , P);				//  X = P^{-1}V
		delete[] P;
		return V;

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
		size_t R = FFPACK::LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasTrans,
					     M, N, A, lda, P, Qt, FFPACK::FfpackLQUP);


		//                write_field (F, std::cout<<"ALU :="<<std::endl, A, M, N, N, true);
		//                PrintLQUP (F,FFLAS::FflasNonUnit,FFLAS::FflasTrans,M,N,A,R,std::cout<<"L.Q.U.P"<<std::endl,Qt,P,true);

		//                std::cout << "ker_dim = " << ker_dim << std::endl;

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
		//                                write_field (F, std::cout<<"V init   ="<<std::endl, V, N, ker_dim, ker_dim,true);
		FFPACK::applyP(F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
			       ker_dim, 0, R, V, ldV, Qt);			//  V = Q V

		//                write_field (F, std::cout<<"V reordered   ="<<std::endl, V, N,ker_dim, ker_dim,true);

		// actually we just select a line in the inverse of L.
		//size_t wda = M ;

		typename Field::Element one, zero;
		F.init(zero,0UL);
		F.init(one,1UL);
		//                F.neg(minus_one, one);

		if (N <= M) { // on a de la place...
			for ( size_t i=0; i< M; ++i ) {
				size_t j = 0 ;
				for ( ; j<std::min(i,N); ++j )
					*(A+i*N+j) = zero;
				if (i==j)
					*(A+i*N+j) = zero;
			}

			//FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans, M,0,M, A, N, Q );
			FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N,0,N, A, N, Qt );
			for ( size_t i=0; i< N; ++i )
				*(A+N*i+i) = one ;

			//write_field (F, std::cout<<"A avant trsm   ="<<std::endl, A, M, N, N,true);
			FFLAS::ftrsm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasUnit,
				     N, ker_dim, one,  A, lda , V, ldV) ;						// V = inv(Lower) V ;
                         //write_field (F, std::cout<<"V if after trsm  ="<<std::endl, V, N,ker_dim, ker_dim,true);
		} else { // N > M we can't ftrsm because we can't add 0's to the lower part...
			typename Field::Element * L = new typename Field::Element[N*N];					// L_inf
			// début de L
			size_t i = 0 ;
			for ( ; i< M; ++i ){
				size_t j=0;
				for (; j< std::min(i,N) ; ++j )
					*(L+i*N+j) = zero ;
				if (i==j) {
					*(L+i*N+j) = zero ;
					j++;
					for (; j<N; ++j )
						*(L+i*N+j) = *(A+N*i+j);
				}
			}
			for ( ; i< N; ++i )
				for (size_t j = 0 ; j<N; ++j )
					*(L+i*N+j) = zero ;

			FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
					N,0,R, L, N, Qt );
			for ( size_t i=0; i< N; ++i )
				*(L+N*i+i) = one ;
			// fin de L.
			//write_field (F, std::cout<<"U_1="<<std::endl, L, M, M, M,true);
			FFLAS::ftrsm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasUnit,
				     N, ker_dim, one,  L, N , V, ldV) ;						// V = inv(Lower) V ;

			//                        write_field (F, std::cout<<"V else after trsm   ="<<std::endl, V, N,ker_dim, ker_dim,true);

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
		size_t *P = new size_t[M];
		size_t *Qt = new size_t[N];
		size_t R = FFPACK::LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasTrans,
					     M, N, A, lda, P, Qt, FFPACK::FfpackLQUP);
		delete[] Qt;

		coker_dim = M-R ;								// dimension of kernel
		if (coker_dim == 0) {
			delete[] P ;
			return NULL ;							// only 0 in kernel
		}
		size_t ldV = M ;
		typename Field::Element * V = new typename Field::Element[coker_dim*M];	// Result here.
		if (R == 0) {
			delete[] P ;
			Identity(F,V,ldV,0,0,coker_dim,M);
			return V ;
		}

		for (size_t i = 0 ; i < coker_dim ; ++i)					// copy U2 to result V before updating with U1
			FFPACK::fcopy (F, R, V + i * ldV, 1, A + (R + i)*lda, 1);
		typename Field::Element one, minus_one ;
		F.init(one,1UL);
		F.neg(minus_one, one);
		FFLAS::ftrsm(F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
			     coker_dim, R, minus_one,  A, lda , V, ldV) ;			// V = U2  ; A = U1 ; U2 <- inv(U1)*U2 ;

		typename Field::Element zero;
		F.init(zero,0UL);
		Identity(F,V,ldV,0,R,coker_dim,M);
		//                for ( size_t i = 0 ; i < coker_dim ; ++i){					// filling the rest of V with minus identity
		//                        for (size_t j = R ; j < M ; ++j)
		//                                *(V+i*M+j) = zero ;
		//                        *(V+i*M+i+R) = one ;
		/* // FIXME : faster ?
		   for (size_t j = R ; j < M ; ++j) {
		   if (j == i+R)
		 *(V+M*i+j) = one ;
		 else
		 *(V+M*i+j) = zero ;
		 }
		 */

		//                }
		FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
			       coker_dim, 0, R, V, ldV , P);				//  X = P^{-1}V
		delete[] P;
		return V;

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
		size_t R = FFPACK::LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans,
					     M, N, A, lda, P, Q, FFPACK::FfpackLQUP);
		assert(R<=std::min(M,N));

		//write_field (F, std::cout<<"ALU :="<<std::endl, A, M, N, N, true);
		//PrintLapackPermutation(P,N,std::cout<<"Permutation P := ");
		//PrintPermutation(F , FFLAS::FflasNoTrans, P, N, 0, N, std::cout, true);
		//PrintLapackPermutation(Q,M,std::cout<<"Permutation Q := ");
		//PrintPermutation(F , FFLAS::FflasTrans, Q, M, 0, M, std::cout, true);
		//PrintLQUP (F,FFLAS::FflasUnit,FFLAS::FflasNoTrans,M,N,A,R,std::cout<<"L.Q.U.P"<<std::endl,Q,P,true);
		coker_dim = M - R ;											// dimension of co-kernel.
		//std::cout << "coker_dim = " << coker_dim << std::endl;
		delete[] P;												// on s'en fout de P !
		if (coker_dim == 0){
			delete[] Q ;
			return NULL;											// CoKernel is \f$\{\mathbf{0}_n\}\f$
		}

		typename Field::Element one, zero ;									// 1,0 dans le corps
		F.init(one,1UL);
		F.init(zero,0UL);

		size_t ldV = M ;
		typename Field::Element * V = new typename Field::Element[coker_dim * M]; // le résultat sera ici.
		if (R == 0) {
			delete[] Q ;
			Identity(F,V,ldV,0,0,coker_dim,M);
			return V ;
		}
		Zero    (F,V,ldV,0,0,coker_dim,R);
		Identity(F,V,ldV,0,R,coker_dim,M);
		//                write_field (F, std::cout<<"V init   ="<<std::endl, V, coker_dim, M, M,true);
		FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
			       coker_dim, 0, R, V, ldV, Q);							//  V = V  tQ FIXME BUG XXX WARNING TODO pourquoi 0..R ??

		//write_field (F, std::cout<<"V reordered   ="<<std::endl, V, coker_dim, M, M,true);

		// actually we just select a line in the inverse of L.
		//size_t wda = M ;
		if (M <= N) {
			for ( size_t i=0; i< M; ++i )
				for (size_t j = i ; j<N; ++j )
					*(A+i*N+j) = zero;

			FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans, M,0,M, A, N, Q );
			for ( size_t i=0; i< M; ++i )
				*(A+N*i+i) = one ;

			//write_field (F, std::cout<<"A avant trsm   ="<<std::endl, A, M, N, N,true);
			FFLAS::ftrsm(F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit,
				     coker_dim , M , one,  A, lda , V, ldV) ;						// V = V inv(Lower) ;
		} else { // M > N we can't ftrsm because we can't add 0's to the lower part...
			typename Field::Element * L = new typename Field::Element[M*M];					// L_inf

			for ( size_t i=0; i< M; ++i ){ // copying A_inf to L
				size_t j=0;
				for (; j< std::min(i,N) ; ++j )
					*(L+i*M+j) = *(A+N*i+j);
				for (; j<M; ++j )
					*(L+i*M+j) = zero;
			}
			FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans, M,0,R, L, M, Q );
			for ( size_t i=0; i< M; ++i )
				*(L+M*i+i) = one ;

			//write_field (F, std::cout<<"U_1="<<std::endl, L, M, M, M,true);
			FFLAS::ftrsm(F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit,
				     coker_dim , M , one,  L, M , V, ldV) ;						// V = V inv(Lower) ;

			delete[] L;
		}
		//write_field (F, std::cout<<"V final   ="<<std::endl, V, coker_dim, M, M,true);
		delete[] Q ;
		return V;
	} 

	/** 
	 * @brief Computes the kernel of a Dense matrix using \c LQUP.
	 *
	 * acccording to the dimensions of the input matrix, we chose different methods.
	 * @warning timings may vary (or may not) and these choices were made on an experimental basis.
	 * 
	 * @param F 
	 * @param Side  left or right from \p FFLAS::FFLAS_SIDE
	 * @param m rows
	 * @param n cols
	 * @param A input matrix
	 * @param lda leading dimension of A
	 * @param Ker Kernel. \c NULL if \c kerdim==0
	 * @param ldk leading dimension of the kernel.
	 * @param kerdim dimension of the kernel.
	 * @return dimension of the kernel. 
	 */
	template<class Field>
	size_t 
	NullSpaceBasis (const Field& F, const FFLAS::FFLAS_SIDE Side,
			const size_t & m, const size_t & n,
			typename Field::Element * A, const size_t & lda,
			typename Field::Element *& Ker, size_t& ldk,
			size_t & kerdim)
	{
		if (Side == FFLAS::FflasRight){
			if (m < n) 
				Ker = RightNullspaceDirect(F,A,m,n,lda,kerdim) ;
			 else
				Ker = RightNullspaceIndirect(F,A,m,n,lda,kerdim) ;
			ldk = kerdim;

		} else {
			if (m < n) 
				Ker = LeftNullspaceDirect(F,A,m,n,lda,kerdim) ;
			 else 
				Ker = LeftNullspaceIndirect(F,A,m,n,lda,kerdim) ;
			ldk = m;
		}
		return kerdim;
	}

	template<class Field>
	size_t& 
	NullSpaceBasis (const Field& F, const FFLAS::FFLAS_SIDE Side,
			const BlasMatrix<typename Field::Element> & A,
			BlasMatrix<typename Field::Element> & Ker, 
                        size_t & kerdim) {

            typename Field::Element * Kert;
            size_t ldk;
            NullSpaceBasis(F,Side,A.rowdim(),A.coldim(), A.getPointer(),A.getStride(), Kert,ldk,kerdim);
            if (Side == FFLAS::FflasRight){
                Ker = BlasMatrix<typename Field::Element>(A.rowdim(),kerdim);
            } else {
                Ker = BlasMatrix<typename Field::Element>(kerdim,A.coldim());
            }
            for(typename BlasMatrix<typename Field::Element>::RawIterator it=Ker.rawBegin(); it!= Ker.rawEnd(); ++it,++Kert)
                *it=*Kert;
           return kerdim;
	}
} // LinBox

#endif // __LINBOX_modulardense_nullspace_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

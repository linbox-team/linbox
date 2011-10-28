/* Copyright (C) 2011 LinBox
 * Written by Brice Boyer <brice.boyer@imag.fr>
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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __LINBOX_matrix_random_matrix_INL
#define __LINBOX_matrix_random_matrix_INL
/** @file matrix/random-matrix.inl
 * @ingroup matrix
 * @brief Implementation of random matrices.
 *
 */

namespace LinBox
{

	//! This is the namespace all LinBox internal code is in.
	namespace Protected {

		template<class Randiter,class Field>
		BlasBlackbox<Field> &
		random_lu_rank(const Field			& F,
			       const Randiter                   & R,
			       BlasBlackbox<Field>              & A,
			       int                              & rank,
			       const RingCategories::ModularTag & tag)
		{
			size_t m = A.rowdim() ;
			size_t n = A.coldim() ;

			linbox_check(m != 0);
			linbox_check(n != 0);
			if (rank == -1)
				rank = (int) ( (double)std::min(m,n)*drand48() );
			linbox_check(!(rank<0 || rank>(int)std::min(m,n)));

			typedef typename Field::Element Element ;

			// be ready for PLUQ
			//size_t  * P = new size_t [m]   ;
			BlasPermutation<size_t> P(m);
			//Element * L = new Element[m*m] ;
			//TriangularBlasMatrix<Element> L(m,m,LinBoxTag::Lower,LinBoxTag::Unit);
			//! @todo !!!
			BlasMatrix<Element> L(m,m);
			// Element * U = new Element[m*n] ;
			// TriangularBlasMatrix<Element> U(m,n,LinBoxTag::Upper,LinBoxTag::NonUnit);
			//
			BlasMatrix<Element> U(m,n);
			//size_t  * Q = new size_t [n]   ;
			BlasPermutation<size_t> Q(n);

			// be ready for random elements
			NonzeroRandIter<Field> Rnz(F,R);
			Element one,zero;
			F.init(one,1UL);
			F.init(zero,0UL);

			/* Create L a random invertible lower unit triangular matrix (m x m format) */
			for (size_t j=0 ; j<m ; ++j)
				for (size_t i=j+1; i<m;++i)
					R.random( L.refEntry( i,j ) );
#if 1
			for (size_t i = 0; i < (size_t) m; ++i)
				Rnz.random( L.refEntry( i,i ) ); // non zero diagonal
			for (size_t j=0 ; j<m ; ++j)
				for (size_t i=0; i<j;++i)
					L.setEntry( i,j,zero );


#endif

			/* Create U a random rank-r upper non-unit triangular matrix (m x n format) */
			for (size_t i = 0 ; i < (size_t) rank; ++i)
				for (size_t j = i+1; j < n; ++j)
					R.random( U.refEntry( i,j ) ); // r random 'triangular' lines
			for (size_t i = 0; i < (size_t) rank; ++i)
				Rnz.random( U.refEntry( i,i ) ); // non zero diagonal on rank first lines
#if 1
			for (size_t i = rank ; i < m ; ++i)
				for (size_t j = i ; j < n ; ++j)
					U.setEntry( i,j,zero ) ; //  zero on remaining 'triangular' lines
#endif

			/**
			 * @todo RandomPermutation avec P de type [Matrix-Blas]Permutation
			 */
			RandomBlasPermutation(P);
			RandomBlasPermutation(Q);

			BlasMatrixDomain< Field > BMD(F) ;
			// LU
			/**
			 * @todo : L = [[L1,0],[A,L2]] ;U = [[U1,B],[0,U2]] ; LU = [[ rec(L1,U1), ftrmm(L1,B)],[ftrmm(A,U1),fgemm(A,B)+rec(L2,U2) ]]
			 * de même UL
			 */

			BlasMatrix<Element> A_ptr(A) ;

			BMD.mul(A_ptr,L,U);

			/*!
			 * @todo create BMD.applyP(A,P,LinBoxTag::Left) ;
			 * avec P : BlasPermutation
			 * ou P : MatrixPermutation
			 * @todo BlasPermutation a un ordre \p p et une taille \p r distinctes !!!
			 */

			BMD.mulin_left(A_ptr,Q);
			BMD.mulin_right(P,A_ptr);

			return A ;

		}


		//!@todo ZZ is A.field() !
		template<class Randiter, class Ring>
		BlasBlackbox<Ring> &
		random_lu_rank(const Ring          & ZZ,
			       const Randiter      & R,
			       BlasBlackbox<Ring>  &A,
			       int                 & rank,
			       const RingCategories::IntegerTag & tag)
		{

			typedef typename Ring::Element Int ;
			size_t m = A.rowdim() ;
			size_t n = A.coldim() ;

			linbox_check(m != 0);
			linbox_check(n != 0);
			if (rank == -1)
				rank = (int) ( (double)std::min(m,n)*drand48() );
			linbox_check(!(rank<0 || rank>(int)std::min(m,n)));


			// be ready for PLUQ
			BlasPermutation<size_t> P(m);
			BlasMatrix<Int> L(m,m);
			BlasMatrix<Int> U(m,n);
			BlasPermutation<size_t> Q(n);

			// be ready for random elements
			Randiter S_(R);
			S_.setBits(R.getBits()-1);
			RandomIntegerIter<false> T_(3);
			NonzeroRandIter<Ring,RandomIntegerIter<false> > U_(ZZ,T_);

			Int one(1),zero(0);

			/* Create L a random invertible lower unit triangular matrix (m x m format) */
			for (size_t j=0 ; j<m ; ++j)
				for (size_t i=j+1; i<m;++i)
					S_.random( L.refEntry( i,j ) );
#if 1
			for (size_t i = 0; i < (size_t) m; ++i)
				U_.random( L.refEntry( i,i ) ); // non zero diagonal
			for (size_t j=0 ; j<m ; ++j)
				for (size_t i=0; i<j;++i)
					L.setEntry( i,j,zero );


#endif

			/* Create U a random rank-r upper non-unit triangular matrix (m x n format) */
			for (size_t i = 0 ; i < (size_t) rank; ++i)
				for (size_t j = i+1; j < n; ++j)
					T_.random( U.refEntry( i,j ) ); // r random 'triangular' lines
			for (size_t i = 0; i < (size_t) rank; ++i)
				U_.random( U.refEntry( i,i ) ); // non zero diagonal on rank first lines
#if 1
			for (size_t i = rank ; i < m ; ++i)
				for (size_t j = i ; j < n ; ++j)
					U.setEntry( i,j,zero ) ; //  zero on remaining 'triangular' lines
#endif

			RandomBlasPermutation(P);
			RandomBlasPermutation(Q);

			BlasMatrixDomain< Ring > BMD(ZZ) ;
			MatrixDomain< Ring > MD(ZZ) ;
			// LU
			// L.write(std::cout << "L:=",true ) << ';' << std::endl;
			// L.write(std::cout << "U:=",true ) << ';' << std::endl;

			typedef typename Ring::Element Element;
			BlasMatrix<Element> A_ptr(A) ;

			MD.mul(A_ptr,L,U);
			// A.write(std::cout << "pre A=",true) << std::endl;

			BMD.mulin_left(A_ptr,Q);
			BMD.mulin_right(P,A_ptr);

			// P.write(std::cout<<"P:=",false) << std::endl;
			// P.write(std::cout<<"Q:=",false) << std::endl;
			// P.write(std::cout<<"Q:=",true) << std::endl;

			return A ;
		}


		template<class Randiter,class Field>
		BlasMatrix<typename Field::Element> &
		random_lu_rank(const Field			   & F,
			       const Randiter                      & R,
			       BlasMatrix<typename Field::Element> & A,
			       int                                 & rank,
			       const RingCategories::IntegerTag    & tag)
		{
			BlasBlackbox<Field> Alink(F,A);
			random_lu_rank(F,R,Alink,rank,RingCategories::IntegerTag());
			return A ;
		}

#if 0 /*  BlasMatrix spec. */
		template<class Randiter, class Field>
		BlasBlackbox<Field> &
		random_rankupdate( const Field        & F
				   ,const Randiter    & R
				   ,BlasBlackbox<Field>& A
				   , int              & rank
				 )
		{
			size_t m = A.rowdim();
			size_t n = A.coldim();


			BlasBlackbox<Field> D(F,m,rank) ;
			BlasBlackbox<Field> G(F,rank,n) ;
			RandomBlasBlackbox<Randiter, Field> RandMatGen(F,R);
			RandMatGen.random(D) ;
			RandMatGen.random(G) ;
			MatrixDomain<Field> MD(F);
			MD.mul(A,D,G);
			return A ;
		}

		template<class Randiter>
		BlasBlackbox<PID_integer> &
		random_rankupdate( PID_integer               & F //!@bug const !
				   ,const Randiter           & R
				   ,BlasBlackbox<PID_integer> & A
				   , int                     & rank
				 )
		{
			//! @todo check randomness
			size_t m = A.rowdim();
			size_t n = A.coldim();

			BlasBlackbox<PID_integer> D(F,(size_t)m,(size_t)rank) ;
			BlasBlackbox<PID_integer> G(F,(size_t)rank,(size_t)n) ;
			Randiter S_(R);
			S_.setBits(R.getBits()-1);
			RandomBlasBlackbox<Randiter,PID_integer > RandMatGen(F,S_);
			RandMatGen.random(D) ;
			RandomIntegerIter<false> T_(3);
			RandomBlasBlackbox<RandomIntegerIter<false>,PID_integer > RandSmallMatGen(F,T_);
			RandMatGen.random(G) ;
			MatrixDomain<PID_integer> MD(F);
			MD.mul(A,D,G);
			return A ;
		}
#endif

		template<class Randiter,class Field>
		BlasMatrix<typename Field::Element> &
		random_rankupdate( Field                & F //!@bug const !
				   ,const Randiter      & R
				   ,BlasMatrix<typename Field::Element> & A
				   , int                & rank
				   , const RingCategories::IntegerTag          &tag
				 )
		{
			typedef typename Field::Element Int ;
			size_t m = A.rowdim();
			size_t n = A.coldim();

			BlasMatrix<Int> D((size_t)m,(size_t)rank) ;
			BlasMatrix<Int> G((size_t)rank,(size_t)n) ;
			Randiter S_(R);
			S_.setBits(R.getBits()-1);
			RandomDenseMatrix<Randiter,Field > RandMatGen(F,S_);
			RandMatGen.random(D) ;
			RandomIntegerIter<false> T_(3);
			RandomDenseMatrix<RandomIntegerIter<false>,Field > RandSmallMatGen(F,T_);
			RandMatGen.random(G) ;

			MatrixDomain<Field> MD(F);
			MD.mul(A,D,G);

			/// @bug do perms ?

#if 0 /*  necessary */
			BlasPermutation<size_t>  P(m);
			BlasPermutation<size_t>  Q(n);

			RandomBlasPermutation(P);
			RandomBlasPermutation(Q);

			PID_integer ZZ ;
			BlasMatrixDomain< PID_integer > BMD(ZZ) ;

			BMD.mulin_left (A,Q);
			BMD.mulin_right(P,A);
#endif

			return A ;
		}
	} // Protected

	//////// Random Matrix
	template<class Randiter, class Field>
	template<class Matrix>
	Matrix &
	RandomDenseMatrix<Randiter, Field>::random( Matrix & A)
	{
		for (size_t i = 0 ; i < A.rowdim() ; ++i)
			for (size_t j = 0 ; j < A.coldim() ; ++j)
				R_.random(A.refEntry(i,j));
		return A;
	}

	//////// Random Matrix with prescribed Rank
#if 0
	template<class Randiter, class Field>
	template<class Matrix,class Method>
	Matrix &
	RandomDenseMatrix<Randiter, Field>::randomRank(Matrix & A,
						       int    rank
						       , const Method & meth )
	{
		throw NotImplementedYet(__func__,__FILE__,__LINE__);
	}
#endif


	template<class Randiter, class Field>
	template<class Matrix>
	Matrix &
	RandomDenseMatrix<Randiter, Field>::randomRank(
						       Matrix & A
						       , int  rank
						      )
	{
		//! @todo use CatergoryTag
		return randomRank(A,rank,RankBuilder::LU_());
	}

	template<class Randiter, class Field>
	template< class Matrix>
	Matrix &
	RandomDenseMatrix<Randiter, Field>::randomRank(
						       Matrix  & A
						       , int   rank
						       , const RankBuilder::LU_ & meth
						      )
	{
		Protected::random_lu_rank( F_,R_,A,rank,
				typename FieldTraits<Field>::categoryTag());
		return A ;
	}



	template<class Randiter, class Field>
	template<class Matrix>
	Matrix &
	RandomDenseMatrix<Randiter, Field>::randomRank(
						       Matrix    & A
						       , int     rank
						       , const RankBuilder::Rank_update_ & meth
						      )
	{
		return Protected::random_rankupdate(F_,R_,A,rank,
						    typename FieldTraits<Field>::categoryTag());
	}

	/*  dense matrix with random entries */

	void RandomBlasPermutation(BlasPermutation<size_t> & P)
	{
		size_t * Pt = P.getWritePointer();
		// size_t n = P.getSize();
		size_t r = P.getOrder();
		size_t n = r ; // no size given ?
		for (size_t i = 0 ; i < r ; ++i) {
			Pt[i] = i + size_t(double (n-i)*( drand48() ) ) ;
		}
		return ;
	}


} // LinBox

#endif // __LINBOX_matrix_random_matrix_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


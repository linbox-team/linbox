/* Copyright (C) 2010 LinBox
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
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/** @file matrix/random-matrix.h
 * @ingroup matrix
 * @brief Implementation of random matrices.
 *
 * We provide function to create random matrices (dense, sparse, structured)
 * on several rings. This header was first introduced to avoid code redundancy in tests/
 * and make it easier to write tests/ examples/.
 *
 * @todo à la vector/stream.h
 */

#include "linbox/matrix/blas-matrix.h"
#include "linbox/randiter/random-integer.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/matrix/matrix-permutation.h"
#include "linbox/algorithms/blas-domain.h"

#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/cra-full-multip-fixed.h"

namespace LinBox
{

	struct zarb {
		zarb() {}
	} ;
	struct	_LU_ : public zarb {
		_LU_(){}
	} ;
	// LU_sparse_,
	// LU_cra_,
	struct _Rank_update_ : public zarb {
		_Rank_update_() {}
	} ;

	/// random method for constructing rank
	struct RankBuilder {
		// private :
		// balancedLC_
		// public:
		typedef _LU_                   LU_ ;
		typedef _Rank_update_ Rank_update_ ;
		RankBuilder(){}
	};

	/// Random Dense Matrix builder.
	template<class Randiter, class Field>
	class RandomDenseMatrix {

	protected:

	private :
		Field      F_ ; //!< The field containing the random entries. @todo is there a copy made ?
		/*! How are entries generated ?
		 * @pre need only provide <code>elmt& random(elmt&);</code>
		 * @see \ref LinBox::RandIterArchetype
		 */
		Randiter   R_ ;
		// Matrix    & A_ ; //!< The resulting random matrix
	public :
		/// constructor
		RandomDenseMatrix(Field & F, Randiter & R) :
			F_(F), R_(R)
		{  }

		RandomDenseMatrix(const Field & F, Randiter & R) :
			F_(F), R_(R)
		{  }

		/// destructor
		~RandomDenseMatrix() {}

		/*! creates a randomly filled matrix.
		 * @param A matrix to be randomized.
		 */
		template<class Matrix>
		Matrix & random(Matrix & A) ;

		/*! provide a matrix with prescribed rank.
		 * @param rank expected rank
		 * @param meth how is the matrix generated ? see \ref RankBuilder.
		 * @warning No certificate yet.
		 */
		template<class Matrix,class Method>
		Matrix & randomRank(Matrix & A, int rank
				    , const Method & meth );

		template<class Matrix>
		Matrix & randomRank(Matrix & A, int rank
				    , const RankBuilder::LU_ & meth );


		template<class Matrix>
		Matrix & randomRank(Matrix & A, int rank
				    , const RankBuilder::Rank_update_ & meth );


		/*! provide a matrix with prescribed rank.
		 * Default method.
		 * @param rank expected rank
		 * @warning No certificate yet.
		 */
		template<class Matrix>
		Matrix & randomRank(Matrix & A, int rank);


		// template<class Matrix>
		// void randomInvertible();

		// void randomNilpotent(int nil k); // P N_k P^(-1)
		//
		// Matrix& getMatrix(Matrix & A) {
		// return A = A_ ;
		// }

	};

	template<class Randiter, class Field>
	template<class Matrix>
	Matrix &
	RandomDenseMatrix<Randiter, Field>::random( Matrix & A)
	{
		// A_ = A ; // no copy ?
		for (size_t i = 0 ; i < A.rowdim() ; ++i)
			for (size_t j = 0 ; j < A.coldim() ; ++j)
				R_.random(A.refEntry(i,j));
		return A;
	}

	template<class Randiter, class Field>
	template<class Matrix,class Method>
	Matrix &
	RandomDenseMatrix<Randiter, Field>::randomRank(Matrix & A,
						       int    rank
						       , const Method & meth )
	{
		throw NotImplementedYet(__func__,__FILE__,__LINE__);
	}


	template<class Randiter, class Field>
	template<class Matrix>
	Matrix &
	RandomDenseMatrix<Randiter, Field>::randomRank(
						       Matrix & A
						       , int  rank
						      )
	{
		return randomRank(A,rank,RankBuilder::LU_());
	}

#if 0 /*  invalid use of incomplete type...  */
	template<class Randiter>
	DenseMatrix<PID_integer> & RandomDenseMatrix<Randiter,PID_integer>::randomRank(
												 DenseMatrix<PID_integer> & A
												 , int rank
												)
	{
		return randomRank(A,rank,RankBuilder::Rank_update_());
	}

#endif

	template<class Randiter, class Field>
	template< class Matrix>
	Matrix &
	RandomDenseMatrix<Randiter, Field>::randomRank(
						       Matrix  & A
						       , int   rank
						       , const RankBuilder::LU_ & meth
						      )
	{
		return random_lu_rank(F_,R_,A,rank);
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
		return random_rankupdate(F_,R_,A,rank);
	}

	/*  dense matrix with random entries */

	/// @todo To be factorized.
	void RandomBlasPermutation(BlasPermutation<size_t> & P)
	{
		size_t * Pt = P.getWritePointer();
		// size_t n = P.getSize();
		size_t r = P.getOrder();
		size_t n = r ; // no size given ?
		for (size_t i = 0 ; i < r ; ++i) {
			Pt[i] = i + (n-i)*( drand48() ) ;
		}
		return ;
	}

	template<class Randiter,class Field>
	DenseMatrix<Field> &
	random_lu_rank(const Field        & F,
		       const Randiter     & R,
		       DenseMatrix<Field> & A,
		       int                & rank)
	{
		size_t m = A.rowdim() ;
		size_t n = A.coldim() ;

		linbox_check(m != 0);
		linbox_check(n != 0);
		if (rank == -1) rank = std::min(m,n)*drand48() ;
		linbox_check(!(rank<0 || rank>std::min(m,n)));

		typedef typename Field::Element Element ;

		// be ready for PLUQ
		//size_t  * P = new size_t [m]   ;
		BlasPermutation<size_t> P(m);
		//Element * L = new Element[m*m] ;
		//TriangularBlasMatrix<Element> L(m,m,BlasTag::low,BlasTag::unit);
		//! @todo !!!
		BlasMatrix<Element> L(m,m);
		// Element * U = new Element[m*n] ;
		// TriangularBlasMatrix<Element> U(m,n,BlasTag::up,BlasTag::nonunit);
		//
		BlasMatrix<Element> U(m,n);
		//size_t  * Q = new size_t [n]   ;
		BlasPermutation<size_t> Q(n);

		// be ready for random elements
		NonzeroRandIter<Randiter, Field> Rnz(F,R);
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
		 * @todo create BMD.applyP(A,P,BlasTag::Left) ;
		 * avec P : BlasPermutation
		 * ou P : MatrixPermutation
		 * @todo BlasPermutation a un ordre \p p et une taille \p r distinctes !!!
		 */

		BMD.mulin_left(A_ptr,Q);
		BMD.mulin_right(P,A_ptr);

		return A ;

	}

#if 0 /*  bad LU mod p + CRA idea */
	template<class Field>
	struct LU_mod_p {
		typedef typename Field::Element Element;
		long int r_  ;
		size_t   m_ ;
		size_t   n_ ;
		LU_mod_p (size_t rows, size_t cols, long int rank) :
			r_(rank),m_(rows),n_(cols)
		{} ;

		BlasMatrix<Element>& operator()(BlasMatrix<Element>& A, const Field& F) const
		{
			BlasMatrix<Element> AA(m_,n_);
			A = AA ; // grrrrrr....
			// std::cout << m_ << ',' << A.rowdim();

			// linbox_check( m_ == A.rowdim() );
			// linbox_check( n_ == A.coldim() ) ;

			linbox_check(!(r_<0 || r_>(long int)std::min(m_,n_)));


			BlasMatrix<Element> L(m_,m_);
			BlasMatrix<Element> U(m_,n_);

			// be ready for random elements
			typedef typename Field::RandIter Randiter  ;
			Randiter R(F) ;
			NonzeroRandIter<Field> Rnz(F,R);
			Element one,zero;
			F.init(one,1UL);
			F.init(zero,0UL);

			/* Create L a random invertible lower unit triangular matrix (m x m format) */
			for (size_t j=0 ; j<m_ ; ++j)
				for (size_t i=j+1; i<m_;++i)
					R.random( L.refEntry( i,j ) );
#if 1
			for (size_t i = 0; i < (size_t) m_; ++i)
				Rnz.random( L.refEntry( i,i ) ); // non zero diagonal
			for (size_t j=0 ; j<m_ ; ++j)
				for (size_t i=0; i<j;++i)
					L.setEntry( i,j,zero );


#endif

			/* Create U a random rank-r upper non-unit triangular matrix (m x n format) */
			for (size_t i = 0 ; i < (size_t) r_; ++i)
				for (size_t j = i+1; j < n_; ++j)
					R.random( U.refEntry( i,j ) ); // r random 'triangular' lines
			for (size_t i = 0; i < (size_t) r_; ++i)
				Rnz.random( U.refEntry( i,i ) ); // non zero diagonal on rank first lines
#if 1
			for (size_t i = r_ ; i < m_ ; ++i)
				for (size_t j = i ; j < n_ ; ++j)
					U.setEntry( i,j,zero ) ; //  zero on remaining 'triangular' lines
#endif

			BlasMatrixDomain< Field > BMD(F) ;
			// LU
			/**
			 * @todo : L = [[L1,0],[A,L2]] ;U = [[U1,B],[0,U2]] ; LU = [[ rec(L1,U1), ftrmm(L1,B)],[ftrmm(A,U1),fgemm(A,B)+rec(L2,U2) ]]
			 * de même UL
			 */

			// BlasMatrix<Element> A_ptr(A) ;

			// std::cout << "A" << F.characteristic() << " := " ;
			// A.write(std::cout,true) << ';' << std::endl;

			// std::cout << "L" << F.characteristic() << " := " ;
			// L.write(std::cout,true) << ';' << std::endl;

			// std::cout << "U" << F.characteristic() << " := " ;
			// U.write(std::cout,true) << ';' << std::endl;
			// exit(-1);
			BMD.mul(A,L,U);


			return A ;

		}

		BlasMatrix<Element> initialize()
		{
			return BlasMatrix<Element>(m_,n_);
		}

	};

	template<class Field> struct CRATemporaryVectorTrait< LU_mod_p<Field>, typename Field::Element>
	{
		//typedef BlasMatrix<double>::pointer Type_t;
		typedef BlasMatrix<typename Field::Element> Type_t;
		// typedef FixedBlasMatrixDouble Type_t;
	};


#define _LB_LOG2 0.6931471807
	template<class Randiter>
	BlasMatrix<integer> &
	random_lu_rank(const PID_integer & F,
					     const Randiter & R,
					     BlasMatrix<integer> &A,
					     int & rank)
	{

		typedef Modular<double> Field ;
		typedef Field::Element Element ;
		/*- CRA -*/
		size_t PrimeSize  = 10;
		double UpperBound = _LB_LOG2 * R.getBits() ;
		RandomPrimeIterator genprime( PrimeSize );
		// ChineseRemainder< FullMultipFixedCRA< Field > > cra( std::pair<size_t,double>(A.rowdim()*A.coldim(), UpperBound) );
		ChineseRemainder< FullMultipBlasMatCRA< Field > > cra( std::pair<size_t,double>(A.rowdim()*A.coldim(), UpperBound) );

		LU_mod_p<Field> Iteration(A.rowdim(),A.coldim(),rank);

		// typename BlasMatrix<integer>::RawIterator A_it = A.rawBegin() ;

		cra(A, Iteration, genprime);

		/*- P A Q -*/
		BlasPermutation<size_t>  P(A.rowdim());
		BlasPermutation<size_t>  Q(A.coldim());

		RandomBlasPermutation(P);
		RandomBlasPermutation(Q);

		PID_integer ZZ ;
		BlasMatrixDomain< PID_integer > BMD(ZZ) ;

		// BMD.mulin_left (A,Q);
		// BMD.mulin_right(P,A);
		// throw(NotImplementedYet(__func__,__FILE__,__LINE__));
		return A ;
	}
#undef _LB_LOG2
#endif

	template<class Randiter>
	BlasMatrix<integer> &
	random_lu_rank(const PID_integer   & ZZ,
		       const Randiter      & R,
		       BlasMatrix<integer> &A,
		       int                 & rank)
	{

		size_t m = A.rowdim() ;
		size_t n = A.coldim() ;

		linbox_check(m != 0);
		linbox_check(n != 0);
		if (rank == -1) rank = std::min(m,n)*drand48() ;
		linbox_check(!(rank<0 || rank>(int)std::min(m,n)));


		// be ready for PLUQ
		BlasPermutation<size_t> P(m);
		BlasMatrix<integer> L(m,m);
		BlasMatrix<integer> U(m,n);
		BlasPermutation<size_t> Q(n);

		// be ready for random elements
		Randiter S_(R);
		S_.setBits(R.getBits()-1);
		RandomIntegerIter<false> T_(3);
		NonzeroRandIter<PID_integer,RandomIntegerIter<false> > U_(ZZ,T_);

		integer one(1),zero(0);

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

		BlasMatrixDomain< PID_integer > BMD(ZZ) ;
		MatrixDomain< PID_integer > MD(ZZ) ;
		// LU
		// L.write(std::cout << "L:=",true ) << ';' << std::endl;
		// L.write(std::cout << "U:=",true ) << ';' << std::endl;

		MD.mul(A,L,U);
		// A.write(std::cout << "pre A=",true) << std::endl;

		BMD.mulin_left(A,Q);
		BMD.mulin_right(P,A);

		// P.write(std::cout<<"P:=",false) << std::endl;
		// P.write(std::cout<<"Q:=",false) << std::endl;
		// P.write(std::cout<<"Q:=",true) << std::endl;

		return A ;
	}

#if 0 /*  DenseMatrixBase spec. */
	template<class Randiter, class Field>
	DenseMatrix<Field> &
	random_rankupdate( const Field        & F
			   ,const Randiter    & R
			   ,DenseMatrix<Field>& A
			   , int              & rank
			 )
	{
		size_t m = A.rowdim();
		size_t n = A.coldim();


		DenseMatrix<Field> D(F,m,rank) ;
		DenseMatrix<Field> G(F,rank,n) ;
		RandomDenseMatrix<Randiter, Field> RandMatGen(F,R);
		RandMatGen.random(D) ;
		RandMatGen.random(G) ;
		MatrixDomain<Field> MD(F);
		MD.mul(A,D,G);
		return A ;
	}

	template<class Randiter>
	DenseMatrix<PID_integer> &
	random_rankupdate( PID_integer               & F //!@bug const !
			   ,const Randiter           & R
			   ,DenseMatrix<PID_integer> & A
			   , int                     & rank
			 )
	{
		//! @todo check randomness
		size_t m = A.rowdim();
		size_t n = A.coldim();

		DenseMatrix<PID_integer> D(F,(size_t)m,(size_t)rank) ;
		DenseMatrix<PID_integer> G(F,(size_t)rank,(size_t)n) ;
		Randiter S_(R);
		S_.setBits(R.getBits()-1);
		RandomDenseMatrix<Randiter,PID_integer > RandMatGen(F,S_);
		RandMatGen.random(D) ;
		RandomIntegerIter<false> T_(3);
		RandomDenseMatrix<RandomIntegerIter<false>,PID_integer > RandSmallMatGen(F,T_);
		RandMatGen.random(G) ;
		MatrixDomain<PID_integer> MD(F);
		MD.mul(A,D,G);
		return A ;
	}
#endif

	template<class Randiter>
	BlasMatrix<integer> &
	random_rankupdate( PID_integer          & F //!@bug const !
			   ,const Randiter      & R
			   ,BlasMatrix<integer> & A
			   , int                & rank
			 )
	{
		size_t m = A.rowdim();
		size_t n = A.coldim();

		BlasMatrix<integer> D((size_t)m,(size_t)rank) ;
		BlasMatrix<integer> G((size_t)rank,(size_t)n) ;
		Randiter S_(R);
		S_.setBits(R.getBits()-1);
		RandomDenseMatrix<Randiter,PID_integer > RandMatGen(F,S_);
		RandMatGen.random(D) ;
		RandomIntegerIter<false> T_(3);
		RandomDenseMatrix<RandomIntegerIter<false>,PID_integer > RandSmallMatGen(F,T_);
		RandMatGen.random(G) ;

		MatrixDomain<PID_integer> MD(F);
		MD.mul(A,D,G);

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

}


/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

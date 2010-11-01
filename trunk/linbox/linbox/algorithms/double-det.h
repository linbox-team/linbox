/* linbox/algorithms/doubledet.h
 * Copyright (C) LinBox
 * 
 *  Written by Clement Pernet <clement.pernet@gmail.com>
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_doubledet_H
#define __LINBOX_doubledet_H

#include "linbox/ffpack/ffpack.h"
#include "linbox/algorithms/matrix-hom.h"	
#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/cra-full-multip.h"
#include "linbox/algorithms/cra-early-multip.h"
#include "linbox/randiter/random-prime.h"	
#include "linbox/solutions/solve.h"
#include "linbox/solutions/methods.h"
#include <vector>
namespace LinBox 
{
	
	/* Given a (n-1) x n full rank matrix A and 2 vectors b,c mod p,
	 * compute
	 * d1 = det( [ A ] ) mod p and d2 = det ( [ A ]) mod p
	 *           [ b ]                        [ c ]
	 *
	 * 
	 * A, b, c will be overwritten.
	 */
	template <class Field>
	void doubleDetModp (const Field& F, const size_t N,
			    typename Field::Element& d1,
			    typename Field::Element& d2,
			    typename Field::Element* A, const size_t lda,
			    typename Field::Element* b, const size_t incb,
			    typename Field::Element* c, const size_t incc){

		size_t* P = new size_t[N];
		size_t* Qt = new size_t[N-1];
		
		FFPACK::LUdivine (F, FFLAS::FflasUnit, FFLAS::FflasNoTrans, N-1, N,
				  A, lda, P, Qt);
		typename Field::Element d;

		// Multiplying all (N-1) first pivots)
		F.init(d, 1UL);
		for (size_t i=0; i<N-1; ++i)
			F.mulin (d, *(A + i*(lda+1)));

		bool count = false;
		for (size_t i=0;i<N-1;++i)
			if (P[i] != i) count = !count;
		if (count)
			F.negin(d);

		// Trick: instead of Right-Trans, do Left-NoTrans in order to use inc*
		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, 1, 0, N-1,
				b, incb, P);
		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, 1, 0, N-1,
				c, incc, P);
		FFLAS::ftrsv (F, FFLAS::FflasUpper, FFLAS::FflasTrans,
			      FFLAS::FflasUnit, N, A, lda, b, incb);
		FFLAS::ftrsv (F, FFLAS::FflasUpper, FFLAS::FflasTrans,
			      FFLAS::FflasUnit, N, A, lda, c, incc);


		F.mul (d1, d, *(b + (N-1) * incb));
		F.mul (d2, d, *(c + (N-1) * incc));
		
		delete[] P;
		delete[] Qt;
	}

	// Iteration class for the Chinese remaindering phase
	template <class BlackBox>
	struct IntegerDoubleDetIteration{

		const BlackBox& _A;
		const typename BlackBox::Field::Element _s1;
		const typename BlackBox::Field::Element _s2;
		IntegerDoubleDetIteration(const BlackBox& A,
					  const typename BlackBox::Field::Element& s1,
					  const typename BlackBox::Field::Element& s2)
			: _A(A), _s1(s1), _s2(s2) {}

		template<class Field>
		std::vector<typename Field::Element>&
		operator () (std::vector<typename Field::Element>& dd,
			     const Field& F) const {

			typedef typename BlackBox::template rebind<Field>::other FBlackbox;
			typename Field::Element s1p, s2p;
			dd.resize(2);
			F.init(s1p, _s1);
			F.init(s2p, _s2);
			FBlackbox Ap(_A, F);
			const size_t N = _A.coldim();
			//Timer tim;
			//tim.clear();
			//tim.start();
			doubleDetModp (F,  N, dd[0], dd[1],
				       Ap.getPointer(), Ap.getStride(),
				       Ap.getPointer() + (N-1) * Ap.getStride(), 1,
				       Ap.getPointer() + N * Ap.getStride(), 1);

			F.divin (dd[0], s1p);
			F.divin (dd[1], s2p);
			//tim.stop();
			//std::cerr<<"doubleDetModp took "<<tim.usertime()<<std::endl;
			return dd;
		}
	};

	/* Computes the actual Hadamard bound of the matrix A by taking the minimum 
	 * of the column-wise and the row-wise euclidean norm.
	 * 
	 * A is supposed to be over integers.
	 *
	 */
	// should use multiprec floating pt arith, wait, maybe not!
	template<class BlackBox>
	integer& HadamardBound (integer& hadamarBound, const BlackBox& A){
		
		
		integer res1 = 1;
		integer res2 = 1;
		integer temp;
		
		typename BlackBox::ConstRowIterator rowIt;
		typename BlackBox::ConstRow::const_iterator col;
		
		for( rowIt = A.rowBegin(); rowIt != A.rowEnd(); ++rowIt ) {
			temp = 0;
			for( col = rowIt->begin(); col != rowIt->end(); ++col )
				temp += static_cast<integer>((*col)) * (*col);
			res1 *= temp;
		}
		res1 = sqrt(res1);

		typename BlackBox::ConstColIterator colIt;
		typename BlackBox::ConstCol::const_iterator row;
		
		for( colIt = A.colBegin(); colIt != A.colEnd(); ++colIt ) {
			temp = 0;
			for( row = colIt->begin(); row != colIt->end(); ++row )
				temp += static_cast<integer>((*row)) * (*row);
			res2 *= temp;
		}
		res2 = sqrt(res2);
		
		return hadamarBound = MIN(res1, res2);
	}
// 	template<class BlackBox>
// 	double& HadamardBound (double& hadamarBound, const BlackBox& A){
		
		
// 		double res1 = 0;
// 		double res2 = 0;
// 		integer temp;
		
// 		typename BlackBox::ConstRowIterator rowIt;
// 		typename BlackBox::ConstRow::const_iterator col;
		
// 		for( rowIt = A.rowBegin(); rowIt != A.rowEnd(); ++rowIt ) {
// 			temp = 0;
// 			for( col = rowIt->begin(); col != rowIt->end(); ++col )
// 				temp += static_cast<integer>((*col)) * (*col);
// 			res1 += log ((double)temp);
// 		}
// 		res1 = res1/2;

// 		typename BlackBox::ConstColIterator colIt;
// 		typename BlackBox::ConstCol::const_iterator row;
		
// 		for( colIt = A.colBegin(); colIt != A.colEnd(); ++colIt ) {
// 			temp = 0;
// 			for( row = colIt->begin(); row != colIt->end(); ++row )
// 				temp += static_cast<integer>((*row)) * (*row);
// 			res2 += log((double)temp);
// 		}
// 		res2 = res2/2;
		
// 		return hadamarBound = MIN(res1, res2);
// 	}

	/* Given a (N+1) x N full rank matrix
	 * [ A ]
	 * [ b ]
	 * [ c ]
	 * where b and c are the last 2 rows,
	 * and given s1 and s2, respectively divisors of d1 and d2 defined as
	 * d1 = det ([ A ])
	 *          ([ b ])
	 * and
	 * d2 = det ([ A ])
	 *          ([ c ]),
	 * compute d1 and d2.
	 * Assumes d1 and d2 are non zero.
	 * Result is probablistic if proof=true
	 */
	template <class BlackBox>
	void doubleDetGivenDivisors (const BlackBox& A,
				     typename BlackBox::Field::Element& d1,
				     typename BlackBox::Field::Element& d2,
				     const typename BlackBox::Field::Element& s1,
				     const typename BlackBox::Field::Element& s2,
				     const bool proof){

		typename BlackBox::Field F = A.field();
		IntegerDoubleDetIteration<BlackBox> iteration(A, s1, s2);
		// 0.7213475205 is an upper approximation of 1/(2log(2))
		RandomPrimeIterator genprime( 25-(int)ceil(log((double)A.rowdim())*0.7213475205)); 
		
		std::vector<typename BlackBox::Field::Element> dd;
		if (proof) {
			integer bound;
			double logbound;
			//Timer t_hd,t_cra;
			//t_hd.clear();
			//t_hd.start();
			HadamardBound (bound, A);
			logbound = (logtwo (bound) - logtwo (MIN(abs(s1),abs(s2))))*0.693147180559945;
			//t_hd.stop();
			//std::cerr<<"Hadamard bound = : "<<logbound<<" in "<<t_hd.usertime()<<"s"<<std::endl;

			ChineseRemainder <FullMultipCRA <Modular <double> > > cra(logbound);

			//t_hd.clear();
			//t_cra.start();
			cra (dd, iteration, genprime);
			//t_cra.stop();
			//std::cerr<<"CRA : "<<t_cra.usertime()<<"s"<<std::endl;

		} else {
			ChineseRemainder <EarlyMultipCRA <Modular<double> > > cra(4UL);
			cra (dd, iteration, genprime);
		}
		F.mul (d1, dd[0], s1);
		F.mul (d2, dd[1], s2);
	}

	/* Given a (n + 1) x n full rank matrix
	 * A = [ B ]
	 *     [ b ]
	 *     [ c ]
	 * over Z, compute
	 * d1 = det( [ A ] ) and d2 = det ( [ A ])
	 *           [ b ]                  [ c ]
	 */
	
	template <class BlackBox>
	void doubleDet (typename BlackBox::Field::Element& d1,
			typename BlackBox::Field::Element& d2,
			const BlackBox& A,
			//const vector<typename BlackBox::Field::Element>& b,
			//const vector<typename BlackBox::Field::Element>& c,
			bool proof){

		linbox_check (A.coldim() == A.rowdim()+1);
		
		const size_t N = A.coldim();
		//		BlasBlackbox<typename BlackBox::Field> B (A,0,0,N,N);
		BlasBlackbox<typename BlackBox::Field> B (A.field(),N,N);
		typename BlackBox::Field::Element den1, den2;
		std::vector<typename BlackBox::Field::Element> x1(N);
		for (size_t i=0; i<N; ++i){
			x1[i] = 0;
			for (size_t j=0; j<N; ++j)
				B.setEntry(j,i, A.getEntry(i,j));
		}
		// 		for (size_t i=0; i<N; ++i)
		// 			B.setEntry (i, N-1, b[i]);

		std::vector<typename BlackBox::Field::Element> c(N);
		for (size_t i=0; i<N; ++i)
			c[i] = A.getEntry (N, i);
		//Timer tim;
		//tim.clear();
		//tim.start();
		solve (x1, den1, B, c);
		//tim.stop();
		//std::cerr<<"Solve took "<<tim.usertime()<<std::endl;
		
		den1 = den1;
		// Should work:
		// den (y[n]) = den (-den1/x[n]) = x[n] 
		den2 = -x1[N-1];
		//tim.clear();
		//tim.start();
		doubleDetGivenDivisors (A, d1, d2, den1, den2, proof);
		//tim.stop();
		//std::cerr<<"CRA "<<(proof?"Determinist":"Probablistic")
		//	 <<" took "<<tim.usertime()<<std::endl;
	}
}


#endif // __LINBOX_doubledet_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

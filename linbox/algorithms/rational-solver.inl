/* linbox/algorithms/rational-solver.inl
 * Copyright (C) 2004 Pascal Giorgi
 *
 * Written by Pascal Giorgi  <pascal.giorgi@ens-lyon.fr>
 * Modified by David Pritchard  <daveagp@mit.edu>
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

#ifndef __LINBOX_rational_solver_INL
#define __LINBOX_rational_solver_INL

#include "linbox/util/debug.h"
#include "linbox/linbox-config.h"

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/blackbox/lambda-sparse.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"
#include "linbox/algorithms/lifting-container.h"
#include "linbox/algorithms/rational-reconstruction.h"
#include "linbox/algorithms/matrix-inverse.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/algorithms/gauss.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/massey-domain.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/block-massey-domain.h"
#include "linbox/algorithms/vector-fraction.h"
#include <fflas-ffpack/ffpack/ffpack.h>
#include <fflas-ffpack/fflas/fflas.h>
#include "linbox/solutions/methods.h"
#include "linbox/blackbox/block-hankel-inverse.h"

#include "linbox/vector/blas-vector.h"

// #ifdef __LINBOX_BLAS_AVAILABLE
#include "linbox/config-blas.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/factorized-matrix.h"
#include "linbox/util/timer.h"
// #endif

//#define DEBUG_DIXON
//#define DEBUG_INC
//#define SKIP_NONSINGULAR

namespace LinBox
{
	// SPECIALIZATION FOR WIEDEMANN

	// note: if Vector1 != Vector2 compilation of solve or solveSingluar will fail (via an invalid call to sparseprecondition)!
	// maybe they should not be templated separately, or sparseprecondition should be rewritten

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>
	SolverReturnStatus
	DixonSolver<Ring,Field,RandomPrime,Method::Wiedemann>::solve (Vector1& num,
								       Integer& den,
								       const IMatrix& A,
								       const Vector2& b,
								       const bool old,
								       int maxPrimes) const
	{
		SolverReturnStatus status=SS_FAILED;

		switch (A.rowdim() == A.coldim() ? solveNonsingular(num, den, A,b) : SS_SINGULAR) {

		case SS_OK:
			status=SS_OK;
			break;

		case SS_SINGULAR:
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING) << "switching to singular" << std::endl;
			//std::cerr<<"switching to singular\n";
			status=solveSingular(num, den,A,b);
			break;

		case SS_FAILED:
			break;

		default:
			throw LinboxError ("Bad return value from solveNonsingular");

		}

		return status;
	}

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>
	SolverReturnStatus
	DixonSolver<Ring, Field, RandomPrime, Method::Wiedemann>::solveNonsingular( Vector1& num,
												      Integer& den,
												      const IMatrix& A,
												      const Vector2& b,
												      int maxPrimes) const
	{
        commentator().start("solve.dixon.integer.nonsingular.wiedemann");
		// checking if matrix is square
		linbox_check(A.rowdim() == A.coldim());

		// checking size of system
		linbox_check(A.rowdim() == b.size());



		SparseMatrix<Field> *Ap;
		FPolynomial MinPoly;
		size_t  deg;
		size_t issingular = SINGULARITY_THRESHOLD;
		static Field *F=NULL;
		Prime prime = _prime;
		do {
#ifdef RSTIMING
			tNonsingularSetup.clear();
			tNonsingularSetup.start();
#endif
			_prime = prime;
			if (F != NULL) delete F;
			F=new Field(prime);
			Ap = new SparseMatrix<Field>(A, *F);
			typename Field::RandIter random(*F);
			BlackboxContainer<Field,SparseMatrix<Field> > Sequence(Ap,*F,random);
			MasseyDomain<Field,BlackboxContainer<Field,SparseMatrix<Field> > > MD(&Sequence);
#ifdef RSTIMING
			tNonsingularSetup.stop();
			ttNonsingularSetup+=tNonsingularSetup;
			tNonsingularMinPoly.clear();
			tNonsingularMinPoly.start();
#endif
			MD.minpoly(MinPoly,deg);
#ifdef RSTIMING
			tNonsingularMinPoly.stop();
			ttNonsingularMinPoly+=tNonsingularMinPoly;
#endif
			prime = *_genprime;
		}
		while(F->isZero(MinPoly.front()) && --issingular );



		if (!issingular){
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING) << "The Matrix is singular" << std::endl;
//			std::cerr<<"The Matrix is singular\n";
			delete Ap;
			return SS_SINGULAR;
		}
		else {

			typedef SparseMatrix<Field> FMatrix;

			typedef WiedemannLiftingContainer<Ring, Field, IMatrix, FMatrix, FPolynomial> LiftingContainer;

			LiftingContainer lc(_ring, *F, A, *Ap, MinPoly, b,_prime);

			RationalReconstruction<LiftingContainer> re(lc);

			re.getRational(num, den, 0);
#ifdef RSTIMING
			ttNonsingularSolve.update(re, lc);
#endif
            commentator().stop("solve.dixon.integer.nonsingular.wiedemann");
			return SS_OK;
		}
	}

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>
	SolverReturnStatus
	DixonSolver<Ring,Field,RandomPrime,Method::Wiedemann>:: solveSingular (Vector1& num,
										     Integer& den,
										     const IMatrix& A,
										     const Vector2& b,
										     int maxPrimes) const
	{
        commentator().start("solve.dixon.integer.singular.wiedemann");

		typedef BlasVector<Ring>  IVector;
		typedef SparseMatrix<Field>                  FMatrix;

		// checking size of system
		linbox_check(A.rowdim() == b.size());

		typedef LambdaSparseMatrix<Ring>  IPreconditioner;
		typedef LambdaSparseMatrix<Field> FPreconditioner;

		typedef Compose<IPreconditioner,Compose<IMatrix,IPreconditioner> > IPrecondMatrix;
		typedef Compose<FPreconditioner,Compose<FMatrix,FPreconditioner> > FPrecondMatrix;

		FMatrix               *Ap;
		IPreconditioner *P     =NULL;
		IPreconditioner *Q     =NULL;
		FPreconditioner *Pmodp =NULL;
		FPreconditioner *Qmodp =NULL;
		IPrecondMatrix  *PAQ   =NULL;
		FPrecondMatrix  *PApQ  =NULL;
		IVector Pb;


		FPolynomial MinPoly;
		size_t  deg;
		size_t badprecondition = BAD_PRECONTITIONER_THRESHOLD;
		Field *F;
		Prime prime = _prime;
		typename Field::Element tmp;
		do {
			if (PApQ != NULL) {
				delete P;
				delete Q;
				delete PApQ;
				delete PAQ;
			}
			_prime = prime;
			F=new Field(prime);//std::cerr<<"here\n";
			Ap = new FMatrix(A, *F);
			sparseprecondition (*F,&A,PAQ,Ap,PApQ,b,Pb,P,Q,Pmodp,Qmodp);
			typename Field::RandIter random(*F);
			BlackboxContainer<Field,FPrecondMatrix> Sequence(PApQ,*F,random);
			MasseyDomain<Field,BlackboxContainer<Field,FPrecondMatrix> > MD(&Sequence);

			MD.minpoly(MinPoly,deg);
			//MinPoly.resize(3);MinPoly[0]=1;MinPoly[1]=2;MinPoly[2]=1;
			prime = _genprime.randomPrime();
			F->add(tmp,MinPoly.at(1),MinPoly.front());
		}
		while(((F->isZero(tmp) || MinPoly.size() <=2) && --badprecondition ));

#ifdef DEBUG_DS_WIEDEMANN
		std::clog<<"minpoly found with size: "<<MinPoly.size()<<std::endl;
		for (size_t i=0;i<MinPoly.size();++i)
			std::clog<<MinPoly[i]<<"*x^"<<i<<"+";
		std::clog<<std::endl;
		std::clog<<"prime is: "<<_prime<<std::endl;
#endif

		if (!badprecondition){
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING) << "Bad Preconditionner" << std::endl;
//			std::cerr<<"Bad Preconditionner\n";

			delete Ap;
			if (PAQ  != NULL) delete PAQ;
			if (PApQ != NULL) delete PApQ;
			if (P    != NULL) delete P;
			if (Q    != NULL) delete Q;

			return SS_BAD_PRECONDITIONER;
		}
		else {

			MinPoly.erase(MinPoly.begin());

			typedef WiedemannLiftingContainer<Ring, Field, IPrecondMatrix, FPrecondMatrix, FPolynomial> LiftingContainer;

#ifdef DEBUG_DS_WIEDEMANN
			std::clog<<"before lc\n";
#endif

			LiftingContainer lc(_ring, *F, *PAQ, *PApQ, MinPoly, Pb, _prime);

#ifdef DEBUG_DS_WIEDEMANN
			std::clog<<"constructing lifting container of length: "<<lc.length()<<std::endl;
#endif

			RationalReconstruction<LiftingContainer> re(lc,_ring,2);

			re.getRational(num, den, 0);


			if (Q    != NULL) {
				IVector Qx(num.size());
				Q->apply(Qx, num);
				num = Qx;
			}


			delete Ap;
			if (PAQ  != NULL) delete PAQ;
			if (PApQ != NULL) delete PApQ;
			if (P    != NULL) delete P;
			if (Q    != NULL) delete Q;

            commentator().stop("solve.dixon.integer.singular.wiedemann");
			return SS_OK;
		}
	}


	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class FMatrix, class IVector>
	void
	DixonSolver<Ring,Field,RandomPrime,Method::Wiedemann>::sparseprecondition (const Field& F,
										     const IMatrix *A,
										     Compose<LambdaSparseMatrix<Ring>, Compose<IMatrix, LambdaSparseMatrix<Ring> > > *&PAQ,
										     const FMatrix *Ap,
										     Compose<LambdaSparseMatrix<Field>, Compose<FMatrix, LambdaSparseMatrix<Field> > > *&PApQ,
										     const IVector& b,
										     IVector& Pb,
										     LambdaSparseMatrix<Ring> *&P,
										     LambdaSparseMatrix<Ring> *&Q,
										     LambdaSparseMatrix<Field> *&Pmodp,
										     LambdaSparseMatrix<Field> *&Qmodp) const
	{
		commentator().start ("Constructing sparse preconditioner");

        VectorDomain<Ring> VD(_ring);
#if DEBUG_DS_WIEDEMANN
		A->write(std::clog<<"A:\n");
		Ap->write(std::clog<<"A mod p:\n");
		VD.write(std::clog<<"b:\n",b)<<std::endl;
#endif


		typedef LambdaSparseMatrix<Ring>  IPreconditioner;
		typedef LambdaSparseMatrix<Field> FPreconditioner;

		size_t min_dim = A->coldim() < A->rowdim() ? A->coldim() : A->rowdim();

		P = new  IPreconditioner(_ring,min_dim,A->rowdim(),2,3.);
		//		std::cerr<<"P:\n";
		//		P->write(std::cerr);

		Q = new  IPreconditioner(_ring,A->coldim(),min_dim,2,3.);
		//		std::cerr<<"Q:\n";
		//		Q->write(std::cerr);

		Compose<IMatrix,IPreconditioner> *AQ;
		AQ = new Compose<IMatrix,IPreconditioner> (A,Q);

		PAQ = new Compose<IPreconditioner, Compose<IMatrix,IPreconditioner> > (P,AQ);		;
		Pb.resize(min_dim);
		P->apply(Pb,b);
		//		std::cerr<<"Pb:\n";
		//		VD.write(std::cerr,Pb)<<std::endl;

		Pmodp = new FPreconditioner(F,*P);
#if DEBUG_DS_WIEDEMANN
		std::clog<<"P mod p completed\n";
#endif

		Qmodp = new FPreconditioner(F,*Q);
#if DEBUG_DS_WIEDEMANN
		std::clog<<"Q mod p completed\n";
#endif

		Compose<FMatrix,FPreconditioner> *ApQ;
		ApQ = new Compose<FMatrix,FPreconditioner> (Ap,Qmodp);

		PApQ = new Compose<FPreconditioner, Compose<FMatrix,FPreconditioner> > (Pmodp, ApQ);

		commentator().stop ("Constructing sparse preconditioner");

	}

	// SPECIALIZATION FOR BLOCK WIEDEMANN

	// note: if Vector1 != Vector2 compilation of solve or solveSingluar will fail (via an invalid call to sparseprecondition)!
	// maybe they should not be templated separately, or sparseprecondition should be rewritten

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>
	SolverReturnStatus
	DixonSolver<Ring,Field,RandomPrime,Method::BlockWiedemann>::solve (Vector1& num, Integer& den,
										     const IMatrix& A,
										     const Vector2& b,
										     const bool old,
										     int maxPrimes) const
	{
		SolverReturnStatus status=SS_FAILED;

		switch (A.rowdim() == A.coldim() ? solveNonsingular(num, den, A,b) : SS_SINGULAR) {

		case SS_OK:
			status=SS_OK;
			break;

		case SS_SINGULAR:
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << "could switch to singular but not doing it(?)" << std::endl;
			//std::cerr<<"switching to singular\n";
			//status=solveSingular(num, den,A,b);
			break;

		case SS_FAILED:
			break;

		default:
			throw LinboxError ("Bad return value from solveNonsingular");

		}

		return status;
	}

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>
	SolverReturnStatus
	DixonSolver<Ring,Field,RandomPrime,Method::BlockWiedemann>::solveNonsingular (Vector1& num,
										       Integer& den,
										       const IMatrix& A,
										       const Vector2& b,
										       int maxPrimes) const
	{
        commentator().start("solve.dixon.integer.nonsingular.blockwiedemann");
		// checking if matrix is square
		linbox_check(A.rowdim() == A.coldim());

		// checking size of system
		linbox_check(A.rowdim() == b.size());

		size_t m,n;
		integer tmp,tmproot;
		tmp=A.coldim();
#if 0
		m = n = tmp.bitsize();
		m = n = sqrt(tmp);
		m = n = root(tmp,3); // wrong # args to root. -bds
#endif
		// m = n =
		root(tmproot, tmp,3);
		m = n = uint32_t(tmproot);
		//		std::cout<<"block factor= "<<m<<"\n";;
		typedef SparseMatrix<Field> FMatrix;

		Field F(_prime);
		FMatrix Ap(A, F);
		Transpose<FMatrix > Bp(Ap);
		//		std::cout<<"Ap:\n";
		//		Ap.write(std::cout);
		typedef BlockWiedemannLiftingContainer<Ring, Field, Transpose<IMatrix >, Transpose<FMatrix > > LiftingContainer;

		Transpose<IMatrix> B(A);

		LiftingContainer lc(_ring, F, B, Bp, b,_prime, m, n);

		RationalReconstruction<LiftingContainer> re(lc);

		re.getRational(num, den, 0);
#ifdef RSTIMING
		ttNonsingularSolve.update(re, lc);
#endif

        commentator().stop("solve.dixon.integer.nonsingular.blockwiedemann");
		return SS_OK;
	}

	// END OF SPECIALIZATION FOR BLOCK WIEDEMANN


	/*
	 * Specialization for Block Hankel method
	 */
	// solve non singular system using block Hankel
	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>
	SolverReturnStatus
	DixonSolver<Ring,Field,RandomPrime,Method::BlockHankel>::solveNonsingular(Vector1& num,
										   Integer& den,
										   const IMatrix& A,
										   const Vector2& b,
										   size_t blocksize,
										   int maxPrimes) const
	{

        commentator().start("solve.dixon.integer.nonsingular.blockhankel");
		linbox_check(A.rowdim() == A.coldim());
		linbox_check(A.rowdim() % blocksize == 0);

		// reduce the matrix mod p
		Field F(_prime);
		typedef typename IMatrix::template rebind<Field>::other FMatrix;
		FMatrix Ap(A, F);

		// precondition Ap  with a random diagonal Matrix
		typename Field::RandIter G(F,0,123456);
		BlasVector<Field> diag(F,Ap.rowdim());

		for(size_t i=0;i<Ap.rowdim();++i){
			do {
				G.random(diag[i]);
			} while(F.isZero(diag[i]));
		}



		Diagonal<Field> D(diag);

		Compose<Diagonal<Field>, FMatrix> DAp(D,Ap);

		size_t n = A.coldim();
		size_t numblock = n/blocksize;

		// generate randomly U and V
		BlasMatrix<Field> U(F,blocksize,A.rowdim()), V(A.coldim(),blocksize);

		for (size_t j=0;j<blocksize; ++j)
			for (size_t i=j*numblock;i<(j+1)*numblock;++i){
				G.random(V.refEntry(i,j));
			}
		for (size_t i=0;i<n;++i)
			G.random(U.refEntry(0,i));

#ifdef RSTIMING
		Timer chrono;
		chrono.clear();
		chrono.start();
#endif

		// compute the block krylov sequence associated to U.A^i.V
		BlackboxBlockContainerRecord<Field, Compose<Diagonal<Field>,FMatrix> >  Seq(&DAp, F, U, V, false);

#ifdef RSTIMING
		chrono.stop();
		std::cout<<"sequence generation: "<<chrono<<"\n";
		chrono.clear();
		chrono.start();
#endif

		// compute the inverse of the Hankel matrix associated with the Krylov Sequence
		BlockHankelInverse<Field> Hinv(F, Seq.getRep());
		BlasVector<Field> y(F,n), x(F,n, F.one);

#ifdef RSTIMING
		chrono.stop();
		std::cout<<"inverse block hankel: "<<chrono<<"\n";
#endif

		typedef BlockHankelLiftingContainer<Ring,Field,IMatrix,Compose<Diagonal<Field>,FMatrix>, BlasMatrix<Field> > LiftingContainer;
		LiftingContainer lc(_ring, F, A, DAp, D, Hinv, U, V, b, _prime);
		RationalReconstruction<LiftingContainer > re(lc);

		if (!re.getRational(num, den, 0)) return SS_FAILED;

#ifdef RSTIMING
		std::cout<<"lifting bound computation : "<<lc.ttSetup<<"\n";
		std::cout<<"residue computation       : "<<lc.ttRingApply<<"\n";
		std::cout<<"rational reconstruction   : "<<re.ttRecon<<"\n";
#endif

        commentator().stop("solve.dixon.integer.nonsingular.blockhankel");
		return SS_OK;
	}


	/*
	 * Specialization for Sparse Elimination method
	 */
	// solve any system using Sparse LU
	// max prime is not use. only check with one prime
	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>
	SolverReturnStatus
	DixonSolver<Ring,Field,RandomPrime,Method::SparseElimination>::solve(Vector1& num,
                                                                          Integer& den,
                                                                          const IMatrix& A,
                                                                          const Vector2& b,
                                                                          int maxPrimes) const
	{

            //linbox_check(A.rowdim() == A.coldim());

		typedef typename Field::Element Element_t;

		// reduce the matrix mod p
		const Field F(_prime);
		typedef typename IMatrix::template rebind<Field>::other FMatrix;
		FMatrix Ap(A, F);

		// compute PLUQ Factorization
		Permutation<Field> P(F,(int)A.coldim()),Q(F,(int)A.rowdim());
		FMatrix L(F, A.rowdim(), A.rowdim());
		size_t rank;
		Element_t det;

		GaussDomain<Field> GD(F);
		GD.QLUPin(rank,det,Q,L,Ap,P,Ap.rowdim(), Ap.coldim());
		if (rank < A.coldim()) {
			// Choose a nonrandom solution with smallest entries:
			// Sets solution values to 0 for coldim()-rank columns
			// Therefore, prune unnecessary elements
			// in those last columns of U
			size_t origNNZ=0,newNNZ=0;
			for(typename FMatrix::RowIterator row=Ap.rowBegin();
			    row != Ap.rowEnd(); ++row) {
				if (row->size()) {
					origNNZ += row->size();
					size_t ns=0;
					for(typename FMatrix::Row::iterator it = row->begin();
					    it != row->end(); ++it, ++ns) {
						if (it->first >= rank) {
							row->resize(ns);
							break;
						}
					}
					newNNZ += row->size();
				}
			}
			commentator().report (Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT) << "Pruned : " << (origNNZ-newNNZ) << " unnecessary elements in upper triangle" << std::endl;
		}

//         A.write(std::cout << "A:=") << ';' << std::endl;
//         Q.write(std::cout << "Q:=") << ';' << std::endl;
//         L.write(std::cout << "L:=") << ';' << std::endl;
//         Ap.write(std::cout << "Ap:=") << ';' << std::endl;
//         P.write(std::cout << "P:=") << ';' << std::endl;

		typedef SparseLULiftingContainer<Ring,Field,IMatrix,FMatrix> LiftingContainer;
		LiftingContainer lc(_ring, F, A, L, Q, Ap, P, rank, b, _prime);
		RationalReconstruction<LiftingContainer > re(lc);

		if (!re.getRational(num, den, 0))
			return SS_FAILED;
		else
			return SS_OK;
	}

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>
	SolverReturnStatus
	DixonSolver<Ring,Field,RandomPrime,Method::SparseElimination>::solveNonsingular(Vector1& num,
                                                                          Integer& den,
                                                                          const IMatrix& A,
                                                                          const Vector2& b,
                                                                          int maxPrimes) const
	{
        commentator().start("solve.dixon.integer.nonsingular.sparselim");
		return solve(num, den, A, b, maxPrimes);
        commentator().stop("solve.dixon.integer.nonsingular.sparselim");
	}

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>
	SolverReturnStatus
	DixonSolver<Ring,Field,RandomPrime,Method::SparseElimination>::solveSingular(Vector1& num,
                                                                          Integer& den,
                                                                          const IMatrix& A,
                                                                          const Vector2& b,
                                                                          int maxPrimes) const
	{
        commentator().start("solve.dixon.integer.singular.sparselim");
		return solve(num, den, A, b, maxPrimes);
        commentator().stop("solve.dixon.integer.singular.sparselim");
	}

} //end of namespace LinBox

//BB : moved the following "guarded" code in a new file, verbatim :
#include "./dixon-solver/dixon-solver-symbolic-numeric.h"

#endif //__LINBOX_rational_solver_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

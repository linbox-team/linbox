/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/algorithms/rational-solver.inl
 * Copyright (C) 2004 Pascal Giorgi
 *
 * Written by Pascal Giorgi  <pascal.giorgi@ens-lyon.fr>
 * Modified by David Pritchard  <daveagp@mit.edu>
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

#ifndef __LINBOX_RATIONAL_SOLVER_INL
#define __LINBOX_RATIONAL_SOLVER_INL

#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/blackbox/lambda-sparse.h>
#include <linbox/algorithms/lifting-container.h>
#include <linbox/algorithms/rational-reconstruction.h>
#include <linbox/algorithms/matrix-inverse.h>
#include <linbox/algorithms/matrix-mod.h>
#include <linbox/algorithms/blackbox-container.h>
#include <linbox/algorithms/massey-domain.h>
#include <linbox/algorithms/vector-fraction.h>
#include <linbox/fflapack/fflapack.h>
#include <linbox/fflas/fflas.h>
#include <linbox/solutions/methods.h>
#include <linbox/util/debug.h>
#include <linbox-config.h>
#ifdef __LINBOX_BLAS_AVAILABLE
#include <linbox/blackbox/blas-blackbox.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/algorithms/blas-domain.h>
#include <linbox/matrix/factorized-matrix.h>
#endif

//#define DEBUG_DIXON
//#define DEBUG_INC

namespace LinBox {

	// SPECIALIZATION FOR WIEDEMANN 	

	// note: if Vector1 != Vector2 compilation of solve or solveSingluar will fail (via an invalid call to sparseprecondition)!
	// maybe they should not be templated separately, or sparseprecondition should be rewritten

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,WiedemannTraits>::solve (Vector1& answer,
											  const IMatrix& A,
											  const Vector2& b,
											  const bool old=false,
											  int maxPrimes) const {

		typename RationalSolver<Ring,Field,RandomPrime,WiedemannTraits>::ReturnStatus status=SS_FAILED;

		switch (A.rowdim() == A.coldim() ? solveNonsingular(answer,A,b) : SS_SINGULAR) {

		case SS_OK:
			status=SS_OK;
			break;

		case SS_SINGULAR:
			std::cerr<<"switching to singular\n";
			status=solveSingular(answer,A,b);
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
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,WiedemannTraits>::solveNonsingular (Vector1& answer,
												     const IMatrix& A,
												     const Vector2& b,
												     int maxPrimes) const {
		// checking if matrix is square
		linbox_check(A.rowdim() == A.coldim());
		
		// checking size of system
		linbox_check(A.rowdim() == b.size());

		SparseMatrix<Field> *Ap;		
		FPolynomial MinPoly;
		unsigned long  deg;
		unsigned long issingular = SINGULARITY_THRESHOLD; 			
		Field *F=NULL;
		Prime prime = _prime;
		do {
			_prime = prime;
			if (F != NULL) delete F;
			F=new Field(prime);				
			MatrixMod::mod (Ap, A, *F);
			typename Field::RandIter random(*F);
			BlackboxContainer<Field,SparseMatrix<Field> > Sequence(Ap,*F,random);
			MasseyDomain<Field,BlackboxContainer<Field,SparseMatrix<Field> > > MD(&Sequence);
			MD.minpoly(MinPoly,deg);
			prime = _genprime.randomPrime();
		}
		while(F->isZero(MinPoly.front()) && --issingular );			
				

		if (!issingular){	
			std::cerr<<"The Matrix is singular\n";
			delete Ap;
			return SS_SINGULAR;
		}			
		else {
 			//std::cerr<<"A:\n";
			//A.write(std::cerr);
 			//std::cerr<<"A mod p:\n";
 			//Ap->write(std::cerr);
			//Ring r;
			//VectorDomain<Ring> VD(r);
			//std::cerr<<"b:\n";		
			//VD.write(std::cerr,b)<<std::endl;
			//std::cerr<<"prime: "<<_prime<<std::endl;
			
			//std::cerr<<"non singular\n";
			typedef WiedemannLiftingContainer<Ring, Field, IMatrix, SparseMatrix<Field>, FPolynomial> LiftingContainer;
			
			LiftingContainer lc(_R, *F, A, *Ap, MinPoly, b,_prime);
			
			RationalReconstruction<LiftingContainer> re(lc);
			
			re.getRational2(answer);
					
			return SS_OK;
		}
	}       

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,WiedemannTraits>::solveSingular (Vector1& answer,
												  const IMatrix& A,
												  const Vector2& b, 
												  int maxPrimes) const {
		std::cerr<<"in singular solver\n";

		typedef std::vector<typename Field::Element> FVector;
		typedef std::vector<typename Ring::Element>  IVector;
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
		unsigned long  deg;
		unsigned long badprecondition = BAD_PRECONTITIONER_THRESHOLD;
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
			MatrixMod::mod (Ap, A, *F);//std::cerr<<"after\n";
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
		std::cerr<<"minpoly found with size: "<<MinPoly.size()<<std::endl;
		for (size_t i=0;i<MinPoly.size();i++)
			std::cerr<<MinPoly[i]<<"*x^"<<i<<"+";
		std::cerr<<std::endl;

		std::cerr<<"prime is: "<<_prime<<std::endl;
		if (!badprecondition){
			std::cerr<<"Bad Preconditionner\n";
		
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
			std::cerr<<"before lc\n";
			LiftingContainer lc(_R, *F, *PAQ, *PApQ, MinPoly, Pb, _prime);
			std::cerr<<"constructing lifting container of length: "<<lc.length()<<std::endl;
			
			RationalReconstruction<LiftingContainer> re(lc,_R,2);
			
			re.getRational2(answer); 


			if (Q    != NULL) {
				typename Ring::Element lden;
				_R. init (lden, 1);
				typename Vector1::iterator p;		
				for (p = answer.begin(); p != answer.end(); ++ p)
					_R. lcm (lden, lden, p->second);


				IVector Qx(answer.size()),x(answer.size());
				typename IVector::iterator p_x;
						
				for (p = answer.begin(), p_x = x. begin(); p != answer.end(); ++ p, ++ p_x) {					
					_R. mul (*p_x, p->first, lden);					
					_R. divin (*p_x, p->second);					
				}

				Q->apply(Qx,x);
				for (p=answer.begin(),p_x=Qx.begin(); p != answer.end();++p,++p_x){
					p->first=*p_x;
					p->second=lden;
				}					
			}


			delete Ap;
			if (PAQ  != NULL) delete PAQ;
			if (PApQ != NULL) delete PApQ;	
			if (P    != NULL) delete P;
			if (Q    != NULL) delete Q;
			
			return SS_OK;
		}
	}       


	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class FMatrix, class IVector>
	void RationalSolver<Ring,Field,RandomPrime,WiedemannTraits>::sparseprecondition (const Field& F,
											 const IMatrix *A,
											 Compose<LambdaSparseMatrix<Ring>,Compose<IMatrix,LambdaSparseMatrix<Ring> > > *&PAQ,
											 const FMatrix *Ap,
											 Compose<LambdaSparseMatrix<Field>,Compose<FMatrix,LambdaSparseMatrix<Field> > > *&PApQ,
											 const IVector& b,
											 IVector& Pb,
											 LambdaSparseMatrix<Ring> *&P,
											 LambdaSparseMatrix<Ring> *&Q,
											 LambdaSparseMatrix<Field> *&Pmodp,
											 LambdaSparseMatrix<Field> *&Qmodp) const
	{
		// 		std::cerr<<"A:\n";
		// 		A->write(std::cerr);
		// 		std::cerr<<"A mod p:\n";
		// 		Ap->write(std::cerr);
		VectorDomain<Ring> VD(_R);
		// 		std::cerr<<"b:\n";		
		// 		VD.write(std::cerr,b)<<std::endl;
		

		commentator.start ("Constructing sparse preconditioner");
		typedef LambdaSparseMatrix<Ring>  IPreconditioner;
		typedef LambdaSparseMatrix<Field> FPreconditioner;

		size_t min_dim = A->coldim() < A->rowdim() ? A->coldim() : A->rowdim();
		
		P = new  IPreconditioner(_R,min_dim,A->rowdim(),2,3.);       
		// 		std::cerr<<"P:\n";
		// 		P->write(std::cerr);
		
		Q = new  IPreconditioner(_R,A->coldim(),min_dim,2,3.);
		// 		std::cerr<<"Q:\n";
		// 		Q->write(std::cerr);

		Compose<IMatrix,IPreconditioner> *AQ;
		AQ = new Compose<IMatrix,IPreconditioner> (A,Q);

		PAQ = new Compose<IPreconditioner, Compose<IMatrix,IPreconditioner> > (P,AQ);		;
		Pb.resize(min_dim);
		P->apply(Pb,b);
		// 		std::cerr<<"Pb:\n";
		// 		VD.write(std::cerr,Pb)<<std::endl;

		Pmodp = new FPreconditioner(F,*P);       
		std::cerr<<"P mod p completed\n";
		Qmodp = new FPreconditioner(F,*Q);
		std::cerr<<"Q mod p completed\n";

		Compose<FMatrix,FPreconditioner> *ApQ;
		ApQ = new Compose<FMatrix,FPreconditioner> (Ap,Qmodp);
		
		PApQ = new Compose<FPreconditioner, Compose<FMatrix,FPreconditioner> > (Pmodp, ApQ);
		std::cerr<<"Preconditioning done\n";
		commentator.stop ("done");
		
	}
											  

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class FMatrix, class IVector,class FVector>
	void RationalSolver<Ring,Field,RandomPrime,WiedemannTraits>::precondition (const Field&                          F,
										   const IMatrix&                        A,
										   BlackboxArchetype<IVector>        *&PAQ,
										   const FMatrix                       *Ap,
										   BlackboxArchetype<FVector>       *&PApQ,
										   const IVector                        &b,
										   IVector                             &Pb,
										   BlackboxArchetype<IVector>          *&P,
										   BlackboxArchetype<IVector>          *&Q) const
	{
		switch (_traits.preconditioner() ) {
			
		case WiedemannTraits::BUTTERFLY:
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<<"ERROR: Butterfly preconditioner not implemented yet. Sorry." << std::endl;		    

		case WiedemannTraits::SPARSE:
			{
				commentator.start ("Constructing sparse preconditioner");
							
				P = new LambdaSparseMatrix<Ring> (_R,Ap->coldim(),Ap->rowdim(),2);
				
				PAQ = new Compose<LambdaSparseMatrix<Ring>, IMatrix> (*P,A);
				
				P->apply(Pb,b);
				
				LambdaSparseMatrix<Field> Pmodp(F,*P);
				
				PApQ = new Compose<LambdaSparseMatrix<Field>, FMatrix> (Pmodp, *Ap);
				
				commentator.stop ("done");
				break;
			}
		
		case WiedemannTraits::TOEPLITZ:
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Toeplitz preconditioner not implemented yet. Sorry." << std::endl;

		case WiedemannTraits::NONE:
			throw PreconditionFailed (__FUNCTION__, __LINE__, "preconditioner is BUTTERFLY, SPARSE, or TOEPLITZ");

		default:
			throw PreconditionFailed (__FUNCTION__, __LINE__, "preconditioner is BUTTERFLY, SPARSE, or TOEPLITZ");
		}



	}
	
			   

	// SPECIALIZATION FOR DIXON
	
	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,DixonTraits>::solve 
	(Vector1& answer, const IMatrix& A, const Vector2& b, const bool old, int maxPrimes, const SolverLevel level ) const {

		SolverReturnStatus status;

		while (maxPrimes > 0){
			switch (A.rowdim() == A.coldim() ? solveNonsingular(answer,A,b,old,1) : SS_SINGULAR) {
					
			case SS_OK:
#ifdef DEBUG_DIXON
				std::cout <<"nonsingular worked\n";
#endif
				return SS_OK;
				break;
					
			case SS_SINGULAR:
#ifdef DEBUG_DIXON
				std::cout<<"switching to singular\n";
#endif
				status = solveSingular(answer,A,b,1,level);
				if (status != SS_FAILED) 
					return status;
				break;
					
			case SS_FAILED:
				std::cout <<"nonsingular failed\n";
				break;
					
			default:
				throw LinboxError ("Bad return value from solveNonsingular");
					
			}
			maxPrimes--;
			if (maxPrimes > 0) chooseNewPrime();
		}
		return SS_FAILED;
	}


	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,DixonTraits>::solveNonsingular 
	(Vector1& answer, const IMatrix& A, const Vector2& b, bool oldMatrix, int maxPrimes) const {

#ifdef DEBUG_DIXON
		cout << "entering nonsingular solver\n";
#endif

		int trials = 0, notfr;

		// history sensitive data for optimal reason
		static const IMatrix* IMP = 0;
		
		static BlasMatrix<Element>* FMP;
		Field *F=NULL;

		do {
			if (trials == maxPrimes) return SS_SINGULAR;			
			if (trials != 0) chooseNewPrime();
			trials++;
#ifdef DEBUG_DIXON
			cout << "_prime: "<<_prime<<"\n";
#endif		       
			typedef typename Field::Element Element;
			typedef typename Ring::Element Integer;

			// checking size of system
			linbox_check(A.rowdim() == A.coldim());
			linbox_check(A.rowdim() == b.size());

			LinBox::integer tmp;		
		
			// if input matrix A is different one.
			if (!oldMatrix) {
		
				//delete IMP;
		
				delete FMP;
		
				IMP = &A;					
		
				F= new Field (_prime);					
		
				FMP = new BlasMatrix<Element>(A.rowdim(),A.coldim());

				typename BlasMatrix<Element>::RawIterator iter_p  = FMP->rawBegin();
				typename IMatrix::ConstRawIterator iter  = A.rawBegin();
				for (;iter != A.rawEnd();++iter,++iter_p)
					F->init(*iter_p, _R.convert(tmp,*iter));
			
				if ( _prime >  Prime(67108863) )
					notfr = MatrixInverse::matrixInverseIn(*F,*FMP);
				else {
					BlasMatrix<Element> *invA = new BlasMatrix<Element>(A.rowdim(),A.coldim());
					BlasMatrixDomain<Field> BMDF(*F);

					BMDF.inv(*invA, *FMP, notfr); //notfr <- nullity
					delete FMP;
					FMP = invA;
// 					cout << "notfr = " << notfr << endl;
// 					cout << "inverse mod p: " << endl;
// 					FMP->write(cout, *F);
				}
			}
			else
				notfr = 0;
		} while (notfr);

		typedef DixonLiftingContainer<Ring,Field,IMatrix,BlasMatrix<Element> > LiftingContainer;
			
		LiftingContainer lc(_R, *F, *IMP, *FMP, b, _prime);		
		RationalReconstruction<LiftingContainer > re(lc);

		if (!re.getRational2(answer)) return SS_FAILED;
		return SS_OK;
	}

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,DixonTraits>::solveSingular 
	(Vector1& answer, const IMatrix& A, const Vector2& b, int maxPrimes, const SolverLevel level) const {

		return monolithicSolve (answer, A, b, false, false, maxPrimes, level);
	}

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,DixonTraits>::findRandomSolution 
	(Vector1& answer, const IMatrix& A, const Vector2& b, int maxPrimes, const SolverLevel level ) const {
		
		return monolithicSolve (answer, A, b, false, true, maxPrimes, level);
	}


	// Most solving is done by the routine below.
	// There used to be one for random and one for deterministic, but they have been merged to ease with 
	//  repeated code (certifying inconsistency, optimization are 2 examples)

	// This function can be thought of broken down into the followin sections:

	// 1. outer loop to cycle through primes
	// 2.   initialize data structures
	// 3.   factor (A mod p) into L.Q.U.P and determine its rank
	// 4.   special case: rank(A mod p) == 0
	// 5.   probabilistically check for consistency: is rank(A mod p) == rank(A|b mod p)?
	// 6.   if ranks are different, return inconsistent, doublechecking over Z first if level >= SL_LASVEGAS
	// 7.   bring A to full row rank form Qt.A, henceforth using only the top (rank) rows
	// 8.   bring A to full column rank form (using Pt if randomSolution=false, random conditioner otherwise)
	// 9.   solve full rank nonsingular system
	// 10.  left-multiply the solution by column conditioner to get solution to original solution
	// 11.  quit if level == SL_MONTECARLO
	// 12.  verify solution, trying new prime if it didn't work
	// 13.  if makeMinDenomCert == true, generate partial certificate for determining solution's minimal denominator

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,DixonTraits>::monolithicSolve 
	(Vector1& answer, const IMatrix& A, const Vector2& b, bool makeMinDenomCert, bool randomSolution,
	 int maxPrimes, const SolverLevel level) const {

		if (level == SL_MONTECARLO && maxPrimes > 1) 
			cout << "WARNING: Even if maxPrimes > 1, SL_MONTECARLO uses just one prime." << endl;
		if (makeMinDenomCert && !randomSolution) 
			cout << "WARNING: Will not compute a certificate of minimal denominator deterministically." << endl;
		if (makeMinDenomCert && level == SL_MONTECARLO) 
			cout << "WARNING: No certificate of min-denominality generated due to  level=SL_MONTECARLO" << endl;

		int trials = 0;
		while (trials < maxPrimes){
			if (trials != 0) chooseNewPrime();
			trials++;
#ifdef DEBUG_DIXON
			cout << "_prime: "<<_prime<<"\n";
#endif		       
			typedef typename Field::Element Element;
			typedef typename Ring::Element Integer;
			typedef DixonLiftingContainer<Ring, Field, 
				BlasBlackbox<Ring>, BlasBlackbox<Field> > LiftingContainer;

			// checking size of system
			linbox_check(A.rowdim() == b.size());
		
			LinBox::integer tmp;
			Integer _rone,_rzero;
			_R.init(_rone,1);
			_R.init(_rzero,0);
			
			Field F (_prime);		
			BlasMatrixDomain<Ring>  BMDI(_R);
			BlasMatrixDomain<Field> BMDF(F);
			BlasApply<Ring> BAR(_R);
			MatrixDomain<Ring> MD(_R);
			VectorDomain<Ring> VDR(_R);

			BlasMatrix<Integer> A_check(A); // used to check answer later

			// TAS_xxx stands for Transpose Augmented System (A|b)t
			// this provides a factorization (A|b) = TAS_Pt . TAS_Ut . TAS_Qt . TAS_Lt
			// such that    
			// - TAS_P . (A|b) . TAS_Q   has nonzero principal minors up to TAS_rank
			// - TAS_Q permutes b to the (TAS_rank)th column of A iff the system is inconsistent mod p

			BlasMatrix<Element>* TAS_factors = new BlasMatrix<Element>(A.coldim()+1, A.rowdim());
			for (size_t i=0;i<A.rowdim();++i)
				for (size_t j=0;j<A.coldim();++j)
					F.init(TAS_factors->refEntry(j,i),_R.convert(tmp,A.getEntry(i,j)));
			for (size_t i=0;i<A.rowdim();++i)
				F.init(TAS_factors->refEntry(A.coldim(),i),_R.convert(tmp,b[i]));
			LQUPMatrix<Field>* TAS_LQUP = new LQUPMatrix<Field>(F, *TAS_factors);
			size_t TAS_rank = TAS_LQUP->getrank();
// 			cout << "tas-rank: " << TAS_rank << endl;
			
			// check consistency. why does getQ return Qt?
			BlasPermutation TAS_P = TAS_LQUP->getP();
			BlasPermutation TAS_Qt = TAS_LQUP->getQ();
			std::vector<size_t> srcRow(A.rowdim()), srcCol(A.coldim()+1);
			std::vector<size_t>::iterator sri = srcRow.begin(), sci = srcCol.begin();
			for (size_t i=0; i<A.rowdim(); i++, sri++) *sri = i;
			for (size_t i=0; i<A.coldim()+1; i++, sci++) *sci = i;
			indexDomain iDom;
			BlasMatrixDomain<indexDomain> BMDs(iDom);
			BMDs.mulin_right(TAS_Qt, srcCol);
			BMDs.mulin_right(TAS_P, srcRow);
#ifdef DEBUG_INC
 			cout << "P takes (0 1 ...) to (";
 			   for (size_t i=0; i<A.rowdim(); i++) cout << srcRow[i] << ' '; cout << ')' << endl;
 			cout << "Q takes (0 1 ...) to (";
 			   for (size_t i=0; i<A.coldim()+1; i++) cout << srcCol[i] << ' '; cout << ')' << endl;
#endif
			
			bool appearsInconsistent = (srcCol[TAS_rank-1] == A.coldim());
			size_t rank = TAS_rank - (appearsInconsistent ? 1 : 0);
#ifdef DIXON_DEBUG
			cout << "TAS_rank, rank: " << TAS_rank << ' ' << rank << endl;
#endif

			if (rank == 0) { 
				delete TAS_LQUP;
				delete TAS_factors;
				//special case when A = 0, mod p. dealt with to avoid crash later
				bool aEmpty = true;
				if (level >= SL_LASVEGAS) { // in monte carlo, we assume A is actually empty
					typename BlasMatrix<Integer>::RawIterator iter = A_check.rawBegin();
					for (; aEmpty && iter != A_check.rawEnd(); ++iter)
						aEmpty &= _R.isZero(*iter);
				}
				if (aEmpty) {
					for (size_t i=0; i<b.size(); i++)
						if (!_R.areEqual(b[i], _rzero)) {
							if (level >= SL_CERTIFIED) {
								lastCertificate.clearAndResize(b.size());
								_R.assign(lastCertificate.numer[i], _rone);
							}
							return SS_INCONSISTENT;
						}
					// both A and b are all zero.
					for (size_t i=0; i<answer.size(); i++) {
						answer[i].first = _rzero;
						answer[i].second = _rone;
					}
					if (level >= SL_LASVEGAS)
						_R.init(lastCertifiedDenFactor, 1);
					if (level == SL_CERTIFIED) {
						_R.init(lastZBNumer, 0);
						lastCertificate.clearAndResize(b.size());
					}
					return SS_OK;
				}
				// so a was empty mod p but not over Z.
				continue; //try new prime
			}
	
			BlasMatrix<Element>* Atp_minor_inv;

			if ((appearsInconsistent || randomSolution == false) && level > SL_MONTECARLO) {
				// take advantage of the (LQUP)t factorization to compute
				// an inverse to the leading minor of (TAS_P . (A|b) . TAS_Q)t 
				
				Atp_minor_inv = new BlasMatrix<Element>(rank, rank);
				typename BlasMatrix<Element>::RawIterator iter = Atp_minor_inv->rawBegin();
				for (; iter != Atp_minor_inv->rawEnd(); iter++)
					F.init(*iter, 0);

				FFLAPACK::InverseInFromLQUP(F, rank, TAS_factors->getPointer(), A.rowdim(), 
							    TAS_Qt.getPointer(), Atp_minor_inv->getPointer(), rank);
			}

			delete TAS_LQUP;
			delete TAS_factors;

			if (appearsInconsistent && level <= SL_MONTECARLO)
				return SS_INCONSISTENT;

			if (appearsInconsistent) {
				std::vector<Integer> zt(rank);
				for (size_t i=0; i<rank; i++)
					_R.assign(zt[i], A.getEntry(srcRow[rank], srcCol[i]));

				BlasMatrix<Integer> At_minor(rank, rank);
				for (size_t i=0; i<rank; i++)
					for (size_t j=0; j<rank; j++)
						_R.assign(At_minor.refEntry(j, i), A.getEntry(srcRow[i], srcCol[j]));
#ifdef DEBUG_INC
 				At_minor.write(cout << "At_minor:" << endl, _R);
 				Atp_minor_inv->write(cout << "Atp_minor_inv:" << endl, F);
				cout << "zt: "; for (size_t i=0; i<rank; i++) cout << zt[i] <<' '; cout << endl;
#endif
				BlasBlackbox<Ring>  BBAt_minor(_R, At_minor);
				BlasBlackbox<Field> BBAtp_minor_inv(F, *Atp_minor_inv);
				LiftingContainer lc(_R, F, BBAt_minor, BBAtp_minor_inv, zt, _prime);
				
				RationalReconstruction<LiftingContainer > re(lc);
				
				Vector1 short_answer(rank);
				
				if (!re.getRational2(short_answer)) return SS_FAILED; 

				delete Atp_minor_inv;
				
				VectorFraction<Ring> cert(_R, short_answer);
				cert.numer.resize(b.size());
				_R.subin(cert.numer[rank], cert.denom);
				_R.init(cert.denom, 1);
				BMDI.mulin_left(cert.numer, TAS_P);
#ifdef DEBUG_INC
 				cert.write(cout << "cert:") << endl;
#endif

				bool certifies = true; //check certificate
				std::vector<Integer> certnumer_A(A.coldim());
				BAR.applyVTrans(certnumer_A, A_check, cert.numer);
				typename std::vector<Integer>::iterator cai = certnumer_A.begin();
				for (size_t i=0; certifies && i<A.coldim(); i++, cai++) 
					certifies &= _R.isZero(*cai);
				if (certifies) { 
					if (level == SL_CERTIFIED) lastCertificate.copy(cert);
					return SS_INCONSISTENT;
				}
				continue; // try new prime. analogous to u.A12 != A22 in Muld.+Storj.
			}

			// we now know system is consistent mod p.
			BlasMatrix<Integer> A_minor(rank, rank);    // -- will have the full rank minor of A
			BlasMatrix<Element> *Ap_minor_inv;          // -- will have inverse mod p of A_minor
			BlasMatrix<Integer> *P, *B;                 // -- only used in random case 

			if (!randomSolution) {
				// use shortcut - transpose Atp_minor_inv to get Ap_minor_inv
				Element _rtmp; 
				Ap_minor_inv = Atp_minor_inv;
				for (size_t i=0; i<rank; i++)
					for (size_t j=0; j<i; j++) {
						Ap_minor_inv->getEntry(_rtmp, i, j);
						Ap_minor_inv->setEntry(i, j, Ap_minor_inv->refEntry(j, i));
						Ap_minor_inv->setEntry(j, i, _rtmp);
					}

				// permute original entries into A_minor
				for (size_t i=0; i<rank; i++)
					for (size_t j=0; j<rank; j++)
						_R.assign(A_minor.refEntry(i, j), A_check.getEntry(srcRow[i], srcCol[j]));
			}
			else {
				P = new BlasMatrix<Integer>(A.coldim(),rank);	
				B = new BlasMatrix<Integer>(rank,A.coldim());
				BlasMatrix<Element> Ap_minor(rank, rank);
				Ap_minor_inv = new BlasMatrix<Element>(rank, rank);
				int nullity;
				
				for (size_t i=0; i<rank; i++)
					for (size_t j=0; j<A.coldim(); j++)
						_R.assign(B->refEntry(i, j), A_check.getEntry(srcRow[i], j));
				
				do { // O(1) loops of this preconditioner expected
					// compute P a n*r random matrix of entry in [0,1]
					typename BlasMatrix<Integer>::RawIterator iter;
					for (iter = P->rawBegin(); iter != P->rawEnd(); ++iter) {
						if (rand() > RAND_MAX/2) 
							_R.assign(*iter, _rone);
						else
							_R.assign(*iter, _rzero);
					}

					// compute A_minor = B.P
					MD.mul(A_minor, *B, *P);
					
					// set Ap_minor = A_minor mod p, try to compute inverse
					for (size_t i=0;i<rank;++i)
						for (size_t j=0;j<rank;++j)
							F.init(Ap_minor.refEntry(i,j),
							       _R.convert(tmp,A_minor.getEntry(i,j)));
					
					BMDF.inv(*Ap_minor_inv, Ap_minor, nullity); 
				} while (nullity > 0);
			}

			// Compute newb = (TAS_P.b)[0..(rank-1)]
			std::vector<Integer> newb(b);
			BMDI.mulin_right(TAS_P, newb);
			newb.resize(rank);

			BlasBlackbox<Ring>  BBA_minor(_R,A_minor);
			BlasBlackbox<Field> BBA_inv(F,*Ap_minor_inv);
			LiftingContainer lc(_R, F, BBA_minor, BBA_inv, newb, _prime);

#ifdef DEBUG_DIXON
			std::cout<<"length of lifting: "<<lc.length()<<std::endl;
#endif
			RationalReconstruction<LiftingContainer > re(lc);

			Vector1 short_answer(rank);

			if (!re.getRational2(short_answer)) 
				return SS_FAILED;

			VectorFraction<Ring> answer_to_vf(_R, short_answer);

			if (!randomSolution) {
				// short_answer = TAS_Q * short_answer
				answer_to_vf.numer.resize(A.coldim()+1,_rzero);
				BMDI.mulin_left(answer_to_vf.numer, TAS_Qt);
				answer_to_vf.numer.resize(A.coldim());
			}
			else {
				// short_answer = P * short_answer
				typename Vector<Ring>::Dense newNumer(A.coldim());
				BAR.applyV(newNumer, *P, answer_to_vf.numer);
				answer_to_vf.numer = newNumer;
			}

			if (level >= SL_LASVEGAS) { //check consistency

				std::vector<Integer> A_times_xnumer(b.size());
	
				BAR.applyV(A_times_xnumer, A_check, answer_to_vf.numer);
				
				Integer tmpi;
				
				typename Vector2::const_iterator ib = b.begin();
				typename std::vector<Integer>::iterator iAx = A_times_xnumer.begin();
				int thisrow = 0;
				bool needNewPrime = false;
				
				for (; !needNewPrime && ib != b.end(); iAx++, ib++, thisrow++)
					if (!_R.areEqual(_R.mul(tmpi, *ib, answer_to_vf.denom), *iAx)) {
						// should attempt to certify inconsistency now
						// as in "if [A31 | A32]y != b3" of step (4)
						needNewPrime = true;
					}
				
				if (needNewPrime) {
					delete Ap_minor_inv;          
					if (randomSolution) {delete P; delete B;}
					continue; //go to start of main loop
				}
			}
			
			answer_to_vf.toFVector(answer);
			
			if (makeMinDenomCert && level >= SL_LASVEGAS && randomSolution) {
				// To make this certificate we solve with the same matrix as to get the 
				// solution, except transposed.

				Integer _rtmp;
				Element _ftmp;
				for (size_t i=0; i<rank; i++)
					for (size_t j=0; j<i; j++) {
						Ap_minor_inv->getEntry(_ftmp, i, j);
						Ap_minor_inv->setEntry(i, j, Ap_minor_inv->refEntry(j, i));
						Ap_minor_inv->setEntry(j, i, _ftmp);
					}

				for (size_t i=0; i<rank; i++)
					for (size_t j=0; j<i; j++) {
						A_minor.getEntry(_rtmp, i, j);
						A_minor.setEntry(i, j, A_minor.refEntry(j, i));
						A_minor.setEntry(j, i, _rtmp);
					}

				// we then try to create a partial certificate
				// the correspondance with Algorithm MinimalSolution from Mulders/Storjohann:
				// paper | here
				// P     | TAS_P
				// Q     | transpose of TAS_Qt
				// B     | *B (== TAS_P . A,  only top #rank rows)
				// c     | newb (== TAS_P . b,   only top #rank rows)
				// P     | P
				// q     | q
				// U     | {0, 1}
				// u     | u
				// zhat  | lastCertificate
				
				// we multiply the certificate by TAS_Pt at the end
				// so it corresponds to b instead of newb

				//q in {0, 1}^rank
				std::vector<Integer> q(rank);
				typename std::vector<Integer>::iterator q_iter;

				bool allzero;
				do {
					allzero = true;
					for (q_iter = q.begin(); q_iter != q.end(); ++q_iter) {
						if (rand() > RAND_MAX/2) {
							(*q_iter) = _rone;
							allzero = false;
						}
						else
							(*q_iter) = _rzero;
					}
				} while (allzero);

				LiftingContainer lc2(_R, F, BBA_minor, BBA_inv, q, _prime);
				RationalReconstruction<LiftingContainer> re(lc2);
				Vector1 u(rank);
				if (!re.getRational2(u)) return SS_FAILED;

				// remainder of code does   z <- denom(partial_cert . Mr) * partial_cert * Qt 
				VectorFraction<Ring> u_to_vf(_R, u);
				std::vector<Integer> uB(A.coldim());
				BAR.applyVTrans(uB, *B, u_to_vf.numer);

// 				cout << "BP: ";
// 				A_minor.write(cout, _R) << endl;
// 				cout << "q: ";
// 				for (size_t i=0; i<rank; i++) cout << q[i]; cout << endl;
// 				u_to_vf.write(cout  << "u: ") << endl;

				Integer numergcd = _rzero;
				vectorGcdIn(numergcd, _R, uB);

				// denom(partial_cert . Mr) = partial_cert_to_vf.denom / numergcd
				VectorFraction<Ring> z(_R, b.size()); //new constructor
				u_to_vf.numer.resize(A.rowdim());

				//TransposedBlasMatrix<BlasPermutation> Ap_Q(Ap_Qt);
				BMDI.mul(z.numer, u_to_vf.numer, TAS_P);

				z.denom = numergcd;

// 				z.write(cout << "z: ") << endl;

				if (level >= SL_CERTIFIED)
					lastCertificate.copy(z);
				
				// output new certified denom factor
				Integer znumer_b, zbgcd;
				VDR.dotprod(znumer_b, z.numer, b);
				_R.gcd(zbgcd, znumer_b, z.denom);
				_R.div(lastCertifiedDenFactor, z.denom, zbgcd);

				if (level >= SL_CERTIFIED) 
					_R.div(lastZBNumer, znumer_b, zbgcd);
			}

			delete Ap_minor_inv;          
			if (randomSolution) {delete P; delete B;}

			// done making certificate, lets blow this popstand
			return SS_OK;
		}
		return SS_FAILED; //all primes were bad
	}

} //end of namespace LinBox
#endif

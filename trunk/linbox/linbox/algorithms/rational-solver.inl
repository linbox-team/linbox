/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/algorithms/rational-solver.inl
 * Copyright (C) 2004 Pascal Giorgi
 *
 * Written by Pascal Giorgi pascal.giorgi@ens-lyon.fr
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
//#define DEBUG_DIXON_MATRICES
//#define DEBUG_DIXON_RANDOM

namespace LinBox {
	
       
	// SPECIALIZATION FOR WIEDEMANN 	

	// note: if Vector1 != Vector2 compilation of solve or solveSingluar will fail (via an invalid call to sparseprecondition)!
	// maybe they should not be templated separately, or sparseprecondition should be rewritten

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	typename RationalSolver<Ring,Field,RandomPrime,WiedemannTraits>::ReturnStatus 
	RationalSolver<Ring,Field,RandomPrime,WiedemannTraits>::solve (Vector1& answer,
								       const IMatrix& A,
								       const Vector2& b,
								       const bool old=false,
								       int maxPrimes) const {

		typename RationalSolver<Ring,Field,RandomPrime,WiedemannTraits>::ReturnStatus status=FAILED;

		switch (solveNonsingular(answer,A,b)) {

		case OK:
			status=OK;
			break;

		case SINGULAR:
			std::cerr<<"switching to singular\n";
			status=solveSingular(answer,A,b);
			break;

		case FAILED:			
			break;

		default:
			throw LinboxError ("Bad return value from solveNonsingular");
			
		}

		return status;
	}

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	typename RationalSolver<Ring,Field,RandomPrime,WiedemannTraits>::ReturnStatus 
	RationalSolver<Ring,Field,RandomPrime,WiedemannTraits>::solveNonsingular (Vector1& answer,
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
		unsigned long issingular=SINGULARITY_THRESHOLD; 			
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
			return SINGULAR;
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
			
			re.getRational(answer);
					
			return OK;
		}
	}       

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	typename RationalSolver<Ring,Field,RandomPrime,WiedemannTraits>::ReturnStatus 
	RationalSolver<Ring,Field,RandomPrime,WiedemannTraits>::solveSingular (Vector1& answer,
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

			return BAD_PRECONDITIONER;
		}				      		      				
		else {	
			
			MinPoly.erase(MinPoly.begin());
			
			typedef WiedemannLiftingContainer<Ring, Field, IPrecondMatrix, FPrecondMatrix, FPolynomial> LiftingContainer;
			std::cerr<<"before lc\n";
			LiftingContainer lc(_R, *F, *PAQ, *PApQ, MinPoly, Pb, _prime);
			std::cerr<<"constructing lifting container of length: "<<lc.length()<<std::endl;
			
			RationalReconstruction<LiftingContainer> re(lc,_R,2);
			
			re.getRational(answer); 


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
			
			return OK;
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
	typename RationalSolver<Ring,Field,RandomPrime,DixonTraits>::ReturnStatus 
	RationalSolver<Ring,Field,RandomPrime,DixonTraits>::solve (Vector1& answer,
								   const IMatrix& A,
								   const Vector2& b,
								   const bool old,
								   int maxPrimes) const {

		typename RationalSolver<Ring,Field,RandomPrime,DixonTraits>::ReturnStatus status=FAILED;
		int inconsistencyCount = 0;
		while (maxPrimes > 0){
			switch (solveNonsingular(answer,A,b,old,1)) {
					
			case OK:
				//cout<<"Okay"<<endl;
				return OK;
				break;
					
			case SINGULAR:
				//cout<<"Singular"<<endl;
#ifdef DEBUG_DIXON
				std::cout<<"switching to singular\n";
#endif
				status=solveSingular(answer,A,b,1);
				if (status==OK) return OK;
				if (status==INCONSISTENT) inconsistencyCount++;
				break;
					
			case FAILED:
				//cout<<"Failed"<<endl;
				break;
					
			default:
				//cout<<"Error"<<endl;
				throw LinboxError ("Bad return value from solveNonsingular");
					
			}

			maxPrimes--;
			if (maxPrimes > 0)
				_prime = _genprime.randomPrime();
		}
		if (inconsistencyCount > 0) return INCONSISTENT;
		return FAILED;
	}


	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	typename RationalSolver<Ring,Field,RandomPrime,DixonTraits>::ReturnStatus 
	RationalSolver<Ring,Field,RandomPrime,DixonTraits>::solveNonsingular (Vector1& answer,
									      const IMatrix& A,
									      const Vector2& b,
									      bool oldMatrix = false,
									      int maxPrimes) const {
		int trials = 0;

		// history sensitive data for optimal reason
		static const IMatrix* IMP = 0;
		
		static BlasMatrix<Element>* FMP;
		Field *F=NULL;
		int notfr;

		do {
			
			if (trials != 0)
				_prime = _genprime.randomPrime();
			trials++;
#ifdef DEBUG_DIXON
			cout << "_prime: "<<_prime<<"\n";
#endif		       
			typedef typename Field::Element Element;
			typedef typename Ring::Element Integer;

			// checking if matrix is square
			linbox_check(A.rowdim() == A.coldim());

			// checking size of system
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
					int nullity;
					//cout << "FMP:\n";
					//FMP->write(cout, *F, true);
					BMDF.inv(*invA, *FMP, nullity);
					//cout << "invA:\n";
					//invA->write(cout, *F, true);
					notfr = nullity;
					delete FMP;
					FMP = invA;  
				}
			}
			else
				notfr = 0;
		} while (notfr && trials < maxPrimes);

		// not full rank
		if(notfr) 
			return SINGULAR;

		typedef DixonLiftingContainer<Ring,Field,IMatrix,BlasMatrix<Element> > LiftingContainer;
			
		LiftingContainer lc(_R, *F, *IMP, *FMP, b, _prime);
			
		RationalReconstruction<LiftingContainer > re(lc);
			
		if (!re.getRational(answer)) {
			std::cerr << "ERROR: lifting failed. Are you using <double> ring with large norm?" << endl;
			return FAILED;
		}
			
		return OK;
	}


	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	typename RationalSolver<Ring,Field,RandomPrime,DixonTraits>::ReturnStatus 
	RationalSolver<Ring,Field,RandomPrime,DixonTraits>::solveSingular (Vector1& answer,
									   const IMatrix& A,
									   const Vector2& b,
									   int maxPrimes) const {
#ifdef DEBUG_DIXON
		cout << "switched to singular" << endl;
#endif
		int inconsistencyCount = 0;
		int trials = 0;
		while (trials < maxPrimes){
			
			if (trials != 0)
				_prime = _genprime.randomPrime();
			trials++;
#ifdef DEBUG_DIXON
			cout << "_prime: "<<_prime<<"\n";
#endif		       

			typedef typename Field::Element Element;
			typedef typename Ring::Element Integer;

			// checking size of system
			linbox_check(A.rowdim() == b.size());
		
			LinBox::integer tmp;

			Field F (_prime);
		
			// set M=A and compute Ap = A mod p
			BlasMatrix<Element> *Ap = new BlasMatrix<Element>(A.rowdim(),A.coldim());
			BlasMatrix<Integer> M(A);

			typename BlasMatrix<Element>::RawIterator iter_p = Ap->rawBegin();
			typename BlasMatrix<Integer>::RawIterator iter   = M.rawBegin();
			for (; iter != M.rawEnd(); ++iter, ++iter_p)
				F.init(*iter_p,_R.convert(tmp,*iter));

			// compute the LQUP factorization of Ap			
			LQUPMatrix<Field>* LQUP = new LQUPMatrix<Field>(F,*Ap);

			// create Qt= Transpose(Q) 
			BlasPermutation Q;
			Q=LQUP->getQ(); 
			size_t rank= LQUP->getrank();
			//TransposedBlasMatrix<BlasPermutation> Qt(Q);
			delete LQUP;
			delete Ap;

	
			// check if system is consistent
			Ap= new BlasMatrix<Element>(A.rowdim(),A.coldim()+1);
			for (size_t i=0;i<A.rowdim();++i)
				for (size_t j=0;j<A.coldim();++j)
					F.init(Ap->refEntry(i,j),_R.convert(tmp,A.getEntry(i,j)));
			for (size_t i=0;i<A.rowdim();++i)
				F.init(Ap->refEntry(i,A.coldim()),_R.convert(tmp,b[i]));
						
			LQUP= new LQUPMatrix<Field>(F,*Ap);
	
#ifdef DEBUG_DIXON
			std::cout<<"rank of A = "<<rank<<", rank of [A|b] is "<<LQUP->getrank()<<endl;
#endif

			if (rank != LQUP->getrank()){
#ifdef DEBUG_DIXON
				std::cout<<"system is inconsistent mod "<<_prime<<endl;
#endif
				inconsistencyCount++;
				if (trials<maxPrimes) continue;
				return INCONSISTENT;
			}
#ifdef DEBUG_DIXON
			else
				std::cout<<"system is consistent\n";
#endif

			if (rank == 0) { 
				//A = b = 0, mod p
				int bEmpty = true;
				int n = b.size();
				Integer zero;
				Integer one;
				_R.init(zero, 0);
				_R.init(one, 1);
				for (int i=0; i<n; i++)
					bEmpty &= _R.areEqual(b[i], zero);
				if (bEmpty) {
					for (int i=0; i<n; i++) {
						answer[i].first = zero;
						answer[i].second = one;
					}
					return OK;
				}
				if (trials < maxPrimes) continue;
				if (inconsistencyCount > 0) return INCONSISTENT;
				return FAILED;
			}

			// compute M = Qt.A 
			BlasMatrixDomain<Ring>  BMDI(_R);
			BlasMatrixDomain<Field> BMDF(F);
			BlasApply<Ring> BAR(_R);
	
			BMDI.mulin_right(Q,M);
					
			// set Ap = Qt.A mod p
			Ap =  new BlasMatrix<Element>(rank,M.coldim());
			for (size_t i=0;i<rank;++i)
				for (size_t j=0;j<M.coldim();++j)
					F.init(Ap->refEntry(i,j),_R.convert(tmp,M.getEntry(i,j)));

			// compute the LQUP factorization of Ap
			LQUP = new LQUPMatrix<Field>(F,*Ap);
       
			// create Pt= Transpose(P)
			BlasPermutation P;
			P=LQUP->getP(); 
			TransposedBlasMatrix<BlasPermutation> Pt(P);

			// compute M= Qt.A.Pt
			BMDI.mulin_left(M,Pt);

			delete Ap;
			delete LQUP;
		
			// set A_minor= the rank*rank leading minor of M
			BlasMatrix<Integer> A_minor (M,0,0,rank,rank);

			// set Ap = A_minor mod p
			Ap =  new BlasMatrix<Element>(rank,rank);
			for (size_t i=0;i<rank;++i)
				for (size_t j=0;j<rank;++j)
					F.init(Ap->refEntry(i,j),_R.convert(tmp,A_minor.getEntry(i,j)));										
		

			// compute Ainv= A^(-1) mod p
			BlasMatrix<Element> Ainv(rank,rank);
			BMDF.invin(Ainv,*Ap); 
			delete Ap;
		
			// Compute Qtb= Qt.b
			std::vector<Integer> Qtb(b);
			BMDI.mul(Qtb,Q,b);
			Qtb.resize(rank);		

			typedef DixonLiftingContainer<Ring, Field, 
				BlasBlackbox<Ring>, BlasBlackbox<Field> > LiftingContainer;
			BlasBlackbox<Ring>  BBA_minor(_R,A_minor);
			BlasBlackbox<Field> BBA_inv(F,Ainv);
			LiftingContainer lc(_R, F, BBA_minor, BBA_inv, Qtb, _prime);

#ifdef DEBUG_DIXON
			std::cout<<"length of lifting: "<<lc.length()<<std::endl;
#endif
			
			RationalReconstruction<LiftingContainer > re(lc);
		
			Vector1 short_answer(rank);

			if (!re.getRational(short_answer)) {
				std::cerr << "ERROR: lifting failed. Maybe <double> field or ring with large norm?";
				if (trials < maxPrimes) {
					std::cerr << " Trying another prime" << endl;
					continue;
				}
				std::cerr << " Failed." << endl;
				return FAILED;
			}

			// test solution of minor
			/*
			  Integer lden;
			  _R.init(lden,1);
		
			  std::vector<std::pair<Integer,Integer> >::iterator pair_iter;
		
			  for (pair_iter= short_answer.begin(); pair_iter != short_answer.end(); ++pair_iter)
			  _R. lcm (lden, lden, pair_iter->second);
		
			  std::vector<Integer> x(rank);
			  std::vector<Integer>::iterator px=x.begin();
		
		
			  for (pair_iter= short_answer.begin(); pair_iter != short_answer.end(); ++pair_iter, ++px){
			  _R. mul (*px, pair_iter->first, lden);
			  _R. divin (*px, pair_iter->second);
			  }

  
			  std::vector<Integer> r(x.size());
	
			  BAR.applyV(r,A_minor,x);
  
			  for (size_t i=0;i< r.size();++i) {
			  if (!_R.isDivisor(r[i], lden)){
			  #ifdef DEBUG_DIXON
			  cout << "not div: "<<r[i]<<" "<<lden<<"\n";
			  #endif
			  cerr << "ring error in normalization" << endl;
			  if (trials < maxPrimes) continue;
			  return FAILED;
			  }
			  _R.divin(r[i],lden);
			  }


			  VectorDomain<Ring> VDR(_R);
			  if (!VDR.areEqual(r,Qtb)) {
			  #ifdef DEBUG_DIXON
			  std::cout<<"\n nonsingular system failed\n";
			  #endif
			  // this branch gets called when reduced system is consistent but not of 
			  // the same rank as the original, i believe
			  if (trials < maxPrimes) continue;
			  if (inconsistencyCount > 0) return INCONSISTENT;
			  return FAILED;
			  }
			  #ifdef DEBUG_DIXON
			  else
			  std::cout<<"\n nonsingular system solved\n";
			  #endif
		

			  // end of checking
			  */
			
			
			Integer _rone,_rzero;
			_R.init(_rone,1);
			_R.init(_rzero,0);
			
			short_answer.resize(A.coldim(),std::pair<Integer,Integer>(_rzero,_rone));
		
			for (int i=Pt.getMatrix().getOrder()-1;i >=0 ;--i){
				if (*(Pt.getMatrix().getPointer()+i) > (size_t) i)
					std::swap( short_answer[i], short_answer[*(Pt.getMatrix().getPointer()+i)]);
			}
			answer=short_answer;

			// check answer

			return OK;
		}
		//should not get here
		cerr << "singular dixon failure\n";
		return FAILED;
	}



	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	typename RationalSolver<Ring,Field,RandomPrime,DixonTraits>::ReturnStatus 
	RationalSolver<Ring,Field,RandomPrime,DixonTraits>::findRandomSolution (Vector1& answer,
										const IMatrix& A,
										const Vector2& b,
										int maxPrimes) const {
		int trials = 0;
		int inconsistencyCount = 0;
		while (trials < maxPrimes){ 
			if (trials != 0)
				_prime = _genprime.randomPrime();
			trials++;
#ifdef DEBUG_DIXON
			cout << "_prime: "<<_prime<<"\n";
#endif		       


			typedef typename Field::Element Element;
			typedef typename Ring::Element Integer;

			// checking size of system
			linbox_check(A.rowdim() == b.size());
		
			LinBox::integer tmp;

			Field F (_prime);
		
			// set M=A and compute Ap = A mod p
			BlasMatrix<Element> *Ap = new BlasMatrix<Element>(A.rowdim(),A.coldim());
			BlasMatrix<Integer> M(A);
			BlasMatrix<Integer> oldA(A); //for checking later

			typename BlasMatrix<Element>::RawIterator iter_p = Ap->rawBegin();
			typename BlasMatrix<Integer>::RawIterator iter   = M.rawBegin();
			for (; iter != M.rawEnd(); ++iter, ++iter_p)
				F.init(*iter_p,_R.convert(tmp,*iter));

			// compute the LQUP factorization of Ap			
			LQUPMatrix<Field>* LQUP = new LQUPMatrix<Field>(F,*Ap);


			// create Qt= Transpose(Q) 
			BlasPermutation Q;
			Q=LQUP->getQ(); 
			size_t rank= LQUP->getrank();
			//TransposedBlasMatrix<BlasPermutation> Qt(Q);
			delete LQUP;
		

			// check if system is consistent
			{
				BlasMatrix<Element> Ap= BlasMatrix<Element>(A.rowdim(),A.coldim()+1);
				for (size_t i=0;i<A.rowdim();++i)
					for (size_t j=0;j<A.coldim();++j)
						F.init(Ap.refEntry(i,j),_R.convert(tmp,A.getEntry(i,j)));
				for (size_t i=0;i<A.rowdim();++i)
					F.init(Ap.refEntry(i,A.coldim()),_R.convert(tmp,b[i]));
				
				LQUP= new LQUPMatrix<Field>(F,Ap);
				
#ifdef DEBUG_DIXON
				std::cout<<"rank of A = "<<rank<<", rank of [A|b] is "<<LQUP->getrank()<<endl;
#endif
				
				if (rank != LQUP->getrank()){
#ifdef DEBUG_DIXON
					std::cout<<"system is inconsistent\n";
#endif
					inconsistencyCount++;
					if (trials<maxPrimes) continue;
					return INCONSISTENT;
				}
#ifdef DEBUG_DIXON
				else
					std::cout<<"system is consistent\n";
#endif
			}

			if (rank == 0) { 
				//A = b = 0, mod p
				int bEmpty = true;
				int n = b.size();
				Integer zero;
				Integer one;
				_R.init(zero, 0);
				_R.init(one, 1);
				for (int i=0; i<n; i++)
					bEmpty &= _R.areEqual(b[i], zero);
				if (bEmpty) {
					for (int i=0; i<n; i++) {
						answer[i].first = zero;
						answer[i].second = one;
					}
					return OK;
				}
				if (trials < maxPrimes) continue;
				if (inconsistencyCount > 0) return INCONSISTENT;
				return FAILED;
			}

			// compute M = Qt.A 
			BlasApply<Ring> BAR(_R);
			MatrixDomain<Ring> MD(_R);
			BlasMatrixDomain<Ring>  BMDI(_R);
			BlasMatrixDomain<Field> BMDF(F);		
			BMDI.mulin_right(Q,M);

			BlasMatrix<Integer> A_minor(rank,rank);
			BlasMatrix<Integer> Mr(M,0,0,rank,M.coldim());
			typename Ring::RandIter G(_R,2,0);
			BlasMatrix<Integer> P(M.coldim(),rank);	
				
			do { // O(1) loops of this preconditioner expected
				//	cout << "P";
				delete Ap;
				// compute P a n*r random matrix of entry in [0,1]
				iter    = P.rawBegin();		
			
				for (; iter != P.rawEnd(); ++iter)
					G.random(*iter);
					
				//cout<<"Mr\n";;
				//Mr.write(cout,_R)<<endl;
				
				// compute A_minor= M[1..rank,...] .P (was BMDI)
				MD.mul(A_minor,Mr,P);
				
				//cout<<"A_minor\n";
				//A_minor.write(cout,_R)<<endl;

				// set Ap = A_minor mod p
				Ap =  new BlasMatrix<Element>(rank,rank);
				for (size_t i=0;i<rank;++i)
					for (size_t j=0;j<rank;++j)
						F.init(Ap->refEntry(i,j),_R.convert(tmp,A_minor.getEntry(i,j)));
			} while (BMDF.rank(*Ap) != rank);

#ifdef DEBUG_DIXON_RANDOM		
			cout << "P:\n";
			P.write(cout, _R, true);
			cout << "A_minor:\n";
			A_minor.write(cout, _R, true);
#endif
			
			// compute Ainv= A^(-1) mod p
			BlasMatrix<Element> Ainv(rank,rank);
			BMDF.invin(Ainv,*Ap);
			delete Ap;

#ifdef DEBUG_DIXON_RANDOM		
			cout << "inverse of A_minor, mod p:\n";
			Ainv.write(cout, F, true);
#endif

			// Compute Qtb= Qt.b
			std::vector<Integer> Qtb(b);
			BMDI.mul(Qtb,Q,b);
			Qtb.resize(rank);	

#ifdef DEBUG_DIXON_RANDOM		
			VectorDomain<Ring> VDR(_R);
			cout << "Qtb:\n";
			VDR.write(cout, Qtb);
			cout << "\n";
#endif
		
			typedef DixonLiftingContainer<Ring, Field, BlasBlackbox<Ring>, BlasBlackbox<Field> > LiftingContainer;

			BlasBlackbox<Ring> BBA_minor(_R,A_minor) ;
			BlasBlackbox<Field> BBA_inv(F,Ainv) ;
			LiftingContainer lc(_R, F, BBA_minor, BBA_inv , Qtb, _prime);
#ifdef DEBUG_DIXON
			std::cout<<"length of lifting: "<<lc.length()<<std::endl;
#endif
			
			RationalReconstruction<LiftingContainer > re(lc);
		
			Vector1 short_answer(rank);

			if (!re.getRational(short_answer)) {
				std::cerr << "ERROR: lifting failed. Maybe <double> field or ring with large norm?";
				if (trials < maxPrimes) {
					std::cerr << " Trying another prime" << endl;
					continue;
				}
				std::cerr << " Failed." << endl;
				return FAILED;
			}

			Integer lden;
			_R.init(lden,1);
		
			typename std::vector<std::pair<Integer,Integer> >::iterator pair_iter;

#ifdef DEBUG_DIXON_RANDOM
			cout << "Solution: [";
			for (pair_iter= short_answer.begin(); pair_iter != short_answer.end(); ++pair_iter) {
				if (pair_iter != short_answer.begin()) cout << ", ";
				cout << pair_iter->first << "/" << pair_iter->second;
			}
			cout << "]" << endl;
#endif
		
			for (pair_iter= short_answer.begin(); pair_iter != short_answer.end(); ++pair_iter)
				_R. lcm (lden, lden, pair_iter->second);
		
			std::vector<Integer> x(rank);
			typename std::vector<Integer>::iterator px=x.begin();
		
			for (pair_iter= short_answer.begin(); pair_iter != short_answer.end(); ++pair_iter, ++px){
				_R. mul (*px, pair_iter->first, lden);
				_R. divin (*px, pair_iter->second);
			}

			Integer _rone,_rzero;
			_R.init(_rone,1);
			_R.init(_rzero,0);
			
			//cout<< "P is "<<P.coldim()<<" by "<<P.rowdim()<<", x has length "<<x.size()<<"\n";
			
			//apply conditioner to get solution to original problem
			std::vector<Integer> Px(P.rowdim());
			
			BAR.applyV(Px,P,x);
			
			//now need to check that A*Px = b*lden
			std::vector<Integer> r(b.size());
			BAR.applyV(r,oldA,Px);
			typename std::vector<Integer>::iterator pr=r.begin();
			typename Vector2::const_iterator pb=b.begin();
			bool okay = true;
			Integer tmp2;
			for (; okay && pr != r.end(); pb++, pr++)
				okay &= _R.areEqual(*pr, _R.mul(tmp2, lden, *pb));
			
#ifdef DEBUG_DIXON
			cout << "lden: "<< lden<<"\n";
#ifdef DEBUG_DIXON_MATRICES
			cout << "M is \n";
			oldA.write(cout, _R);

			std::cout<<"Px:\n";
			for (int i=0;i<Px.size();++i)
				std::cout<<Px[i]<<",";
			std::cout<<std::endl;		
			std::cout<<"b:\n";
			for (int i=0;i<b.size();++i)
				std::cout<<b[i]<<",";
			std::cout<<std::endl;		
			std::cout<<"r:\n";
			for (int i=0;i<r.size();++i)
				std::cout<<r[i]<<",";
			std::cout<<std::endl;		
#endif
#endif
			
			if (!okay){
#ifdef DEBUG_DIXON
				std::cout<<"nonsingular system FAILED\n";
#endif
				// this branch gets called when reduced system is consistent but not of 
				// the same rank as the original, i believe
				if (trials < maxPrimes) continue;
				if (inconsistencyCount > 0) return INCONSISTENT;
				return FAILED;
			}
			
			short_answer.resize(A.coldim());
			pair_iter = short_answer.begin();
			px = Px.begin();
			for (; pair_iter != short_answer.end(); pair_iter++, px++)
				(*pair_iter) = std::pair<Integer, Integer>((*px), lden);
			
			answer=short_answer;
			return OK;
			

		}
		return OK;
	}


} //end of namespace LinBox
#endif

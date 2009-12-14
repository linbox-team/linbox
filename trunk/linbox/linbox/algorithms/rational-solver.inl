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
#include <linbox/blackbox/lambda-sparse.h>
#include <linbox/blackbox/transpose.h>
#include <linbox/blackbox/diagonal.h>
#include <linbox/blackbox/compose.h>
#include <linbox/algorithms/lifting-container.h>
#include <linbox/algorithms/rational-reconstruction.h>
#include <linbox/algorithms/matrix-inverse.h>
#include <linbox/algorithms/matrix-hom.h>
#include <linbox/algorithms/blackbox-container.h>
#include <linbox/algorithms/massey-domain.h>
#include <linbox/algorithms/blackbox-block-container.h>
#include <linbox/algorithms/block-massey-domain.h>
#include <linbox/algorithms/vector-fraction.h>
#include <linbox/ffpack/ffpack.h>
#include <linbox/fflas/fflas.h>
#include <linbox/solutions/methods.h>
#include <linbox/util/debug.h>
#include <linbox/linbox-config.h>
#include <linbox/field/multimod-field.h>
#include <linbox/blackbox/block-hankel-inverse.h>



#ifdef __LINBOX_BLAS_AVAILABLE
#include <linbox/config-blas.h>
#include <linbox/blackbox/blas-blackbox.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/algorithms/blas-domain.h>
#include <linbox/matrix/factorized-matrix.h>
#include <linbox/util/timer.h>
#endif

//#define DEBUG_DIXON 
//#define DEBUG_INC
//#define SKIP_NONSINGULAR

namespace LinBox {
	
	template <class Prime>
	inline bool checkBlasPrime(const Prime p) {
		return p < Prime(67108863);
	}

	template<>
	inline bool checkBlasPrime(const std::vector<integer> p){
		bool tmp=true;
		for (size_t i=0;i<p.size();++i)
			if  (p[i] >= integer(67108863)) {tmp=false;break;}
		
		return tmp;
	}


	// SPECIALIZATION FOR WIEDEMANN 	

	// note: if Vector1 != Vector2 compilation of solve or solveSingluar will fail (via an invalid call to sparseprecondition)!
	// maybe they should not be templated separately, or sparseprecondition should be rewritten

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,WiedemannTraits>::solve (Vector1& num, Integer& den,
											  const IMatrix& A,
											  const Vector2& b,
											  const bool old=false,
											  int maxPrimes) const {
		SolverReturnStatus status=SS_FAILED;

		switch (A.rowdim() == A.coldim() ? solveNonsingular(num, den, A,b) : SS_SINGULAR) {

		case SS_OK:
			status=SS_OK;
			break;

		case SS_SINGULAR:
			std::cerr<<"switching to singular\n";
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
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,WiedemannTraits>::solveNonsingular (Vector1& num, Integer& den,
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
#ifdef RSTIMING
			tNonsingularSetup.clear();
			tNonsingularSetup.start();
#endif			
			_prime = prime;
			if (F != NULL) delete F;
			F=new Field(prime);				
			MatrixHom::map (Ap, A, *F);
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

			//CSRSparseMatrix<Field> csr_Ap(*F,*Ap);

			//typedef CSRSparseMatrix<Field> FMatrix;
			typedef SparseMatrix<Field> FMatrix;

			typedef WiedemannLiftingContainer<Ring, Field, IMatrix, FMatrix, FPolynomial> LiftingContainer;
			
			LiftingContainer lc(_R, *F, A, *Ap, MinPoly, b,_prime);
			
			RationalReconstruction<LiftingContainer> re(lc);
			
			re.getRational(num, den, 0);
#ifdef RSTIMING
			ttNonsingularSolve.update(re, lc);
#endif
			return SS_OK;
		}
	}       

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,WiedemannTraits>::solveSingular (Vector1& num, Integer& den,
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
			MatrixHom::map (Ap, A, *F);//std::cerr<<"after\n";
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
			
			re.getRational(num, den, 0); 


			if (Q    != NULL) {

				/*
				  typename Ring::Element lden;
				  _R. init (lden, 1);
				  typename Vector1::iterator p;		
				  for (p = answer.begin(); p != answer.end(); ++ p)
				  _R. lcm (lden, lden, p->second);

				*/

				IVector Qx(num.size());

				/*
				  typename IVector::iterator p_x;
						
				  for (p = answer.begin(), p_x = x. begin(); p != answer.end(); ++ p, ++ p_x) {					
				  _R. mul (*p_x, p->first, lden);					
				  _R. divin (*p_x, p->second);					
				  }
				*/

				Q->apply(Qx, num);
				/*
				  for (p=answer.begin(),p_x=Qx.begin(); p != answer.end();++p,++p_x){
				  p->first=*p_x;
				  p->second=lden;
				  }					
				*/
				num = Qx;
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
											  

	/*
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
	*/
	
	
	// SPECIALIZATION FOR BLOCK WIEDEMANN 	

	// note: if Vector1 != Vector2 compilation of solve or solveSingluar will fail (via an invalid call to sparseprecondition)!
	// maybe they should not be templated separately, or sparseprecondition should be rewritten

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,BlockWiedemannTraits>::solve (Vector1& num, Integer& den,
											       const IMatrix& A,
											       const Vector2& b,
											       const bool old=false,
											       int maxPrimes) const {
		SolverReturnStatus status=SS_FAILED;
		
		switch (A.rowdim() == A.coldim() ? solveNonsingular(num, den, A,b) : SS_SINGULAR) {
			
		case SS_OK:
			status=SS_OK;
			break;
			
		case SS_SINGULAR:
			std::cerr<<"switching to singular\n";
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
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,BlockWiedemannTraits>::solveNonsingular (Vector1& num, Integer& den,
													  const IMatrix& A,
													  const Vector2& b,
													  int maxPrimes) const {
		// checking if matrix is square
		linbox_check(A.rowdim() == A.coldim());
	
		// checking size of system
		linbox_check(A.rowdim() == b.size());

		size_t m,n;
		integer tmp,tmproot;
		tmp=A.coldim();
		//m = n = tmp.bitsize();
		//m = n = sqrt(tmp);
		//m = n = root(tmp,3); // wrong # args to root. -bds
		m = n = root(tmproot, tmp,3);
		m = n = tmproot;
		std::cout<<"block factor= "<<m<<"\n";;
		typedef SparseMatrix<Field> FMatrix;		
		FMatrix *Ap;

		Field F(_prime);
		MatrixHom::map (Ap, A, F);
		Transpose<FMatrix > Bp(*Ap);
		std::cout<<"Ap:\n";
		Ap->write(std::cout);
		typedef BlockWiedemannLiftingContainer<Ring, Field, Transpose<IMatrix >, Transpose<FMatrix > > LiftingContainer;
		
		Transpose<IMatrix> B(A);
		
		LiftingContainer lc(_R, F, B, Bp, b,_prime, m, n);
	
		RationalReconstruction<LiftingContainer> re(lc);
		
		re.getRational(num, den, 0);
#ifdef RSTIMING		
		ttNonsingularSolve.update(re, lc);
#endif
		delete Ap;		
		
		return SS_OK;	
	}       

	// END OF SPECIALIZATION FOR BLOCK WIEDEMANN

		   

	// SPECIALIZATION FOR DIXON
	
	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,DixonTraits>::solve 
	(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, const bool old, int maxP, const SolverLevel level ) const {

		SolverReturnStatus status;	
		int maxPrimes=maxP;
		while (maxPrimes > 0){
#ifdef SKIP_NONSINGULAR
			switch (SS_SINGULAR) {
#else
				switch (A.rowdim() == A.coldim() ? solveNonsingular(num, den,A,b,old,maxPrimes) : SS_SINGULAR) {
#endif
					
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
					status = solveSingular(num, den,A,b,maxPrimes,level);
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
	(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, bool oldMatrix, int maxPrimes) const {

	
#ifdef DEBUG_DIXON
		std::cout << "entering nonsingular solver\n";
#endif
		int trials = 0, notfr;

		// history sensitive data for optimal reason
		static const IMatrix* IMP;
		
		static BlasBlackbox<Field>* FMP = NULL;
		static Field *F=NULL;

		do {

			//if (trials == maxPrimes) return SS_SINGULAR;			
			//if (trials != 0) chooseNewPrime();
			//trials++;
#ifdef DEBUG_DIXON
			//std::cout << "_prime: "<<_prime<<"\n";
			std::cout<<"A:=\n";
			A.write(std::cout);
			std::cout<<"b:=\n";
			for (size_t i=0;i<b.size();++i) std::cout<<b[i]<<" , ";
			std::cout<<std::endl;			
#endif		       
#ifdef RSTIMING
			tNonsingularSetup.start();
#endif
			typedef typename Field::Element Element;
			typedef typename Ring::Element Integer;

			// checking size of system
			linbox_check(A.rowdim() == A.coldim());
			linbox_check(A.rowdim() == b.size());

			LinBox::integer tmp;		
		
			// if input matrix A is different one.
			if (!oldMatrix) {
				if (trials == maxPrimes) return SS_SINGULAR;
                                if (trials != 0) chooseNewPrime();
                                trials++;

				//delete IMP;
		
				// Could delete a non allocated matrix -> segfault
				//delete FMP;
		
				IMP = &A;					
		
				if (F != NULL) delete F;
				F= new Field (_prime);					
		
				//FMP = new BlasBlackbox<Field>(*F, A.rowdim(),A.coldim());
			
				MatrixHom::map (FMP, A, *F); // use MatrixHom to reduce matrix PG 2005-06-16
				//typename BlasBlackbox<Field>::RawIterator iter_p  = FMP->rawBegin();
				//typename IMatrix::ConstRawIterator iter  = A.rawBegin();
				//for (;iter != A.rawEnd();++iter,++iter_p)
				//	F->init(*iter_p, _R.convert(tmp,*iter));

#ifdef DEBUG_DIXON
				std::cout<< "p = ";
				F->write(std::cout);
				std::cout<<" A mod p :=\n";
				FMP->write(std::cout);
#endif				
			
				if (!checkBlasPrime(_prime)){
					if (FMP != NULL) delete FMP;
					FMP = new BlasBlackbox<Field>(*F, A.rowdim(),A.coldim());
					notfr = MatrixInverse::matrixInverseIn(*F,*FMP);
				}
				else {
					BlasBlackbox<Field> *invA = new BlasBlackbox<Field>(*F, A.rowdim(),A.coldim());
					BlasMatrixDomain<Field> BMDF(*F);
#ifdef RSTIMING
					tNonsingularSetup.stop();
					ttNonsingularSetup += tNonsingularSetup;
					tNonsingularInv.start();
#endif				
					BMDF.invin(*invA, *FMP, notfr); //notfr <- nullity
					if (FMP != NULL) delete FMP;
					FMP = invA;
					//std::cout << "notfr = " << notfr << std::endl;
					//std::cout << "inverse mod p: " << std::endl;
					//FMP->write(std::cout, *F);
#ifdef RSTIMING
					tNonsingularInv.stop();
					ttNonsingularInv += tNonsingularInv;
#endif
				}
			}
			else {
#ifdef RSTIMING
				tNonsingularSetup.stop();
				ttNonsingularSetup += tNonsingularSetup;
#endif
				notfr = 0;
			}
		} while (notfr);

#ifdef DEBUG_DIXON
		std::cout<<"A^-1 mod p :=\n";
		FMP->write(std::cout);
#endif		

		typedef DixonLiftingContainer<Ring,Field,IMatrix,BlasBlackbox<Field> > LiftingContainer;		
		LiftingContainer lc(_R, *F, A, *FMP, b, _prime);
		RationalReconstruction<LiftingContainer > re(lc);
		if (!re.getRational(num, den, 0)){
			//delete FMP;
			return SS_FAILED;
		}
#ifdef RSTIMING
		ttNonsingularSolve.update(re, lc);
#endif	
		//delete FMP;
		//if (F!=NULL)
		//	delete F;
		return SS_OK;
	}

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,DixonTraits>::solveSingular 
	(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, int maxPrimes, const SolverLevel level) const {	
		return monolithicSolve (num, den, A, b, false, false, maxPrimes, level);
	}

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,DixonTraits>::findRandomSolution 
	(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, int maxPrimes, const SolverLevel level ) const {
		
		return monolithicSolve (num, den, A, b, false, true, maxPrimes, level);
	}


	// Most solving is done by the routine below.
	// There used to be one for random and one for deterministic, but they have been merged to ease with 
	//  repeated code (certifying inconsistency, optimization are 2 examples)

	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,DixonTraits>::monolithicSolve 
	(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, bool makeMinDenomCert, bool randomSolution,
	 int maxPrimes, const SolverLevel level) const {

		if (level == SL_MONTECARLO && maxPrimes > 1) 
			std::cout << "WARNING: Even if maxPrimes > 1, SL_MONTECARLO uses just one prime." << std::endl;
		//if (makeMinDenomCert && !randomSolution) 
		//	std::cout << "WARNING: Will not compute a certificate of minimal denominator deterministically." << std::endl;
		if (makeMinDenomCert && level == SL_MONTECARLO) 
			std::cout << "WARNING: No certificate of min-denominality generated due to  level=SL_MONTECARLO" << std::endl;
		int trials = 0;
		while (trials < maxPrimes){
			if (trials != 0) chooseNewPrime();
			trials++;
#ifdef DEBUG_DIXON
			std::cout << "_prime: "<<_prime<<"\n";
#endif		       
#ifdef RSTIMING
			tSetup.start();
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

			BlasBlackbox<Ring> A_check(_R, A); // used to check answer later
		
			// TAS_xxx stands for Transpose Augmented System (A|b)t
			// this provides a factorization (A|b) = TAS_Pt . TAS_Ut . TAS_Qt . TAS_Lt
			// such that    
			// - TAS_P . (A|b) . TAS_Q   has nonzero principal minors up to TAS_rank
			// - TAS_Q permutes b to the (TAS_rank)th column of A iff the system is inconsistent mod p

			BlasBlackbox<Field>* TAS_factors = new BlasBlackbox<Field>(F, A.coldim()+1, A.rowdim());
			Hom<Ring, Field> Hmap(_R, F);

			BlasBlackbox<Field> *Ap;
			MatrixHom::map(Ap, A, F);

			for (size_t i=0;i<A.rowdim();++i)
				for (size_t j=0;j<A.coldim();++j)
					TAS_factors->setEntry(j,i, Ap->getEntry(i,j));
					
			delete Ap;

			for (size_t i=0;i<A.rowdim();++i){
				typename Field::Element tmpe;
				F.init(tmpe);
				F.init(tmpe,_R.convert(tmp,b[i]));
				TAS_factors->setEntry(A.coldim(),i, tmpe);
			}
#ifdef RSTIMING
			tSetup.stop();
			ttSetup += tSetup;
			tLQUP.start();
#endif
			LQUPMatrix<Field>* TAS_LQUP = new LQUPMatrix<Field>(F, *TAS_factors);
			size_t TAS_rank = TAS_LQUP->getrank();
						
			// check consistency. note, getQ returns Qt.
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
			std::cout << "P takes (0 1 ...) to (";
			for (size_t i=0; i<A.rowdim(); i++) std::cout << srcRow[i] << ' '; std::cout << ')' << std::endl;
			std::cout << "Q takes (0 1 ...) to (";
			for (size_t i=0; i<A.coldim()+1; i++) std::cout << srcCol[i] << ' '; std::cout << ')' << std::endl;
#endif
			
			bool appearsInconsistent = (srcCol[TAS_rank-1] == A.coldim());
			size_t rank = TAS_rank - (appearsInconsistent ? 1 : 0);
#ifdef DIXON_DEBUG
			std::cout << "TAS_rank, rank: " << TAS_rank << ' ' << rank << std::endl;
#endif
#ifdef RSTIMING
			tLQUP.stop();
			ttLQUP += tLQUP;
#endif
			if (rank == 0) { 
				delete TAS_LQUP;
				delete TAS_factors;
				//special case when A = 0, mod p. dealt with to avoid crash later
				bool aEmpty = true;
				if (level >= SL_LASVEGAS) { // in monte carlo, we assume A is actually empty
					typename BlasBlackbox<Ring>::RawIterator iter = A_check.rawBegin();
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
					/*
					// both A and b are all zero.
					for (size_t i=0; i<answer.size(); i++) {
					answer[i].first = _rzero;
					answer[i].second = _rone;
					}
					*/
					_R. assign (den, _rone);
					for (typename Vector1::iterator p = num. begin(); p != num. end(); ++ p)
						_R. assign (*p, _rzero);

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
	
			BlasBlackbox<Field>* Atp_minor_inv = NULL;

			if ((appearsInconsistent && level > SL_MONTECARLO) || randomSolution == false) {
				// take advantage of the (LQUP)t factorization to compute
				// an inverse to the leading minor of (TAS_P . (A|b) . TAS_Q) 
#ifdef RSTIMING
				tFastInvert.start();
#endif
				Atp_minor_inv = new BlasBlackbox<Field>(F, rank, rank);
				

				FFPACK::LQUPtoInverseOfFullRankMinor(F, rank, TAS_factors->getPointer(), A.rowdim(), 
								     TAS_Qt.getPointer(), 
								     Atp_minor_inv->getPointer(), rank);
#ifdef RSTIMING
				tFastInvert.stop();
				ttFastInvert += tFastInvert;
#endif
			}

			delete TAS_LQUP;
			delete TAS_factors;

			if (appearsInconsistent && level <= SL_MONTECARLO)
				return SS_INCONSISTENT;

			if (appearsInconsistent) {
#ifdef RSTIMING
				tCheckConsistency.start();
#endif			
				std::vector<Integer> zt(rank);
				for (size_t i=0; i<rank; i++)
					_R.assign(zt[i], A.getEntry(srcRow[rank], srcCol[i]));

				BlasBlackbox<Ring> At_minor(_R, rank, rank);
				for (size_t i=0; i<rank; i++)
					for (size_t j=0; j<rank; j++)
						_R.assign(At_minor.refEntry(j, i), A.getEntry(srcRow[i], srcCol[j]));
#ifdef DEBUG_INC
				At_minor.write(std::cout << "At_minor:" << std::endl);//, _R);
				Atp_minor_inv->write(std::cout << "Atp_minor_inv:" << std::endl);//, F);
				std::cout << "zt: "; for (size_t i=0; i<rank; i++) std::cout << zt[i] <<' '; std::cout << std::endl;
#endif
#ifdef RSTIMING
				tCheckConsistency.stop();
				ttCheckConsistency += tCheckConsistency;
#endif			

				LiftingContainer lc(_R, F, At_minor, *Atp_minor_inv, zt, _prime);
				
				RationalReconstruction<LiftingContainer > re(lc);
				
				Vector1 short_num(rank); Integer short_den;
				
				if (!re.getRational(short_num, short_den,0)) 
					return SS_FAILED;    // dirty, but should not be called
				// under normal circumstances
#ifdef RSTIMING
				ttConsistencySolve.update(re, lc);
				tCheckConsistency.start();
#endif
				VectorFraction<Ring> cert(_R, short_num. size());
				cert. numer = short_num;
				cert. denom = short_den;
				cert.numer.resize(b.size());
				_R.subin(cert.numer[rank], cert.denom);
				_R.init(cert.denom, 1);
				BMDI.mulin_left(cert.numer, TAS_P);
#ifdef DEBUG_INC
				cert.write(std::cout << "cert:") << std::endl;
#endif

				bool certifies = true; //check certificate
				std::vector<Integer> certnumer_A(A.coldim());
				BAR.applyVTrans(certnumer_A, A_check, cert.numer);
				typename std::vector<Integer>::iterator cai = certnumer_A.begin();
				for (size_t i=0; certifies && i<A.coldim(); i++, cai++) 
					certifies &= _R.isZero(*cai);				
#ifdef RSTIMING
				tCheckConsistency.stop();
				ttCheckConsistency += tCheckConsistency;
#endif
				if (certifies) { 
					if (level == SL_CERTIFIED) lastCertificate.copy(cert);
					return SS_INCONSISTENT;
				}
				std::cout<<"system is suspected to be inconsistent but it was only a bad prime\n";
				continue; // try new prime. analogous to u.A12 != A22 in Muld.+Storj.
			}
			
#ifdef RSTIMING
			tMakeConditioner.start();
#endif
			// we now know system is consistent mod p.
			BlasBlackbox<Ring> A_minor(_R, rank, rank);    // -- will have the full rank minor of A
			BlasBlackbox<Field> *Ap_minor_inv;          // -- will have inverse mod p of A_minor
			BlasBlackbox<Ring> *P = NULL, *B = NULL;   // -- only used in random case 

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
#ifdef RSTIMING
				tMakeConditioner.stop();
				ttMakeConditioner += tMakeConditioner;
#endif
			
				if (makeMinDenomCert && level >= SL_LASVEGAS){
					B = new BlasBlackbox<Ring>(_R, rank, A.coldim());
					for (size_t i=0; i<rank; i++)
						for (size_t j=0; j<A.coldim(); j++)
							_R.assign(B->refEntry(i, j), A_check.getEntry(srcRow[i],j));
				}					
			}
			else {
				P = new BlasBlackbox<Ring>(_R, A.coldim(), rank);	
				B = new BlasBlackbox<Ring>(_R, rank,A.coldim());
				BlasBlackbox<Field> Ap_minor(F, rank, rank);
				Ap_minor_inv = new BlasBlackbox<Field>(F, rank, rank);
				int nullity;
				
				LinBox::integer tmp=0;
				size_t maxBitSize = 0;
				for (size_t i=0; i<rank; i++)
					for (size_t j=0; j<A.coldim(); j++){
						_R.assign(B->refEntry(i, j), A_check.getEntry(srcRow[i], j));
						_R.convert(tmp, A_check.getEntry(srcRow[i], j));
						maxBitSize = std::max(maxBitSize, tmp.bitsize());
					}
#ifdef RSTIMING
				bool firstLoop = true;
#endif
				// prepare B to be preconditionned through BLAS matrix mul
				MatrixApplyDomain<Ring, BlasBlackbox<Ring> > MAD(_R,*B);
				MAD.setup(2);

				do { // O(1) loops of this preconditioner expected
#ifdef RSTIMING
					if (firstLoop) 
						firstLoop = false;
					else
						tMakeConditioner.start();
#endif
					// compute P a n*r random matrix of entry in [0,1]
					typename BlasBlackbox<Ring>::RawIterator iter;
					for (iter = P->rawBegin(); iter != P->rawEnd(); ++iter) {
						if (rand() > RAND_MAX/2) 
							_R.assign(*iter, _rone);
						else
							_R.assign(*iter, _rzero);
					}

					// compute A_minor = B.P
					/*
					  if (maxBitSize * log((double)A.coldim()) > 53) 
					  MD.mul(A_minor, *B, *P);
					  else {
					  double *B_dbl= new double[rank*A.coldim()]; 
					  double *P_dbl= new double[A.coldim()*rank]; 
					  double *A_minor_dbl = new double[rank*rank];
					  for (size_t i=0;i<rank;++i)
					  for (size_t j=0;j<A.coldim(); j++){
					  _R.convert(B_dbl[j+i*A.coldim()], B->getEntry(i,j));
					  _R.convert(P_dbl[i+j*rank], P->getEntry(j,i));
					  }
					  cblas_dgemm(CblasRowMajor, CblasNoTrans,
					  CblasNoTrans,
					  rank, rank, A.coldim(), 1,
					  B_dbl, A.coldim(), P_dbl, rank, 0,A_minor_dbl, rank);

					  for (size_t i=0;i<rank;++i)
					  for (size_t j=0;j<rank;++j)
					  _R.init(A_minor.refEntry(i,j),A_minor_dbl[j+i*rank]);

					  delete[] B_dbl;	
					  delete[] P_dbl;
					  delete[] A_minor_dbl;
					  }
					*/
								
					MAD.applyM(A_minor,*P);



					// set Ap_minor = A_minor mod p, try to compute inverse
					for (size_t i=0;i<rank;++i)
						for (size_t j=0;j<rank;++j)
							F.init(Ap_minor.refEntry(i,j),
							       _R.convert(tmp,A_minor.getEntry(i,j)));
#ifdef RSTIMING
					tMakeConditioner.stop();
					ttMakeConditioner += tMakeConditioner;
					tInvertBP.start();
#endif
					BMDF.inv(*Ap_minor_inv, Ap_minor, nullity); 
#ifdef RSTIMING
					tInvertBP.stop();
					ttInvertBP += tInvertBP;
#endif
				} while (nullity > 0);
			}
			// Compute newb = (TAS_P.b)[0..(rank-1)]
			std::vector<Integer> newb(b);
			BMDI.mulin_right(TAS_P, newb);
			newb.resize(rank);

			BlasBlackbox<Ring>  BBA_minor(_R,A_minor);
			//BlasBlackbox<Field> BBA_inv(F,*Ap_minor_inv);
			//BlasMatrix<Integer>  BBA_minor(A_minor);
			//BlasMatrix<Element> BBA_inv(*Ap_minor_inv);
			//LiftingContainer lc(_R, F, BBA_minor, BBA_inv, newb, _prime);
			LiftingContainer lc(_R, F, BBA_minor, *Ap_minor_inv, newb, _prime);

#ifdef DEBUG_DIXON
			std::cout<<"length of lifting: "<<lc.length()<<std::endl;
#endif
			RationalReconstruction<LiftingContainer > re(lc);

			Vector1 short_num(rank); Integer short_den;
			
			if (!re.getRational(short_num, short_den,0))
				return SS_FAILED;    // dirty, but should not be called
			// under normal circumstances		
#ifdef RSTIMING
			ttSystemSolve.update(re, lc);
			tCheckAnswer.start();
#endif		
			VectorFraction<Ring> answer_to_vf(_R, short_num. size());
			answer_to_vf. numer = short_num;
			answer_to_vf. denom = short_den;

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
				//BAR.applyVspecial(newNumer, *P, answer_to_vf.numer);
				
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
#ifdef RSTIMING
					tCheckAnswer.stop();
					ttCheckAnswer += tCheckAnswer;
#endif
					continue; //go to start of main loop
				}
			}
			
			//answer_to_vf.toFVector(answer);
			num = answer_to_vf. numer;
			den = answer_to_vf. denom;
#ifdef RSTIMING  
			tCheckAnswer.stop();
			ttCheckAnswer += tCheckAnswer;
#endif			
			if (makeMinDenomCert && level >= SL_LASVEGAS){  // && randomSolution) {
				// To make this certificate we solve with the same matrix as to get the 
				// solution, except transposed.
#ifdef RSTIMING
				tCertSetup.start();
#endif
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
				// B     | *B (== TAS_P . A,  but only top #rank rows)
				// c     | newb (== TAS_P . b,   but only top #rank rows)
				// P     | P
				// q     | q
				// U     | {0, 1}
				// u     | u
				// z-hat | lastCertificate
				
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
#ifdef RSTIMING
				tCertSetup.stop();
				ttCertSetup += tCertSetup;
#endif
				//LiftingContainer lc2(_R, F, BBA_minor, BBA_inv, q, _prime);
				LiftingContainer lc2(_R, F, A_minor, *Ap_minor_inv, q, _prime);

				RationalReconstruction<LiftingContainer> re(lc2);
				Vector1 u_num(rank); Integer u_den;
				if (!re.getRational(u_num, u_den,0)) return SS_FAILED;

#ifdef RSTIMING
				ttCertSolve.update(re, lc2);
				tCertMaking.start();
#endif
				// remainder of code does   z <- denom(partial_cert . Mr) * partial_cert * Qt 
				VectorFraction<Ring> u_to_vf(_R, u_num.size());
				u_to_vf. numer = u_num;
				u_to_vf. denom = u_den;
				std::vector<Integer> uB(A.coldim());
				BAR.applyVTrans(uB, *B, u_to_vf.numer);

				// 				std::cout << "BP: ";
				// 				A_minor.write(std::cout, _R) << std::endl;
				// 				std::cout << "q: ";
				// 				for (size_t i=0; i<rank; i++) std::cout << q[i]; std::cout << std::endl;
				// 				u_to_vf.write(std::cout  << "u: ") << std::endl;

				Integer numergcd = _rzero;
				vectorGcdIn(numergcd, _R, uB);

				// denom(partial_cert . Mr) = partial_cert_to_vf.denom / numergcd
				VectorFraction<Ring> z(_R, b.size()); //new constructor
				u_to_vf.numer.resize(A.rowdim());

				BMDI.mul(z.numer, u_to_vf.numer, TAS_P);

				z.denom = numergcd;

				// 				z.write(std::cout << "z: ") << std::endl;

				if (level >= SL_CERTIFIED)
					lastCertificate.copy(z);
				
				// output new certified denom factor
				Integer znumer_b, zbgcd;
				VDR.dotprod(znumer_b, z.numer, b);
				_R.gcd(zbgcd, znumer_b, z.denom);
				_R.div(lastCertifiedDenFactor, z.denom, zbgcd);

				if (level >= SL_CERTIFIED) 
					_R.div(lastZBNumer, znumer_b, zbgcd);
#ifdef RSTIMING
				tCertMaking.stop();
				ttCertMaking += tCertMaking;
#endif
			}

			delete Ap_minor_inv;			
			delete B;
			
			if (randomSolution) {delete P;}

			// done making certificate, lets blow this popstand
			return SS_OK;
		}
		return SS_FAILED; //all primes were bad
	}




	/*
	 * Specialization for Block Hankel method
	 */
	
	// solve non singular system using block Hankel
	template <class Ring, class Field, class RandomPrime>
	template <class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus RationalSolver<Ring,Field,RandomPrime,BlockHankelTraits>::solveNonsingular 
	(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, size_t blocksize, int maxPrimes) const {

		linbox_check(A.rowdim() == A.coldim());
		linbox_check(A.rowdim() % blocksize == 0);
		
		typedef typename Field::Element Element;

		

		// reduce the matrix mod p
		Field F(_prime);
		typedef typename IMatrix::template rebind<Field>::other FMatrix;
		FMatrix *Ap;
		typename IMatrix::template rebind<Field>()( Ap, A, F);

		// precondition Ap  with a random diagonal Matrix
		typename Field::RandIter G(F,0,123456);
		std::vector<Element> diag(Ap->rowdim());
		
		for(size_t i=0;i<Ap->rowdim();++i){
			do {
				G.random(diag[i]);
			} while(F.isZero(diag[i]));
		}
	


		Diagonal<Field> D(F, diag);
		
		Compose<Diagonal<Field>, FMatrix> DAp(D,*Ap);

		size_t n = A.coldim();
		size_t numblock = n/blocksize;

		// generate randomly U and V
		BlasMatrix<Element> U(blocksize,A.rowdim()), V(A.coldim(),blocksize);
	
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
		std::vector<Element> y(n), x(n, 1);

#ifdef RSTIMING
		chrono.stop();
		std::cout<<"inverse block hankel: "<<chrono<<"\n";
#endif

		typedef BlockHankelLiftingContainer<Ring,Field,IMatrix,Compose<Diagonal<Field>,FMatrix>, BlasMatrix<Element> > LiftingContainer;		
		LiftingContainer lc(_R, F, A, DAp, D, Hinv, U, V, b, _prime);		
		RationalReconstruction<LiftingContainer > re(lc);
		
		if (!re.getRational(num, den, 0)) return SS_FAILED;

#ifdef RSTIMING
		std::cout<<"lifting bound computation : "<<lc.ttSetup<<"\n";
		std::cout<<"residue computation       : "<<lc.ttRingApply<<"\n";
		std::cout<<"rational reconstruction   : "<<re.ttRecon<<"\n";
#endif
		
		return SS_OK;
	}





} //end of namespace LinBox




/* Author Z. Wan
 * Modified to fit in linbox
 * Implementation the algorithm in manuscript, available at http://www.cis.udel.edu/~wan/jsc_wan.ps
 */



#ifndef __RATIONAL_SOLVER2__H__
#define __RATIONAL_SOLVER2__H__

#include <memory.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <linbox/integer.h>
#include <linbox/algorithms/rational-reconstruction2.h>

namespace LinBox {

#if __LINBOX_HAVE_DGETRF && __LINBOX_HAVE_DGETRI
	template <class Ring, class Field, class RandomPrime>
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::cblas_dgeinv(double* M, int n) {
		enum CBLAS_ORDER order = CblasRowMajor;
		int lda = n;
		int P[n];
		int ierr = clapack_dgetrf (order, n, n, M, lda, P);
		if (ierr != 0) {
                    std::cerr << "In RationalSolver::cblas_dgeinv Matrix is not full rank" << std::endl;
                    return -1;
		}
		clapack_dgetri (order, n, M, lda, P);
		return 0;
	}

	template <class Ring, class Field, class RandomPrime>
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::cblas_rsol (int n, const double* M, integer* numx, integer& denx, double* b) {
		if (n < 1) return 0;
		double* IM = new double[n * n];
		memcpy ((void*)IM, (const void*)M, sizeof(double)*n*n);
		int ret;
		//compute the inverse by flops
		ret = cblas_dgeinv (IM, n); 
		if (ret != 0) {delete[] IM; return 1;}

		double mnorm = cblas_dOOnorm(M, n, n);
		// residual
		double* r = new double [n];
		// A^{-1}r
		double* x = new double [n];
		//ax  = A x
		double* ax = new double [n];
		// a digit, d \approx \alpha x
		double* d = new double [n];

		const double* p2;
		double* pd;
		const double T = 1 << 30;
	
		integer* num = new integer [n];
		integer* p_mpz;
		integer tmp_mpz, den, denB, B;

		den = 1;
		// compute the hadamard bound
		cblas_hbound (denB, n, n, M);      
		B = denB * denB;
		// shouble be a check for tmp_mpz
		tmp_mpz = 2 * mnorm + cblas_dmax (n, b, 1);
		B <<= 1; B *= tmp_mpz; //B *= tmp_mpz;
	
		//double log2 = log (2.0);
		double log2 = M_LN2;
		// r = b
		memcpy ((void*) r, (const void*) b, sizeof(double)*n);

		do  {
			cblas_dapply (n, n, IM, r, x);
			// compute ax
			cblas_dapply (n, n, M, x, ax);
			// compute ax = ax -r, the negative of residual
			cblas_daxpy (n, -1, r, 1, ax, 1);
			// compute possible shift
			double normr1, normr2, normr3, shift1, shift2;
			normr1 = cblas_dmax(n, r, 1);
			normr2 = cblas_dmax(n, ax, 1);
			normr3 = cblas_dmax(n, x, 1);
			//try to find a good scalar
			int shift = 30;
			if (normr2 <.0000000001) 
				shift = 30;
			else {
				shift1 = floor(log (normr1 / normr2) / log2) - 2;
				shift = (int)(30 < shift1 ? 30 : shift1);
			}

			normr3 = normr3 > 2 ? normr3 : 2;
			shift2 = floor(53. * log2 / log (normr3));
			shift = (int)(shift < shift2 ? shift : shift2);

			if (shift <= 0) {
#ifdef DEBUGRC
				printf ("%s", "Bad scalar \n");
				printf("%f, %f\n", normr1, normr2);
				printf ("%d, shift = ", shift);
				printf ("OO-norm of matrix: %f\n", cblas_dOOnorm(M, n, n));
				printf ("OO-norm of inverse: %f\n", cblas_dOOnorm(IM, n, n));
				printf ("Error, abort\n"); 
#endif
				delete[] IM; delete[] r; delete[] x; delete[] ax; delete[] d; delete[] num; 
				return 2;
			}

			int scalar = ((long long int)1 << shift);
			for (pd = d, p2 = x; pd != d + n; ++ pd, ++ p2) 
				//better use round, but sun sparc machine doesnot supprot it
				*pd = floor (*p2 * scalar);

			// update den
			den <<= shift; 
			//update num
			update_num (num, n, d, shift);
		
#ifdef DEBUGRC
			printf ("in iteration\n");
			printf ("residual=\n");
			printvec (r,  n);
			printf ("A^(-1) r\n");
			printvec (x,  n);
			printf ("scalar= ");
			printf ("%d \n", scalar);
			printf ("One digit=\n");
			printvec (d, n);
			printf ("Current bound= \n");
			std::cout << B;
			printf ("den= \n");
			std::cout << den;
			printf ("accumulate numerator=\n");
			printvec (num, n);
#endif
			// update r = r * shift - M d
			double tmp = 2 * mnorm + cblas_dmax (n, r, 1);
			if (tmp < T) update_r_int (r, n, M, d, shift);
			else update_r_ll (r, n, M, d, shift);
			//update_r_ll (r, n, M, d, shift);
		} while (den < B);
	
		integer q, rem, den_lcm, tmp_den;
		integer* p_x, * p_x1;
		p_mpz = num;
		p_x = numx;
		// construct first answer
		rational_reconstruction (*p_x, denx, *p_mpz, den, denB);
		++ p_mpz;
		++ p_x;

		int sgn;
		for (; p_mpz != num + n; ++ p_mpz, ++ p_x)  {
			sgn = sign (*p_mpz);
			tmp_mpz = denx * (*p_mpz);
			tmp_mpz = abs (tmp_mpz);
			integer::divmod (q, rem, tmp_mpz, den);
		
			if ( rem < denx)  {
				if (sgn >= 0)
					*p_x = q;
				else 
					*p_x = -q;
			}
			else {
				rem = den - rem;
				q += 1;
				if (rem < denx) {
					if (sgn >= 0)
						*p_x = q;
					else
						*p_x = -q;
				}
				else {
					rational_reconstruction (*p_x, tmp_den, *p_mpz, den, denB);
					lcm (den_lcm, tmp_den, denx);
					integer::divexact (tmp_mpz, den_lcm, tmp_den);
					integer::mul (*p_x, *p_x, tmp_mpz);
					integer::divexact (tmp_mpz, den_lcm, denx);
					denx = den_lcm;
					for (p_x1 = numx; p_x1 != p_x; ++ p_x1)
						integer::mul (*p_x1, *p_x1, tmp_mpz);
				}
			}	
		}
#ifdef DEBUGRC
		std::cout << "raiotanl answer\nCommon den = ";
		std::cout << denx;
		std::cout << "\nNumerator= \n";
		printvec (numx, n);
#endif
	
		//normalize the answer
		if (denx != 0) {	
			integer g; g = denx;
			for (p_x = numx; p_x != numx + n; ++ p_x)
				g = gcd (g, *p_x);
			for (p_x = numx; p_x != numx + n; ++ p_x)
				integer::divexact (*p_x, *p_x, g);
			integer::divexact (denx, denx, g);
		}
	
		//check if the answer is correct, not necessary
		cblas_mpzapply (n, n, M, (const integer*)numx, num);
		integer* sb = new integer [n];
		double* p;
		for (p_mpz = sb, p = b; p_mpz != sb + n; ++ p_mpz, ++ p) {
			*p_mpz = *p;
			integer::mulin(*p_mpz, denx);
		}
		ret = 0;
		for (p_mpz = sb, p_x = num; p_mpz != sb + n; ++ p_mpz, ++ p_x) 
			if (*p_mpz != *p_x) {
				ret = 3;
				break;
			}
#ifdef DEBUGRC
		if (ret == 3) {
	
			std::cout << "Input matrix:\n";
			for (int i = 0; i < n; ++ i) {
				const double* p = M + (i * n);
				printvec (p, n);
			}
			std::cout << "Input rhs:\n";
			printvec (b, n);
			std::cout << "Common den: " << denx << '\n';
			std::cout << "Numerator: ";
			printvec (numx, n);
			std::cout << "A num: ";
			printvec (num, n);
			std::cout << "denx rhs: ";
			printvec (sb, n);
		}
#endif

		// garbage collector
		delete[] IM; delete[] r; delete[] x; delete[] ax; delete[] d; delete[] num; delete[] sb;

		return ret;
	}

#endif


	template <class Ring, class Field, class RandomPrime>
	/* apply  y <- Ax */
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::cblas_dapply (int m, int n, const double* A, const double* x, double* y) {
		cblas_dgemv (CblasRowMajor, CblasNoTrans, m, n, 1, A, n, x, 1, 0, y, 1);
		return 0;
	}

	template <class Ring, class Field, class RandomPrime>
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::cblas_mpzapply (int m, int n, const double* A, const integer* x, integer* y) {
		const double* p_A;
		const integer* p_x;
		integer* p_y;
		integer tmp;
		for (p_A = A, p_y = y; p_y != y + m; ++ p_y) {
			*p_y = 0;
			for (p_x = x; p_x != x + n; ++ p_x, ++ p_A) {
				//mpz_set_d (tmp, *p_A);
				//mpz_addmul_si (*p_y, *p_x, (int)(*p_A));
				tmp = *p_x  * (long long int)(*p_A);
				integer::addin (*p_y, tmp);
			}
		}
		return 0;
	}

	template <class Ring, class Field, class RandomPrime>
	template <class Elt>
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::printvec (const Elt* v, int n) {
		const Elt* p;
		std::cout << "\[";
		for (p = v; p != v + n; ++ p)
			std::cout << *p << ' ';
		std::cout << "]\n";
		return 0;
	}

	template <class Ring, class Field, class RandomPrime>
	//update num, *num <- *num * 2^shift + d
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::update_num (integer* num, int n, const double* d, int shift) {
		integer* p_mpz;
		integer tmp_mpz;
		const double* pd;
		for (p_mpz = num, pd = d; p_mpz != num + n; ++ p_mpz, ++ pd) {
			(*p_mpz) = (*p_mpz) << shift;
			tmp_mpz = *pd;
			integer::add (*p_mpz, *p_mpz, tmp_mpz);
		} 
		return 0;
	}

	template <class Ring, class Field, class RandomPrime>
	//update r = r * shift - M d
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::update_r_int (double* r, int n, const double* M, const double* d, int shift) {
		int tmp;
		double* p1;
		const double* p2;
		const double* pd;
		for (p1 = r, p2 = M; p1 != r + n; ++ p1) {
			tmp = (int)(long long int) *p1;
			tmp <<= shift;
			for (pd = d; pd != d + n; ++ pd, ++ p2) {
				tmp -= (int)(long long int)*pd * (int)(long long int)*p2;
			}
			*p1 = (double)tmp;
		}
		return 0;
	}

	template <class Ring, class Field, class RandomPrime>
	//update r = r * shift - M d
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::update_r_ll (double* r, int n, const double* M, const double* d, int shift) {
		long long int tmp;
		double* p1;
		const double* p2;
		const double* pd;
		for (p1 = r, p2 = M; p1 != r + n; ++ p1) {
			tmp = (long long int) *p1;
			tmp <<= shift;
			for (pd = d; pd != d + n; ++ pd, ++ p2) {
				tmp -= (long long int)*pd * (long long int) *p2;
			}
			*p1 = tmp;
		}
		return 0;
	}

	template <class Ring, class Field, class RandomPrime>
	inline double RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::cblas_dOOnorm(const double* M, int m, int n) {
		double norm = 0;
		double old = 0;
		const double* p;
		for (p = M; p != M + (m * n); ) {
			old = norm;
			norm = cblas_dasum (n, p ,1);
			if (norm < old) norm = old;
			p += n;
		}	
		return norm;
	}

	template <class Ring, class Field, class RandomPrime>
	inline double RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::cblas_dmax (const int N, const double* a, const int inc) {
		return fabs(a[cblas_idamax (N, a, inc)]);
	}

	template <class Ring, class Field, class RandomPrime>
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::cblas_hbound (integer& b, int m, int n, const double* M){
		double norm = 0;
		const  double* p;
		integer tmp;
		b = 1;
		for (p = M; p != M + (m * n); ) {
			norm = cblas_dnrm2 (n, p ,1);
			tmp =  norm;
			integer::mulin (b, tmp);
			p += n;
		}

		return 0;
	}
}//LinBox

#endif

#endif

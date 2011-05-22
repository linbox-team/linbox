/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/algorithms/block-massey-domain.h
 * Copyright (C) 2002 Pascal Giorgi
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



#ifndef __LINBOX_massey_block_domain_H
#define __LINBOX_massey_block_domain_H

#include <vector>
#include <iostream>
#include <iomanip>

#include <linbox/util/commentator.h>
#include <linbox/util/timer.h>
#include <linbox/blackbox/dense.h>
#include <linbox/field/unparametric.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/matrix/factorized-matrix.h>
#include <linbox/algorithms/blas-domain.h>

#include <linbox/util/timer.h>

//#define  __CHECK_RESULT
//#define __DEBUG_MAPLE
//#define __CHECK_LOOP
//#define __PRINT_MINPOLY
//#define __CHECK_DISCREPANCY
//#define __CHECK_TRANSFORMATION
//#define __CHECK_SIGMA_RESULT
//#define __PRINT_SEQUENCE
#define _BM_TIMING

namespace LinBox 
{


#define DEFAULT_EARLY_TERM_THRESHOLD 20


	/** 
	    \brief Compute the linear generator of a sequence of matrices

	    * Giorgi, Jeannerod Villard algorithm from ISSAC'03
	    * This class encapsulates the functionality required for computing 
	    * the block minimal polynomial of a matrix.
	    */
	template<class _Field, class _Sequence>
	class BlockMasseyDomain {

	public:
		typedef _Field                           Field;
		typedef typename Field::Element        Element;
		typedef _Sequence                     Sequence;
		typedef BlasMatrix<Element>        Coefficient;
		
		
	private:
		Sequence                          *_container;
		Field                                      _F;
		BlasMatrixDomain<Field>                  _BMD;
		MatrixDomain<Field>                       _MD;
		unsigned long            EARLY_TERM_THRESHOLD;
		

	public:

#ifdef _BM_TIMING
		mutable Timer   ttGetMinPoly;			mutable Timer     tGetMinPoly;
		mutable Timer	ttNewDiscrepancy;		mutable Timer tNewDiscrepancy;
		mutable Timer	ttShiftSigma;			mutable Timer tShiftSigma;
		mutable Timer   ttApplyPerm;			mutable Timer   tApplyPerm; 
		mutable Timer   ttUpdateSigma;			mutable Timer tUpdateSigma;
		mutable Timer   ttInverseL;				mutable Timer tInverseL;
		mutable Timer   ttGetPermutation;		mutable Timer tGetPermutation;
		mutable Timer   ttLQUP;					mutable Timer          tLQUP;
		mutable Timer   ttDiscrepancy;			mutable Timer tDiscrepancy;
		mutable Timer   ttGetCoeff;				mutable Timer tGetCoeff;
		mutable Timer   ttCheckSequence;		mutable Timer tCheckSequence;
		mutable Timer   ttSetup;				mutable Timer tSetup;
		mutable Timer   ttMBasis;				mutable Timer tMBasis;
		mutable Timer   ttUpdateSerie;			mutable Timer tUpdateSerie;
		mutable Timer   ttBasisMultiplication;	mutable Timer tBasisMultiplication;
		mutable Timer   ttCopyingData;			mutable Timer tCopyingData;
		mutable Timer   Total; 

		void clearTimer() {
			 ttGetMinPoly.clear();     
			 ttNewDiscrepancy.clear(); 
			 ttShiftSigma.clear();     
			 ttApplyPerm.clear();      
			 ttUpdateSigma.clear();    
			 ttInverseL.clear();       
			 ttGetPermutation.clear(); 
			 ttLQUP.clear();           
			 ttDiscrepancy.clear();    
			 ttGetCoeff.clear();       
			 ttCheckSequence.clear();  
			 ttSetup.clear();   
			 ttMBasis.clear();
			 ttUpdateSerie.clear();
			 ttBasisMultiplication.clear();
			 ttCopyingData.clear(),
			 Total.clear();
		}
		
		void print(Timer& T, const  char* timer, const char* title) {
			if (&T != &Total)
				Total+=T;
			if (T.count() > 0) {
				std::cout<<title<<": "<<timer;
				for (int i=strlen(timer); i<28; i++) 
					std::cout << ' ';
				std::cout<<T<<std::endl;
			}
		}

		void printTimer() {
			print(ttSetup, "Setup", "direct");
			print(ttCheckSequence, "Rank of Seq[0]", "direct");
			print(ttGetCoeff, "Compute sequence", "direct");
			print(ttDiscrepancy, "Compute Discrepancy", "direct");
			print(ttLQUP, "LQUP","direct");
			print(ttGetPermutation, "Compute Permutation", "direct");
			print(ttApplyPerm, "Apply Permutation", "direct");
			print(ttInverseL, "Inverse of L", "direct");
			print(ttUpdateSigma, "Update Sigma", "direct");
			print(ttShiftSigma, "Shift Sigma by x", "direct");
			print(ttNewDiscrepancy, "Keep half Discrepancy", "direct");
			print(ttMBasis, "MBasis computation", "recursive");
			print(ttUpdateSerie, "Updating Power Serie", "recursive");
			print(ttBasisMultiplication, "Basis Multiplication", "recursive");
			print(ttCopyingData, "Copying Data", "recursive");
			print(Total, "Total", "");
			std::cout<<std::endl<<std::endl;
		}
#endif
	

		BlockMasseyDomain (const BlockMasseyDomain<Field, Sequence> &M, unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
			: _container(M._container), _F(M._F), _BMD(M._F), _MD(M._F),  EARLY_TERM_THRESHOLD (ett_default) {
#ifdef _BM_TIMING
			clearTimer();
#endif			
		}

		BlockMasseyDomain (Sequence *D, unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
			: _container(D), _F(D->getField ()), _BMD(D->getField ()), _MD(D->getField ()), EARLY_TERM_THRESHOLD (ett_default) {
#ifdef _BM_TIMING
			clearTimer();
#endif
		}
  
			
		// field of the domain
		const Field &getField    () const { return _F; }
		
		// sequence of the domain
		Sequence *getSequence () const { return _container; }

		// left minimal generating polynomial of the sequence
		void left_minpoly  (std::vector<Coefficient> &P) { 
			masseyblock_left(P);
		}
		
		void left_minpoly_rec  (std::vector<Coefficient> &P) { 
			masseyblock_left_rec(P);
		}


		// left minimal generating polynomial  of the sequence, keep track on degree
		void left_minpoly (std::vector<Coefficient> &phi, std::vector<size_t> &degree) {
			degree = masseyblock_left(phi);
		}		
		
		void left_minpoly_rec  (std::vector<Coefficient> &P, std::vector<size_t> &degree) { 
			degree = masseyblock_left_rec(P);
		}


		// right minimal generating polynomial of the sequence
		void right_minpoly (std::vector<Coefficient> &P) { masseyblock_right(P);}
		
		
	private:
		
		template<class Field>
		void write_maple(const Field& F, const std::vector<Coefficient> & P) {
			std::cout<<"Matrix([";
			for (size_t i=0;i< P[0].rowdim();++i){
				std::cout<<"[";
				for (size_t j=0;j< P[0].coldim();++j){
					F.write(std::cout,P[0].getEntry(i,j));
					for (size_t k=1;k<P.size();++k){
						std::cout<<"+ x^"<<k<<"*";
						F.write(std::cout,P[k].getEntry(i,j));
					}
					if (j != P[0].coldim()-1)
						std::cout<<",";
				}
				if (i != P[0].rowdim()-1)
					std::cout<<"],";
				else
					std::cout<<"]";						
			}
			std::cout<<"]);\n";
		}


		std::vector<size_t> masseyblock_left (std::vector<Coefficient> &P) {
			
#ifdef _BM_TIMING
			tSetup.clear();
			tSetup.start();
#endif			
			const size_t length = _container->size ();
		
			const size_t m = _container->rowdim();
			const size_t n = _container->coldim();
		
			// ====================================================
			// Sequence and iterator initialization
			// ====================================================
		
			// Initialization of the sequence iterator
			typename Sequence::const_iterator _iter (_container->begin ());

			// Reservation of memory for the entire sequence
			std::vector<Coefficient> S (length); //,Coefficient(m,n));
					
			Coefficient Unit(m+n,m);
			const Coefficient Zero(m+n,m);
			Element one,zero,mone;
			_F.init(one,1L);
			_F.init(zero,0L);
			_F.init(mone,-1L);
			for (size_t i=0;i<m;i++) 			
				Unit.setEntry(i,i,one);							
			size_t min_mn=(m <n)? m :n;
	
			

			// initialization of discrepancy
			Coefficient Discrepancy(m+n,n);
			for (size_t i=0;i<n;i++)
				Discrepancy.setEntry(i+m,i,one);

			
			// initialization of sigma base
			std::vector<Coefficient> SigmaBase(length);
			SigmaBase.resize(1);
			SigmaBase[0]=Unit;
								
			// initialization of order of sigma base's rows
			std::vector<long> order(m+n,1);
			for (size_t i=0;i<m;++i)
				order[i]=0;

			// initialisation of degree of sigma base's rows
			std::vector<long> degree(m+n,0);
			for (size_t i=0;i<m;++i)
				degree[i]=0;
#ifdef _BM_TIMING
			tSetup.stop();
			ttSetup += tSetup;
			tCheckSequence.clear();
			tCheckSequence.start();
#endif		



			// The first sequence element should be of full rank
			// this is due to the strategy which say that we can compute 
			// only the first column of the approximation of [ S(x) Id]^T
			// since the other colums have always lower degree.			
			if (_BMD.rank(*_iter)< min_mn) 
				throw PreconditionFailed (__FUNCTION__, __LINE__, "Bad random Blocks, abort\n");
			//	cerr<<"\n**************************************************\n";
			//	cerr<<"*** THE FIRST ELEMENT OF SEQUENCE IS SINGULAR  ***\n";
			//	cerr<<"***            ALGORTIHM ABORTED               ***\n";
			//	cerr<<"**************************************************\n";
			//}
			
		
#ifdef _BM_TIMING
			tCheckSequence.stop();
			ttCheckSequence += tCheckSequence;
#endif

			unsigned long early_stop=0;
			long N;
		
			for (N = 0; (N < (long)length) && (early_stop < EARLY_TERM_THRESHOLD) ; ++N, ++_iter) {
									
				// Get the next coefficient in the sequence
				S[N]=*_iter;													
			
#ifdef  _BM_TIMING
				if (N != 0){
					tGetCoeff.stop();
					ttGetCoeff += tGetCoeff;
				}

				tDiscrepancy.clear();
				tDiscrepancy.start();
#endif				

				/*	
				 * Compute the new discrepancy (just updating the first m rows)					
				 */				
				// view of m first rows of SigmaBasis[0]
				Coefficient Sigma(SigmaBase[0],0,0,m,m);
				
				// view of m first rows of Discrepancy
				Coefficient Discr(Discrepancy,0,0,m,n);
								
				_BMD.mul(Discr,Sigma,S[N]);
				for (size_t i=1;i<SigmaBase.size();i++){
					Coefficient  Sigmaview(SigmaBase[i],0,0,m,m);
					_BMD.axpyin(Discr,Sigmaview,S[N-i]);
				}
					
#ifdef _BM_TIMING
				tDiscrepancy.stop();
				ttDiscrepancy += tDiscrepancy;				
#endif	
				
				typename Coefficient::RawIterator _iter_Discr = Discr.rawBegin();

				while ((_F.isZero(*_iter_Discr) && _iter_Discr != Discr.rawEnd()))
					++_iter_Discr;
					
				// maybe there is something to do here
				// increase the last n rows of orders
				// multiply by X the last n rows of SigmaBase
				if (_iter_Discr != Discr.rawEnd())
					early_stop=0;
				else {
					early_stop++;
				}
					
			
#ifdef _BM_TIMING
				tGetPermutation.clear();
				tGetPermutation.start();
#endif						 			
				// Computation of the permutation BPerm1 such that BPerm1.order is in increasing order.
				// order=Perm.order			   			
				std::vector<size_t> Perm1(m+n);			
				for (size_t i=0;i<m+n;++i)
					Perm1[i]=i;	
				if (N>=2) {
					for (size_t i=0;i<m+n;++i) {
						size_t idx_min=i;
						for (size_t j=i+1;j<m+n;++j) 
							if (order[j]< order[idx_min]) 
								idx_min=j;												
						std::swap(order[i],order[idx_min]);				
						Perm1[i]=idx_min;
					}	
				}
				BlasPermutation BPerm1(Perm1);
					
#ifdef _BM_TIMING
				tGetPermutation.stop();
				ttGetPermutation += tGetPermutation;
				tApplyPerm.clear();
				tApplyPerm.start();				
				
#endif		
				// Discrepancy= BPerm1.Discrepancy						
				_BMD.mulin_right(BPerm1,Discrepancy);

#ifdef _BM_TIMING
				tApplyPerm.stop();
				ttApplyPerm += tApplyPerm;
				tLQUP.clear();
				tLQUP.start();				
#endif

#ifdef __CHECK_DISCREPANCY
				std::cout<<"Discrepancy"<<N<<":=Matrix(";
				Discrepancy.write(std::cout,_F,true)<<");"<<std::endl;
#endif
			

				
				// Computation of the LQUP decomposition of the discrepancy
				Coefficient CopyDiscr;
				CopyDiscr=Discrepancy;
				LQUPMatrix<Field> LQUP(_F, CopyDiscr);
				
#ifdef _BM_TIMING
				tLQUP.stop();	
				ttLQUP += tLQUP;

#endif
				// Get the matrix L of LQUP decomposition
				TriangularBlasMatrix<Element> L(m+n,m+n, BlasTag::low, BlasTag::unit );
				LQUP.getL(L);
				
				// Get the tranposed  permutation of Q from LQUP
				BlasPermutation Qt=LQUP.getQ();
			

#ifdef _BM_TIMING
				tGetPermutation.clear();
				tGetPermutation.start();
#endif
				// Computation of permutations BPerm2 such that the last n rows of BPerm2.Qt.Discrepancy are non zero.
				std::vector<size_t> Perm2(m+n);	
				for (size_t i=0;i<n;++i)
					Perm2[i]=m+i;
				for (size_t i=n;i<m+n;++i)
					Perm2[i]=i;					
				BlasPermutation BPerm2(Perm2);			
				
#ifdef _BM_TIMING
				tGetPermutation.stop();
				ttGetPermutation += tGetPermutation;
				tInverseL.clear();
				tInverseL.start();
#endif			
				// compute the inverse of L
				TriangularBlasMatrix<Element> invL (m+n,m+n, BlasTag::low,BlasTag::unit); 
				FFPACK::trinv_left(_F,m+n,L.getPointer(),L.getStride(),invL.getWritePointer(),invL.getStride());

#ifdef _BM_TIMING
				tInverseL.stop();
				ttInverseL += tInverseL;
#endif					

#ifdef 	__CHECK_TRANSFORMATION
				std::cout<<"invL"<<N<<":=Matrix(";
				invL.write(std::cout,_F,true)<<");"<<std::endl;

#endif
				// SigmaBase =  BPerm2.Qt. L^(-1) . BPerm1 . SigmaBase
				for (size_t i=0;i<SigmaBase.size();i++) {
#ifdef _BM_TIMING
					tApplyPerm.clear();
					tApplyPerm.start();					
#endif
					_BMD.mulin_right(BPerm1,SigmaBase[i]);

#ifdef _BM_TIMING
					tApplyPerm.stop();					
					ttApplyPerm +=tApplyPerm;

					tUpdateSigma.clear();
					tUpdateSigma.start();					
#endif
					_BMD.mulin_right(invL,SigmaBase[i]);
#ifdef _BM_TIMING
					tUpdateSigma.stop();
					ttUpdateSigma += tUpdateSigma;
					tApplyPerm.clear();
					tApplyPerm.start(); 										
#endif
					_BMD.mulin_right(Qt,SigmaBase[i]);
					_BMD.mulin_right(BPerm2,SigmaBase[i]);
#ifdef _BM_TIMING
					tApplyPerm.stop(); 
					ttApplyPerm +=tApplyPerm;
#endif
				}							
		

#ifdef _BM_TIMING				
				tApplyPerm.clear();
				tApplyPerm.start();
#endif

				// Apply BPerm2 and Qt to the vector of order and increase by 1 the last n rows
				UnparametricField<long> UF;
				BlasMatrixDomain<UnparametricField<long> > BMDUF(UF);
				BMDUF.mulin_right(Qt,order);
				BMDUF.mulin_right(BPerm2,order);
				BMDUF.mulin_right(BPerm1,degree);
				BMDUF.mulin_right(Qt,degree);
				BMDUF.mulin_right(BPerm2,degree);				
				for (size_t i=m;i<m+n;++i){
					order[i]++;		
					degree[i]++;
				}

#ifdef _BM_TIMING
				tApplyPerm.stop();
				ttApplyPerm += tApplyPerm;
				tShiftSigma.clear();
				tShiftSigma.start();
#endif
				// Multiplying the last n row of SigmaBase by x.
				long max_degree=degree[m];
				for (size_t i=m+1;i<m+n;++i) {
					if (degree[i]>max_degree)
						max_degree=degree[i];
				}			
				size_t size=SigmaBase.size();			
				if (SigmaBase.size()<= (size_t)max_degree)
					{			
						SigmaBase.resize(size+1,Zero);					
						size++;
					}		
				for (int i= (int)size-2;i>=0;i--)
					for (size_t j=0;j<n;j++)
						for (size_t k=0;k<n;++k)
							_F.assign(SigmaBase[i+1].refEntry(m+j,k), SigmaBase[i].getEntry(m+j,k));			

				for (size_t j=0;j<n;j++)
					for (size_t k=0;k<n;++k)
						_F.assign(SigmaBase[0].refEntry(m+j,k),zero);


#ifdef _BM_TIMING
				tShiftSigma.stop();
				ttShiftSigma += tShiftSigma;
#endif


#ifdef __DEBUG_MAPLE	
				std::cout<<"\n\nSigmaBase"<<N<<":= ";
				write_maple(_F,SigmaBase);
				
				std::cout<<"order"<<N<<":=<";
				for (size_t i=0;i<m+n;++i){
					std::cout<<order[i];
					if (i!=m+n-1) std::cout<<",";					
				}
				std::cout<<">;"<<std::endl;
				std::cout<<"degree"<<N<<":=<";
				for (size_t i=0;i<m+n;++i){
					std::cout<<degree[i];
					if (i!=m+n-1) std::cout<<",";					
				}
				std::cout<<">;"<<std::endl;
				
#endif
			
#ifdef __CHECK_LOOP
				std::cout<<"\nCheck validity of current SigmaBase\n";
				std::cout<<"SigmaBase size: "<<SigmaBase.size()<<std::endl;
				std::cout<<"Sequence size:  "<<N+1<<std::endl;
				size_t min_t = (SigmaBase.size() > N+1)? N+1: SigmaBase.size();
				for (size_t i=min_t - 1 ; i<N+1; ++i){
					Coefficient Disc(m+n,n);
					_BMD.mul(Disc,SigmaBase[0],S[i]);
					for (size_t j=1;j<min_t -1;++j)
						_BMD.axpyin(Disc,SigmaBase[j],S[i-j]);
					Disc.write(std::cout,_F)<<std::endl;
				}				
#endif


#ifdef _BM_TIMING
				tNewDiscrepancy.clear();
				tNewDiscrepancy.start();
#endif		
				// Discrepancy= BPerm2.U.P from LQUP
				Coefficient U(m+n,n);
				TriangularBlasMatrix<Element> trU(U,BlasTag::up,BlasTag::nonunit);
				LQUP.getU(trU);	 
				Discrepancy=U;
				BlasPermutation P= LQUP.getP();
				_BMD.mulin_left(Discrepancy,P);
				_BMD.mulin_right(BPerm2,Discrepancy);

#ifdef _BM_TIMING
				tNewDiscrepancy.stop();
				ttNewDiscrepancy+=tNewDiscrepancy;

				// timer in the loop 
				tGetCoeff.clear();	
				tGetCoeff.start();
#endif	

			}
			if ( early_stop == EARLY_TERM_THRESHOLD)
				std::cout<<"Early termination is used: stop at "<<N<<" from "<<length<<" iterations\n\n";
			
#ifdef __PRINT_SEQUENCE	
			std::cout<<"\n\nSequence:= ";
			write_maple(_F,S);
#endif

		

#ifdef __CHECK_SIGMA_RESULT
			std::cout<<"Check SigmaBase application\n";
			for (size_t i=SigmaBase.size()-1 ;i< length ;++i){
				Coefficient res(m+n,n);
				for (size_t k=0;k<SigmaBase.size();++k)
					_BMD.axpyin(res,SigmaBase[k],S[i-k]);
				res.write(std::cout,_F)<<std::endl;
			}

#endif

#ifdef _BM_TIMING
			tGetMinPoly.clear();
			tGetMinPoly.start();
#endif
			// Get the reverse matrix polynomial of the forst m rows of SigmaBase according to degree.
			degree=order;
			long max=degree[0];
			for (size_t i=1;i<m;i++) {
				if (degree[i]>max)
					max=degree[i];
			}
			P = std::vector<Coefficient> (max+1);
			Coefficient tmp(m,m);
			for (long i=0;i< max+1;++i)
				P[i]=tmp;
			
			for (size_t i=0;i<m;i++) 
				for (long j=0;j<=degree[i];j++) 
					for (size_t k=0;k<m;k++) 
						_F.assign(P[degree[i]-j].refEntry(i,k), SigmaBase[j].getEntry(i,k));
#ifdef _BM_TIMING
			tGetMinPoly.stop();
			ttGetMinPoly +=tGetMinPoly;
#endif


#ifdef __CHECK_RESULT
			std::cout<<"Check minimal polynomial application\n";
			bool valid=true;
			for (size_t i=0;i< N - P.size();++i){
				Coefficient res(m,n);
				_BMD.mul(res,P[0],S[i]);
				for (size_t k=1,j=i+1;k<P.size();++k,++j)
					_BMD.axpyin(res,P[k],S[j]);
				for (size_t j=0;j<m*n;++j)
					if (!_F.isZero(*(res.getPointer()+j)))
						valid= false;
				//res.write(std::cout,_F)<<std::endl;				
			}
			if (valid)
				std::cout<<"minpoly is correct\n";
			else
				std::cout<<"minpoly is wrong\n";
#endif

#ifdef __PRINT_MINPOLY
			std::cout<<"MinPoly:=";
			write_maple(_F,P);
			//Coefficient Mat(*_container->getBB());
			//std::cout<<"A:=Matrix(";
			//Mat.write(std::cout,_F,true);
#endif	       

			std::vector<size_t> deg(m);
			for (size_t i=0;i<m;++i)
				deg[i]=degree[i];

			return deg;
		}


		std::vector<size_t> masseyblock_left_rec (std::vector<Coefficient> &P) {

			// Get information of the Sequence (U.A^i.V)
			size_t length = _container->size();
			size_t m, n;
			m = _container->rowdim();
			n = _container->coldim();

			// Set some useful constant
			Element one;
			_F.init(one,1UL);
			const Coefficient Zero(2*m,2*m);

			// Make the Power Serie from  Sequence (U.A^i.V) and Identity
			_container->recompute(); // make sure sequence is already computed
			std::vector<Coefficient> PowerSerie(length);
			typename Sequence::const_iterator _iter (_container->begin ());
			for (size_t i=0;i< length; ++i, ++_iter){
				Coefficient value(2*m,n);
				PowerSerie[i] = value;	
				for (size_t j=0;j<m;++j)
					for (size_t k=0;k<n;++k)
						PowerSerie[i].setEntry(j,k, (*_iter).getEntry(j,k));
			}
			for (size_t j=0;j<n;++j)
				PowerSerie[0].setEntry(m+j, j, one);
#ifdef __PRINT_SEQUENCE	
			std::cout<<"PowerSerie:=";
			write_maple(_F,PowerSerie);		
#endif

			
			// Set the defect to [0 ... 0 1 ... 1]^T
			std::vector<size_t> defect(2*m,0);
			for (size_t i=m;i< 2*m;++i)
				defect[i]=1;
									
			// Prepare SigmaBase
			std::vector<Coefficient> SigmaBase(length,Zero);
			
			// Compute Sigma Base up to the order length - 1
			PM_Basis(SigmaBase, PowerSerie, length-1, defect);
		
			// take the m rows which have lowest defect
			// compute permutation such that first m rows have lowest defect						
			std::vector<size_t> Perm(2*m);			
			for (size_t i=0;i<2*m;++i)
				Perm[i]=i;			
			for (size_t i=0;i<2*m;++i) {
				size_t idx_min=i;
				for (size_t j=i+1;j<2*m;++j) 
					if (defect[j]< defect[idx_min]) 
						idx_min=j;												
				std::swap(defect[i],defect[idx_min]);				
				Perm[i]=idx_min;
			}	
			BlasPermutation BPerm(Perm);
			
			// Apply BPerm to the Sigma Base
			for (size_t i=0;i<SigmaBase.size();++i)
				_BMD.mulin_right(BPerm,SigmaBase[i]);
			
			//std::cout<<"SigmaBase:=";
			//write_maple(_F,SigmaBase);
						
			// Compute the reverse polynomial of SigmaBase according to defect of each row
			size_t max=defect[0];
			for (size_t i=0;i<m;++i)
				if (defect[i] > max)
					max=defect[i];
			
			P = std::vector<Coefficient> (max+1);
			Coefficient tmp(m,m);
			for (size_t i=0;i< max+1;++i)
				P[i]=tmp;
			for (size_t i=0;i<m;i++) 
				for (size_t j=0;j<=defect[i];j++) 
					for (size_t k=0;k<m;k++) 
						_F.assign(P[defect[i]-j].refEntry(i,k), SigmaBase[j].getEntry(i,k));

#ifdef __CHECK_RESULT
			std::cout<<"Check minimal polynomial application\n";
			_container->recompute();
			typename Sequence::const_iterator _ptr (_container->begin ());
			for (size_t i=0;i< length; ++i, ++_ptr){				
				PowerSerie[i] = *_ptr;	
			}			
			bool valid=true;
			for (size_t i=0;i< length - P.size();++i){
				Coefficient res(m,n);
				Coefficient Power(PowerSerie[i],0,0,m,n);
				_BMD.mul(res,P[0],Power);
				for (size_t k=1,j=i+1;k<P.size();++k,++j){
					Coefficient Powerview(PowerSerie[j],0,0,m,n);
					_BMD.axpyin(res,P[k],Powerview);
				}
				for (size_t j=0;j<m*n;++j)
					if (!_F.isZero(*(res.getPointer()+j)))
						valid= false;
				//res.write(std::cout,_F)<<std::endl;				
			}
			if (valid)
				std::cout<<"minpoly is correct\n";
			else
				std::cout<<"minpoly is wrong\n";
#endif

#ifdef __PRINT_MINPOLY
			std::cout<<"MinPoly:=";
			write_maple(_F,P);
			//Coefficient Mat(*_container->getBB());
			//std::cout<<"A:=Matrix(";
			//Mat.write(std::cout,_F,true);
#endif		
			std::vector<size_t> degree(m);
			for (size_t i=0;i<m;++i)
				degree[i] = defect[i];
			return degree;
		}

		
		// Computation of a minimal Sigma Base of a Power Serie up to a degree
		// according to a vector of defect.
		// algorithm is from Giorgi, Jeannerod and Villard  ISSAC'03
		//
		// SigmaBase must be already allocated with degree+1 elements
		
		void PM_Basis(std::vector<Coefficient>     &SigmaBase,
			      std::vector<Coefficient>    &PowerSerie, 
			      size_t                           degree, 
			      std::vector<size_t>             &defect) {
						
			size_t m,n;
			m = PowerSerie[0].rowdim();
			n = PowerSerie[0].coldim();
			Element one;
			_F.init(one,1UL);
			const Coefficient ZeroSigma(m,m);
			const Coefficient ZeroSerie(m,n);

			if (degree == 0) {
				Coefficient Identity(m,m);
				for (size_t i=0;i< m;++i)
					Identity.setEntry(i,i,one);
				SigmaBase[0]=Identity;
			}
			else {
				if (degree == 1) {
#ifdef _BM_TIMING				
					tMBasis.clear();
					tMBasis.start();
#endif
					M_Basis(SigmaBase, PowerSerie, degree, defect);
#ifdef _BM_TIMING
					tMBasis.stop();
					ttMBasis += tMBasis;
#endif			
				}
				else {
					size_t degree1,degree2;
					degree1 = (degree >> 1) + (degree & 1);
					degree2 = degree - degree1;									

					// Compute Sigma Base of half degree
					std::vector<Coefficient> Sigma1(degree1+1,ZeroSigma);
					std::vector<Coefficient> Serie1(degree1+1);
					for (size_t i=0;i< degree1+1;++i)
						Serie1[i] = PowerSerie[i];
					PM_Basis(Sigma1, Serie1, degree1, defect);
#ifdef _BM_TIMING				
					tUpdateSerie.clear();
					tUpdateSerie.start();
#endif
					// Compute Serie2 = x^(-degree1).Sigma.PowerSerie mod x^degree2
					std::vector<Coefficient> Serie2(degree1+1,ZeroSerie);										
					ComputeNewSerie(Serie2,Sigma1,PowerSerie, degree1, degree2);
#ifdef _BM_TIMING				
					tUpdateSerie.stop();
					ttUpdateSerie += tUpdateSerie;
#endif
					// Compute Sigma Base of half degree from updated Power Serie					
					std::vector<Coefficient> Sigma2(degree2+1,ZeroSigma);
					PM_Basis(Sigma2, Serie2, degree2, defect);
						
#ifdef _BM_TIMING				
					tBasisMultiplication.clear();
					tBasisMultiplication.start();
#endif
					// Compute the whole Sigma Base through the product 
					// of the Sigma Basis Sigma1 x Sigma2										
					MulSigmaBasis(SigmaBase,Sigma2,Sigma1);						
#ifdef _BM_TIMING				
					tBasisMultiplication.stop();
					ttBasisMultiplication += tBasisMultiplication;
#endif

				}
			}
		}


		// Computation of a minimal Sigma Base of a Power Serie up to length
		// algorithm is from Giorgi, Jeannerod and Villard  ISSAC'03		
		void M_Basis(std::vector<Coefficient>     &SigmaBase,
			     std::vector<Coefficient>    &PowerSerie, 
			     size_t                           length, 
			     std::vector<size_t>             &defect) {

			// Get the dimension of matrices inside 
			// the Matrix Power Serie
			size_t m,n;
			m = PowerSerie[0].rowdim();
			n = PowerSerie[0].coldim();

			// Set some useful constants
			const Coefficient Zero(m,m);
			Element one, zero;
			_F.init(one,1UL);
			_F.init(zero,0UL);
			
			// Reserve memory for the Sigma Base and set SigmaBase[0] to Identity
			SigmaBase.reserve(length+1);
			SigmaBase.resize(1);
			Coefficient Identity(m,m);
			for (size_t i=0;i< m;++i)
				Identity.setEntry(i,i,one);
			SigmaBase[0]=Identity;
			
			// Keep track on Sigma Base's row degree
			std::vector<size_t> degree(m,0);
			for (size_t i=0;i<n;++i)
				degree[i]=0;
			

			// Compute the minimal Sigma Base of the PowerSerie up to length
			for (size_t k=0; k< length; ++k) {

				// compute BPerm1 such that BPerm1.defect is in increasing order
				std::vector<size_t> Perm1(m);			
				for (size_t i=0;i<m;++i)
					Perm1[i]=i;			
				for (size_t i=0;i<m;++i) {
					size_t idx_min=i;
					for (size_t j=i+1;j<m;++j) 
						if (defect[j]< defect[idx_min]) 
							idx_min=j;												
					std::swap(defect[i], defect[idx_min]);				
					Perm1[i]=idx_min;
				}					
				BlasPermutation BPerm1(Perm1);

				// Apply Bperm1 to the current SigmaBase
				for (size_t i=0;i<SigmaBase.size();++i)
					_BMD.mulin_right(BPerm1,SigmaBase[i]);

				// Compute Discrepancy
				Coefficient Discrepancy(m,n);								
				_BMD.mul(Discrepancy,SigmaBase[0],PowerSerie[k]);
				for (size_t i=1;i<SigmaBase.size();i++){
					_BMD.axpyin(Discrepancy,SigmaBase[i],PowerSerie[k-i]);
				}
				
				// Compute LQUP of Discrepancy
				LQUPMatrix<Field> LQUP(_F,Discrepancy);

				// Get L from LQUP
				TriangularBlasMatrix<Element> L(m, m, BlasTag::low, BlasTag::unit);
				LQUP.getL(L);

				// get the transposed permutation of Q from LQUP
				BlasPermutation Qt =LQUP.getQ();

				// Compute the inverse of L
				TriangularBlasMatrix<Element> invL(m, m, BlasTag::low, BlasTag::unit);
				FFPACK::trinv_left(_F,m,L.getPointer(),L.getStride(),invL.getWritePointer(),invL.getStride());

				// Update Sigma by L^(-1)
				// Sigma = L^(-1) . Sigma
				for (size_t i=0;i<SigmaBase.size();++i) 
					_BMD.mulin_right(invL,SigmaBase[i]);

				//std::cout<<"BaseBis"<<k<<":=";
				//write_maple(_F,SigmaBase);
				// Increase by degree and defect according to row choosen as pivot in LQUP
				for (size_t i=0;i<n;++i){
					defect[*(Qt.getPointer()+i)]++;
					degree[*(Qt.getPointer()+i)]++;
				}
								
				size_t max_degree=degree[*(Qt.getPointer())];
				for (size_t i=0;i<n;++i) {
					if (degree[*(Qt.getPointer()+i)]>max_degree)
						max_degree=degree[*(Qt.getPointer()+i)];
				}	
			
		
				size_t size=SigmaBase.size();
				if (SigmaBase.size()<= max_degree+1)
					{					
						SigmaBase.resize(size+1,Zero);					
						size++;
					}				
				// Mulitply by x the rows of Sigma involved as pivot in LQUP
				for (size_t i=0;i<n;++i){
					for (int j= (int) size-2;j>=0; --j){						
						for (size_t l=0;l<m;++l)
							_F.assign(SigmaBase[j+1].refEntry(*(Qt.getPointer()+i),l),
								  SigmaBase[j].getEntry(*(Qt.getPointer()+i),l));			
					}
					for (size_t l=0;l<m;++l)
						_F.assign(SigmaBase[0].refEntry(*(Qt.getPointer()+i),l),zero);
				}

				//std::cout<<"Base"<<k<<":=";
				//write_maple(_F,SigmaBase);
			}
			//std::cout<<"defect: ";
			//for (size_t i=0;i<m;++i)
			//	std::cout<<defect[i]<<" ";
			//std::cout<<std::endl;
			
			//std::cout<<"SigmaBase"<<length<<":=";
			//write_maple(_F,SigmaBase);

		}


		// compute the middle product of A [1..n].B[1..2n]
		// using Karatsuba multiplication
		// algorithm is that of Hanrot, Quercia and Zimmermann 2002

		void MP_Karatsuba(std::vector<Coefficient> &C, const std::vector<Coefficient> &A, const std::vector<Coefficient> &B){
		
			if (A.size() == 1)
				_BMD.mul(C[0],A[0],B[0]);
			else {				
				size_t k0= A.size()>>1;
				size_t k1= A.size()-k0;

				size_t m = B[0].rowdim();
				size_t n = B[0].coldim();

				const Coefficient Zero(m,n);
				std::vector<Coefficient> alpha(k1,Zero), beta(k1,Zero), gamma(k0,Zero);
				std::vector<Coefficient> A_low(k0), A_high(k1), B1(2*k1-1), B2(2*k1-1);
				
				for (size_t i=0;i<k0;++i)
					A_low[i] = A[i];

				for (size_t i=k0;i<A.size();++i)
					A_high[i-k0] = A[i];

				for (size_t i=0;i<2*k1-1;++i){
					B1[i] = B[i];
					B2[i] = B[i+k1];
					_MD.addin(B1[i],B2[i]);
				}		
				MP_Karatsuba(alpha, A_high, B1);			

				if (k0 == k1) {
					for (size_t i=0;i<k1;++i)
						_MD.subin(A_high[i],A_low[i]);
					MP_Karatsuba(beta, A_high, B2);					
				}
				else {
					for (size_t i=1;i<k1;++i)
						_MD.subin(A_high[i],A_low[i-1]);
					MP_Karatsuba(beta, A_high, B2);
				}				
				
				std::vector<Coefficient> B3(2*k0-1,Zero);
				for (size_t i=0;i<2*k0-1;++i)
					_MD.add(B3[i],B[i+2*k1],B[i+k1]);
				
				MP_Karatsuba(gamma, A_low, B3);

				for (size_t i=0;i<k1;++i)
					_MD.sub(C[i],alpha[i],beta[i]);
				
				for (size_t i=0;i<k0;++i){
					C[k1+i]=gamma[i];
					_MD.addin(C[k1+i],beta[i]);								
				}
			}
		}


		// Multiply a Power Serie by a Sigma Base.
		// only affect coefficients of the Power Serie between degree1 and degree2
		void ComputeNewSerie(std::vector<Coefficient>          &NewSerie, 
				     const std::vector<Coefficient>   &SigmaBase, 
				     const std::vector<Coefficient>    &OldSerie,
				     size_t                              degree1,
				     size_t                              degree2){						
			

			// degree1 >= degree2
			//size_t size = 2*degree1 + 1;
					
			const Coefficient ZeroSerie (OldSerie[0].rowdim(), OldSerie[0].coldim());
			const Coefficient ZeroBase  (SigmaBase[0].rowdim(), SigmaBase[0].coldim());

			// Work on a copy of the old  Serie (increase size by one for computation of middle product)
			std::vector<Coefficient> Serie(OldSerie.size()+1,ZeroSerie);
			for (size_t i=0;i< OldSerie.size();++i)
				Serie[i] = OldSerie[i];

			// Work on a copy of the Sigma Base 
			std::vector<Coefficient> Sigma(SigmaBase.size());
			for (size_t i=0;i<SigmaBase.size();++i){
				Sigma[i] = SigmaBase[i];
			}

			MP_Karatsuba(NewSerie, Sigma, Serie);

			//std::vector<Coefficient> NewPowerSerie(SigmaBase.size()+OldSerie.size(), Zero);
			//MulSigmaBasis(NewPowerSerie, Sigma, Serie);		       
			//for (size_t i=0;i<degree2;++i)
			//	NewSerie[i] = NewPowerSerie[i+degree1];	

			
		}
		

		// matrix polynomial multiplication
		// using Karatsuba's algorithm		
		void MulPolyMatrix(std::vector<Coefficient> &C, size_t shiftC,
				   std::vector<Coefficient> &A, size_t shiftA, size_t degA,
				   std::vector<Coefficient> &B, size_t shiftB, size_t degB){
			
			const Coefficient ZeroC(C[0].rowdim(), C[0].coldim());
			const Coefficient ZeroA(A[0].rowdim(), A[0].coldim());
			const Coefficient ZeroB(B[0].rowdim(), B[0].coldim());
			
			if ((degA == 1) || (degB == 1)) {
				
				if ((degA == 1) && (degB == 1))
					_BMD.mul(C[shiftC],A[shiftA],B[shiftB]); 
				else 
					if (degA == 1) 
						for (size_t i=0;i< degB;++i)
							_BMD.mul(C[shiftC+i],A[shiftA],B[shiftB+i]);
					else 
						for (size_t i=0;i< degA;++i)
							_BMD.mul(C[shiftC+i],A[shiftA+i],B[shiftB]);
			}
			else {
				size_t degA_low, degA_high, degB_low, degB_high, half_degA, half_degB, degSplit;
				half_degA= (degA & 1) + degA >>1;
				half_degB= (degB & 1) + degB >>1;
				degSplit= (half_degA > half_degB) ? half_degA : half_degB;
				
				degB_low = (degB < degSplit) ? degB : degSplit;
				degA_low = (degA < degSplit) ? degA : degSplit;
				degA_high= degA - degA_low;
				degB_high= degB - degB_low;
				
				// multiply low degrees
				MulPolyMatrix(C, shiftC, A, shiftA, degA_low, B, shiftB, degB_low);   
				
				// multiply high degrees (only if they are both different from zero)
				if ((degA_high !=0) && (degB_high != 0)) {	
					MulPolyMatrix(C, shiftC+(degSplit << 1), A, shiftA+degSplit, degA_high, B, shiftB+degSplit, degB_high);
				}
				
				// allocate space for summation of low and high degrees
				std::vector<Coefficient> A_tmp(degA_low,ZeroA);
				std::vector<Coefficient> B_tmp(degB_low,ZeroB);
				std::vector<Coefficient> C_tmp(degA_low+degB_low-1,ZeroC);
				
				// add low and high degrees of A
				for (size_t i=0;i<degA_low;++i)
					A_tmp[i]=A[shiftA+i];
				if ( degA_high != 0) 
					for (size_t i=0;i<degA_high;++i)
						_MD.addin(A_tmp[i],A[shiftA+degSplit+i]);	
				
				// add low and high degrees of B
				for (size_t i=0;i<degB_low;++i)
					B_tmp[i]=B[shiftA+i];
				if ( degB_high != 0)
					for (size_t i=0;i<degB_high;++i)
						_MD.addin(B_tmp[i],B[shiftB+degSplit+i]);
				
				//  multiply the sums
				MulPolyMatrix(C_tmp, 0, A_tmp, 0, degA_low, B_tmp, 0, degB_low);
				
				// subtract the low product from the product of sums
				for (size_t i=0;i< C_tmp.size();++i)
					_MD.subin(C_tmp[i], C[shiftC+i]);	
				
				// subtract the high product from the product of sums
				if ((degA_high !=0) && (degB_high != 0))
					for (size_t i=0;i< degA_high+degB_high-1; ++i)
						_MD.subin(C_tmp[i], C[shiftC+(degSplit << 1)+i]);
				
				// add the middle term of the product
				size_t mid= (degA_low+degB_high > degB_low+degA_high)? degA_low+degB_high :degB_low+degA_high;
				for (size_t i=0;i< mid-1; ++i)
					_MD.addin(C[shiftC+degSplit+i], C_tmp[i]);											
			}		    
		}				
			
		// Multiply two Sigma Basis
		// in fact this is matrix polynomial multiplication
		// we assume that we can modify each operand 
		// since only result will be used
		void MulSigmaBasis(std::vector<Coefficient> &C, 
				   std::vector<Coefficient> &A,
				   std::vector<Coefficient> &B){
			//std::cout<<"C=A*B: "<<C.size()<<" "<<A.size()<<" "<<B.size()<<std::endl;
			MulPolyMatrix(C, 0, A, 0, A.size(), B, 0, B.size());
			
			//for (size_t i=0;i<A.size();++i)
			//	for (size_t j=0;j<B.size();++j)
			//		_BMD.axpyin(C[i+j],A[i],B[j]);
		}
					
	}; //end of class BlockMasseyDomain
	
} // end of namespace LinBox
	
#endif // __LINBOX_massey_block_domain_H
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

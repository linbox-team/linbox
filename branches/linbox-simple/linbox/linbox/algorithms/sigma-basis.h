/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/algorithms/sigma-basis.h
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by Pascal Giorgi pgiorgi@uwaterlo.ca
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



#ifndef __SIGMA_BASIS_H
#define __SIGMA_BASIS_H

#include <vector>
#include <iostream>
#include <iomanip>

#include <linbox/util/commentator.h>
#include <linbox/util/timer.h>
#include <linbox/algorithms/blas-domain.h>
#include <linbox/blackbox/dense.h>
#include <linbox/field/unparametric.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/matrix/factorized-matrix.h>

#include <linbox/algorithms/matpoly-mult.h>

#include <linbox/util/timer.h>


namespace LinBox {


	template<class _Field>
	class SigmaBasis {
		
	public:
		typedef _Field                           Field;
		typedef typename Field::Element        Element;
		typedef BlasMatrix<Element>        Coefficient;

	private:
		Field                             _F;
		BlasMatrixDomain<Field>         _BMD;
		MatrixDomain<Field>              _MD;
		std::vector<Coefficient>     &_Serie;
		
#ifdef _BM_TIMING
		mutable Timer   
		       ttNewDiscrepancy,  tNewDiscrepancy,
			ttShiftSigma,      tShiftSigma,
			ttApplyPerm,       tApplyPerm, 
			ttUpdateSigma,     tUpdateSigma,
			ttInverseL,        tInverseL,
			ttGetPermutation,  tGetPermutation,
			ttLQUP,            tLQUP,
			ttDiscrepancy,     tDiscrepancy,
			ttGetCoeff,        tGetCoeff,
			ttCheckSequence,   tCheckSequence,
			ttMBasis,          tMBasis,
			ttUpdateSerie,     tUpdateSerie,
			ttBasisMultiplication, tBasisMultiplication,
			ttCopyingData,     tCopyingData,
			Total;

		void clearTimer() {
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
	public:
		void printTimer() {
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
	public:

		SigmaBasis(const Field &F, std::vector<Coefficient> &PowerSerie) : _F(F), _BMD(F), _MD(F), _Serie(PowerSerie) 
		{
#ifdef  _BM_TIMING
			clearTimer();
#endif
		}

		void right_basis(std::vector<Coefficient>     &SigmaBase,
				size_t                           degree, 
				std::vector<size_t>             &defect) {

			
			size_t length=_Serie.size();
			size_t m = _Serie[0].rowdim();
			size_t n = _Serie[0].coldim();  
			
			const Coefficient TransposedZero(n,m);

			// take the transposed of the Serie and use PM_Basis
			std::vector<Coefficient> TransposedSerie (length, TransposedZero);

			for (size_t k=0;k<length; ++k)
				for (size_t i=0;i<m;++i)
					for (size_t j=0;j<n;++j)
						TransposedSerie[k].setEntry(j,i, _Serie[k].getEntry(i,j));
			
			const Coefficient Zero(n, n);	
		
			// Prepare SigmaBase
			SigmaBase.resize(degree+1, Zero);

			// Compute Sigma Base up to the order degree
			PM_Basis(SigmaBase, TransposedSerie, degree, defect);					
		}

		void left_basis(std::vector<Coefficient>      &SigmaBase,
				 size_t                           degree, 
				 std::vector<size_t>             &defect) {

			const size_t m = _Serie[0].rowdim();
			const Coefficient Zero(m, m);

			// Prepare SigmaBase
			SigmaBase.resize(degree+1, Zero);
			
			// Compute Sigma Base up to the order degree
			PM_Basis(SigmaBase, _Serie, degree, defect);			
		}


		void multi_left_basis(std::vector<Coefficient>      &SigmaBase1,
				      size_t                            degree1, 
				      std::vector<size_t>              &defect1,
				      std::vector<Coefficient>      &SigmaBase2,
				      size_t                            degree2, 
				      std::vector<size_t>              &defect2) {
			
			linbox_check(degree1 < degree2);

			const size_t m = _Serie[0].rowdim();
			const size_t n = _Serie[0].coldim();
			const Coefficient Zero(m, m);
			const Coefficient ZeroSerie(m,n);

			
			// Prepare SigmaBases
			SigmaBase1.resize(degree1+1, Zero);
			SigmaBase2.resize(degree2+1, Zero);
			
			// Compute Sigma Base up to degree1
			PM_Basis(SigmaBase1, _Serie, degree1, defect1);	

			// copy defect
			defect2 = defect1;
						
			// Update Serie to compute Sigma base up to degree1 - degree2			
			std::vector<Coefficient> Serie2(degree2-degree1+1,ZeroSerie);	
			ComputeNewSerieClassic(Serie2,SigmaBase1, _Serie, degree1, degree2-degree1);
		
			// Compute Sigma Base up to degree2 
			std::vector<Coefficient> Sigma2(degree2-degree1+1, Zero);
			PM_Basis(Sigma2, Serie2, degree2-degree1, defect2);	

			// multily SigmaBase1 by Sigma
			MulSigmaBasis(SigmaBase2,Sigma2,SigmaBase1);		
		}


		void multi_right_basis(std::vector<Coefficient>      &SigmaBase1,
				       size_t                            degree1, 
				       std::vector<size_t>              &defect1,
				       std::vector<Coefficient>      &SigmaBase2,
				       size_t                            degree2, 
				       std::vector<size_t>              &defect2) {
			
			linbox_check(degree1 < degree2);

			size_t length=_Serie.size();
			size_t m = _Serie[0].rowdim();
			size_t n = _Serie[0].coldim();  			
			const Coefficient TransposedZero(n,m);
			const Coefficient Zero(n, n);

			// take the transposed of the Serie and use PM_Basis
			std::vector<Coefficient> TransposedSerie (length, TransposedZero);

			for (size_t k=0;k<length; ++k)
				for (size_t i=0;i<m;++i)
					for (size_t j=0;j<n;++j)
						TransposedSerie[k].setEntry(j,i, _Serie[k].getEntry(i,j));
							
			// Prepare SigmaBases
			SigmaBase1.resize(degree1+1, Zero);
			SigmaBase2.resize(degree2+1, Zero);
			
			// Compute Sigma Base up to degree1
			PM_Basis(SigmaBase1, TransposedSerie, degree1, defect1);
			
			// copy defect
			defect2 = defect1;

			// Update Serie to compute Sigma base up to degree1 - degree2			
			std::vector<Coefficient> Serie2(degree2-degree1+1,TransposedZero);	
			ComputeNewSerieClassic(Serie2,SigmaBase1, TransposedSerie, degree1, degree2-degree1);
		
			// Compute Sigma Base up to degree2 
			std::vector<Coefficient> Sigma2(degree2-degree1+1, Zero);
			PM_Basis(Sigma2, Serie2, degree2-degree1, defect2);	

			// multily SigmaBase1 by Sigma
			MulSigmaBasis(SigmaBase2,Sigma2,SigmaBase1);	
		}



		// function to compute the left denominator from Matrix Pade Approximant
		// compute Q(x) in Q(x).S(x) - R(x) = O(x^degree).

		void left_PadeMatrix (std::vector<Coefficient>            &Approx, 			
				      size_t                               degree, 
				      std::vector<size_t>                 &defect) {

			PadeApproximant(Approx, _Serie, degree, defect);
		}
	
		// function to compute the right denominator from Matrix Pade Approximant
		// compute Q(x) in S(x).Q(x) - R(x) = O(x^degree).
		
		void right_PadeMatrix (std::vector<Coefficient>            &Approx, 
				       size_t                               degree, 
				       std::vector<size_t>                 &defect) {

			size_t deg = _Serie.size();
			size_t m   = _Serie[0].rowdim();
			size_t n   = _Serie[0].coldim();
			

			// transpose the PowerSerie
			const Coefficient Zero(m,n);
			std::vector<Coefficient> TransPowerSerie(deg,Zero);
			for (size_t i=0;i<deg;++i){
				for (size_t j=0;j<m;++j)
					for (size_t k=0;k<n;++k)
						TransPowerSerie[i].setEntry(j,k, _Serie[i].getEntry(k,j));
			}

			// call left pade
			PadeApproximant(Approx, TransPowerSerie, degree, defect);

			// transpose the result
			size_t d= Approx.size();
			for (size_t k=0;k<d;++k)
				for (size_t i=0;i<m;++i)
					for (size_t j=i;j<n;++j)
						std::swap(Approx[k].refEntry(i,j), Approx[k].refEntry(j,i));				
		}

			
		// function to compute two right denominators from Matrix Pade Approximant at two different degrees
		// compute Q1(x), Q2(x) such that  Q1(x).S(x) - R1(x) = O(x^degree1) and Q2(x).S(x) - R2(x) = O(x^degree2) 
		
		void multi_left_PadeMatrix (std::vector<Coefficient>            &Approx1, 			
					    size_t                               degree1,
					    std::vector<Coefficient>            &Approx2, 			
					    size_t                               degree2,
					    std::vector<size_t>                  &defect) {
			
			MultiPadeApproximant(Approx1, degree1, Approx2, degree2, _Serie, defect);
		}
		
		// function to compute two right denominators from Matrix Pade Approximant at two different degrees
		// compute Q1(x), Q2(x) such that  S(x).Q1(x) - R1(x) = O(x^degree1) and S(x).Q2(x) - R2(x) = O(x^degree2) 
		
		void multi_right_PadeMatrix (std::vector<Coefficient>            &Approx1, 			
					     size_t                               degree1,
					     std::vector<Coefficient>            &Approx2, 			
					     size_t                               degree2,
					     std::vector<size_t>                  &defect) {

			size_t deg = _Serie.size();
			size_t m   = _Serie[0].rowdim();
			size_t n   = _Serie[0].coldim();
			

			// transpose the PowerSerie
			const Coefficient Zero(m,n);
			std::vector<Coefficient> TransPowerSerie(deg,Zero);
			for (size_t i=0;i<deg;++i){
				for (size_t j=0;j<m;++j)
					for (size_t k=0;k<n;++k)
						TransPowerSerie[i].setEntry(j,k, _Serie[i].getEntry(k,j));
			}

			// call multi pade
			MultiPadeApproximant(Approx1, degree1, Approx2, degree2, TransPowerSerie, defect);

			// transpose the results
			size_t d;
			d = Approx1.size();
			for (size_t k=0;k<d;++k)
				for (size_t i=0;i<m;++i)
					for (size_t j=i;j<n;++j)
						std::swap(Approx1[k].refEntry(i,j), Approx1[k].refEntry(j,i));				
			d = Approx2.size();
			for (size_t k=0;k<d;++k)
				for (size_t i=0;i<m;++i)
					for (size_t j=i;j<n;++j)
						std::swap(Approx2[k].refEntry(i,j), Approx2[k].refEntry(j,i));				
		}

		
	public:
		//protected:

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
				//write_maple(F,SigmaBase);
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
			
			PolynomialMatrixDomain<Field, std::vector<Coefficient> > PM_domain(_F);
			PM_domain.midproduct(NewSerie, Sigma, Serie);
		}
		

		void ComputeNewSerieClassic(std::vector<Coefficient>          &NewSerie, 
					    const std::vector<Coefficient>   &SigmaBase, 
					    const std::vector<Coefficient>    &OldSerie,
					    size_t                              degree1,
					    size_t                              degree2){						
						
			for (size_t i=0;i<degree2+1;++i)
				for (size_t j=0;j<degree1+1;++j)
					_BMD.axpyin(NewSerie[i],SigmaBase[j],OldSerie[degree1+i-j]);
		}


		


		// Multiply two Sigma Basis
		void MulSigmaBasis(std::vector<Coefficient> &C, 
				   std::vector<Coefficient> &A,
				   std::vector<Coefficient> &B){

			PolynomialMatrixDomain<Field, std::vector<Coefficient> > PM_domain(_F);
			PM_domain.mul(C,A,B);			
		}			
		
		
		void PadeApproximant (std::vector<Coefficient>            &Approx, 
				      const std::vector<Coefficient>  &PowerSerie, 
				      size_t                               length, 
				      std::vector<size_t>                 &defect) {//pph
			

			const size_t m = PowerSerie[0].rowdim();
			const size_t n = PowerSerie[0].coldim();
		
				
			// Initialization of the Serie iterator
			typename std::vector<Coefficient>::const_iterator _iter (PowerSerie.begin ());

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

		
			// The first sequence element should be of full rank
			// this is due to the strategy which say that we can compute 
			// only the first column of the approximation of [ S(x) Id]^T
			// since the other colums have always lower degree.			
			if (_BMD.rank(*_iter)< min_mn) {
				std::cout<<"LinBox ERROR: constant term in the Power Serie is singular\n";
				throw PreconditionFailed (__FUNCTION__, __LINE__, "Bad random Blocks, abort\n");
			}
			
			unsigned long early_stop=0;
			long N;
			
			for (N = 0; (N < (long)length) && (early_stop < 20) ; ++N) {
									
			
				/*	
				 * Compute the new discrepancy (just updating the first m rows)					
				 */				
				// view of m first rows of SigmaBasis[0]
				Coefficient Sigma(SigmaBase[0],0,0,m,m);
				
				// view of m first rows of Discrepancy
				Coefficient Discr(Discrepancy,0,0,m,n);
								
				_BMD.mul(Discr,Sigma, PowerSerie[N]);
				for (size_t i=1;i<SigmaBase.size();i++){
					Coefficient  Sigmaview(SigmaBase[i],0,0,m,m);
					_BMD.axpyin(Discr,Sigmaview,PowerSerie[N-i]);
				}
										
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
					
					 			
				// Computation of the permutation BPerm1 such that BPerm1.order is in increasing order.
				// order=Perm.order			   			
				std::vector<size_t> Perm1(m+n);			
				for (size_t i=0;i<m+n;++i)
					Perm1[i]=i;	
				if (N>=1) {
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
					

				// Discrepancy= BPerm1.Discrepancy						
				_BMD.mulin_right(BPerm1,Discrepancy);


				// Computation of the LQUP decomposition of the discrepancy
				Coefficient CopyDiscr;
				CopyDiscr=Discrepancy;
				LQUPMatrix<Field> LQUP(_F, CopyDiscr);

				// Get the matrix L of LQUP decomposition
				TriangularBlasMatrix<Element> L(m+n,m+n, BlasTag::low, BlasTag::unit );
				LQUP.getL(L);
				
				// Get the tranposed  permutation of Q from LQUP
				BlasPermutation Qt=LQUP.getQ();
			
				// Computation of permutations BPerm2 such that the last n rows of BPerm2.Qt.Discrepancy are non zero.
				std::vector<size_t> Perm2(m+n);	
				for (size_t i=0;i<n;++i)
					Perm2[i]=m+i;
				for (size_t i=n;i<m+n;++i)
					Perm2[i]=i;					
				BlasPermutation BPerm2(Perm2);			
				
				// compute the inverse of L
				TriangularBlasMatrix<Element> invL (m+n,m+n, BlasTag::low,BlasTag::unit); 
				FFPACK::trinv_left(_F,m+n,L.getPointer(),L.getStride(),invL.getWritePointer(),invL.getStride());

				
				// SigmaBase =  BPerm2.Qt. L^(-1) . BPerm1 . SigmaBase
				for (size_t i=0;i<SigmaBase.size();i++) {
					_BMD.mulin_right(BPerm1,SigmaBase[i]);
					_BMD.mulin_right(invL,SigmaBase[i]);
					_BMD.mulin_right(Qt,SigmaBase[i]);
					_BMD.mulin_right(BPerm2,SigmaBase[i]);
				}							
		
				
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


				// Discrepancy= BPerm2.U.P from LQUP
				Coefficient U(m+n,n);
				TriangularBlasMatrix<Element> trU(U,BlasTag::up,BlasTag::nonunit);
				LQUP.getU(trU);	 
				Discrepancy=U;
				BlasPermutation P= LQUP.getP();
				_BMD.mulin_left(Discrepancy,P);
				_BMD.mulin_right(BPerm2,Discrepancy);
								
			}
			if ( early_stop == 20)
				std::cout<<"Early termination is used: stop at "<<N<<" from "<<length<<" iterations\n\n";
			
			// extract the first m rows of SigmaBase
			degree=order;
			long max=degree[0];
			for (size_t i=1;i<m;i++) {
				if (degree[i]>max)
					max=degree[i];
			}
			
			const Coefficient AZero(m,m);
			Approx.resize(max+1, AZero);		
			
			for (size_t i=0;i<m;i++) 
				for (long j=0;j<=degree[i];j++) 
					for (size_t k=0;k<m;k++) 
						_F.assign(Approx[j].refEntry(i,k), SigmaBase[j].getEntry(i,k));	

			
		}

		void MultiPadeApproximant (std::vector<Coefficient>            &Approx1,
					   size_t                               degree1,
					   std::vector<Coefficient>            &Approx2,
					   size_t                               degree2,
					   const std::vector<Coefficient>   &PowerSerie,
					   std::vector<size_t>                  &defect) {
					   

			linbox_check(degree1 < degree2);

			size_t length = degree2;
			const size_t m = PowerSerie[0].rowdim();
			const size_t n = PowerSerie[0].coldim();
		
				
			// Initialization of the Serie iterator
			typename std::vector<Coefficient>::const_iterator _iter (PowerSerie.begin ());

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

		
			// The first sequence element should be of full rank
			// this is due to the strategy which say that we can compute 
			// only the first column of the approximation of [ S(x) Id]^T
			// since the other colums have always lower degree.			
			if (_BMD.rank(*_iter)< min_mn) {
				std::cout<<"LinBox ERROR: constant term in the Power Serie is singular\n";
				throw PreconditionFailed (__FUNCTION__, __LINE__, "Bad random Blocks, abort\n");
			}
			
			unsigned long early_stop=0;
			long N;
			
			for (N = 0; (N < (long)length) && (early_stop < 20) ; ++N) {
									
			
				/*	
				 * Compute the new discrepancy (just updating the first m rows)					
				 */				
				// view of m first rows of SigmaBasis[0]
				Coefficient Sigma(SigmaBase[0],0,0,m,m);
				
				// view of m first rows of Discrepancy
				Coefficient Discr(Discrepancy,0,0,m,n);
								
				_BMD.mul(Discr,Sigma, PowerSerie[N]);
				for (size_t i=1;i<SigmaBase.size();i++){
					Coefficient  Sigmaview(SigmaBase[i],0,0,m,m);
					_BMD.axpyin(Discr,Sigmaview,PowerSerie[N-i]);
				}
										
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
					
					 			
				// Computation of the permutation BPerm1 such that BPerm1.order is in increasing order.
				// order=Perm.order			   			
				std::vector<size_t> Perm1(m+n);			
				for (size_t i=0;i<m+n;++i)
					Perm1[i]=i;	
				if (N>=1) {
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
					

				// Discrepancy= BPerm1.Discrepancy						
				_BMD.mulin_right(BPerm1,Discrepancy);


				// Computation of the LQUP decomposition of the discrepancy
				Coefficient CopyDiscr;
				CopyDiscr=Discrepancy;
				LQUPMatrix<Field> LQUP(_F, CopyDiscr);

				// Get the matrix L of LQUP decomposition
				TriangularBlasMatrix<Element> L(m+n,m+n, BlasTag::low, BlasTag::unit );
				LQUP.getL(L);
				
				// Get the tranposed  permutation of Q from LQUP
				BlasPermutation Qt=LQUP.getQ();
			
				// Computation of permutations BPerm2 such that the last n rows of BPerm2.Qt.Discrepancy are non zero.
				std::vector<size_t> Perm2(m+n);	
				for (size_t i=0;i<n;++i)
					Perm2[i]=m+i;
				for (size_t i=n;i<m+n;++i)
					Perm2[i]=i;					
				BlasPermutation BPerm2(Perm2);			
				
				// compute the inverse of L
				TriangularBlasMatrix<Element> invL (m+n,m+n, BlasTag::low,BlasTag::unit); 
				FFPACK::trinv_left(_F,m+n,L.getPointer(),L.getStride(),invL.getWritePointer(),invL.getStride());

				
				// SigmaBase =  BPerm2.Qt. L^(-1) . BPerm1 . SigmaBase
				for (size_t i=0;i<SigmaBase.size();i++) {
					_BMD.mulin_right(BPerm1,SigmaBase[i]);
					_BMD.mulin_right(invL,SigmaBase[i]);
					_BMD.mulin_right(Qt,SigmaBase[i]);
					_BMD.mulin_right(BPerm2,SigmaBase[i]);
				}							
		
				
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


				// Discrepancy= BPerm2.U.P from LQUP
				Coefficient U(m+n,n);
				TriangularBlasMatrix<Element> trU(U,BlasTag::up,BlasTag::nonunit);
				LQUP.getU(trU);	 
				Discrepancy=U;
				BlasPermutation P= LQUP.getP();
				_BMD.mulin_left(Discrepancy,P);
				_BMD.mulin_right(BPerm2,Discrepancy);
				
				// save the first pade matrix
				if (N == degree1 -1) {
					// extract the first m rows of SigmaBase
					long max=order[0];
					for (size_t i=1;i<m;i++) {
						if (order[i]>max)
							max=order[i];
					}
					
					const Coefficient AZero(m,m);
					Approx1.resize(max+1, AZero);		
					
					for (size_t i=0;i<m;i++) 
						for (long j=0;j<=order[i];j++) 
							for (size_t k=0;k<m;k++) 
								_F.assign(Approx1[j].refEntry(i,k), SigmaBase[j].getEntry(i,k));
				}
			
			}

			if ( early_stop == 20)
				std::cout<<"Early termination is used: stop at "<<N<<" from "<<length<<" iterations\n\n";
			
			// extract the first m rows of SigmaBase
			degree=order;
			long max=degree[0];
			for (size_t i=1;i<m;i++) {
				if (degree[i]>max)
					max=degree[i];
			}
			
			const Coefficient AZero(m,m);
			Approx2.resize(max+1, AZero);		
			
			for (size_t i=0;i<m;i++) 
				for (long j=0;j<=degree[i];j++) 
					for (size_t k=0;k<m;k++) 
						_F.assign(Approx2[j].refEntry(i,k), SigmaBase[j].getEntry(i,k));				
		}


	}; // end of class SigmaBasis
	






} // end of namespace LinBox

#endif

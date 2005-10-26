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
				      size_t                           degree1, 
				      std::vector<size_t>             &defect1,
				      std::vector<Coefficient>      &SigmaBase2,
				      size_t                           degree2, 
				      std::vector<size_t>             &defect2) {
		
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
						
			// Update Serie to compute Sigma base up to degree1 - degree2			
			std::vector<Coefficient> Serie2(degree2-degree1+1,ZeroSerie);	
			ComputeNewSerieClassic(Serie2,SigmaBase1, _Serie, degree1, degree2-degree1);
		
			// copy defect
			defect2 = defect1;
			
			// Compute Sigma Base up to degree2 
			std::vector<Coefficient> Sigma2(degree2-degree1+1, Zero);
			PM_Basis(Sigma2, Serie2, degree2-degree1, defect2);	

			// multily SigmaBase1 by Sigma
			MulSigmaBasis(SigmaBase2,Sigma2,SigmaBase1);		
		}


		void multi_right_basis(std::vector<Coefficient>      &SigmaBase1,
				       size_t                           degree1, 
				       std::vector<size_t>             &defect1,
				       std::vector<Coefficient>      &SigmaBase2,
				       size_t                           degree2, 
				       std::vector<size_t>             &defect2) {
			
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

			// Update Serie to compute Sigma base up to degree1 - degree2			
			std::vector<Coefficient> Serie2(degree2-degree1+1,TransposedZero);	
			ComputeNewSerieClassic(Serie2,SigmaBase1, TransposedSerie, degree1, degree2-degree1);
		

			// copy defect
			defect2 = defect1;
			
			// Compute Sigma Base up to degree2 
			std::vector<Coefficient> Sigma2(degree2-degree1+1, Zero);
			PM_Basis(Sigma2, Serie2, degree2-degree1, defect2);	

			// multily SigmaBase1 by Sigma
			MulSigmaBasis(SigmaBase2,Sigma2,SigmaBase1);	

		}


	protected:

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

				//std::cout<<"Base"<<k<<":=";
				//write_maple(F,SigmaBase);
			}
			//std::cout<<"defect: ";
			//for (size_t i=0;i<m;++i)
			//	std::cout<<defect[i]<<" ";
			//std::cout<<std::endl;
			
			//std::cout<<"SigmaBase"<<length<<":=";
			//write_maple(F,SigmaBase);

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
		

		void ComputeNewSerieClassic(std::vector<Coefficient>          &NewSerie, 
					    const std::vector<Coefficient>   &SigmaBase, 
					    const std::vector<Coefficient>    &OldSerie,
					    size_t                              degree1,
					    size_t                              degree2){						
						
			for (size_t i=0;i<degree2+1;++i)
				for (size_t j=0;j<degree1+1;++j)
					_BMD.axpyin(NewSerie[i],SigmaBase[j],OldSerie[degree1+i-j]);
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
					{_BMD.mul(C[shiftC],A[shiftA],B[shiftB]); num_multi++;}
				else 
					if (degA == 1) 
						for (size_t i=0;i< degB;++i)
							{_BMD.mul(C[shiftC+i],A[shiftA],B[shiftB+i]);num_multi++;}
					else 
						for (size_t i=0;i< degA;++i)
							{_BMD.mul(C[shiftC+i],A[shiftA+i],B[shiftB]);num_multi++;}
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
			
			//FFTMulDomain<Field, std::vector<Coefficient> > fft_domain(_F);
			//fft_domain.mul(C,A,B);

			//for (size_t i=0;i<A.size();++i)
			//	for (size_t j=0;j<B.size();++j)
			//		_BMD.axpyin(C[i+j],A[i],B[j]);
		}	


		
	};
	

} // end of namespace LinBox

#endif

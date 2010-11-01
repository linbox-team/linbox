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

#ifndef __LINBOX_sigma_basis_H
#define __LINBOX_sigma_basis_H

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
#include <linbox/algorithms/echelon-form.h>
#include <linbox/vector/subvector.h>
#include <linbox/util/timer.h>


//#define OPTMIZED_SIGMA_UPDATE 


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
		PolynomialMatrixDomain<Field > PM_domain;
		
#ifdef _BM_TIMING
		mutable Timer ttMBasis              , tMBasis,
			ttUpdateSerie         , tUpdateSerie,
			ttBasisMultiplication , tBasisMultiplication,
			ttPermutation         , tPermutation,
			ttTransformation      , tTransformation,
			ttSigmaUp             , tSigmaUp,
			ttResidueUp           , tResidueUp,
			ttSigmaSh             , tSigmaSh,
			ttResidueSh           , tResidueSh,
			Total;

		void print(Timer& T, const  char* timer, const char* title) 
		{
			if (&T != &Total)
				Total+=T;
			if (T.count() > 0) {
				std::cout<<title<<": "<<timer;
				for (int i=strlen(timer); i<28; i++) 
					std::cout << ' ';
				double tt = (T.usertime()<0.01)? 0.01 : T.usertime();
				std::cout<<tt<<"s. (cpu)  ["<<T.count()<<"]"<<std::endl;
			}
		}
	public:
		void printTimer() 
		{
		
			print(ttMBasis              , "MBasis computation", "recursive");
			print(ttUpdateSerie         , "Updating Power Serie", "recursive");
			print(ttBasisMultiplication , "Basis Multiplication", "recursive");
			print(ttPermutation         , "Permutation"   , "iterative"); 
			print(ttTransformation      , "Transformation", "iterative");
			print(ttSigmaUp             , "Sigma Update"  , "iterative");        			
			print(ttResidueUp           , "Residue Update", "iterative");     
			print(ttSigmaSh             , "Sigma Shifting", "iterative");
			print(ttResidueSh           , "Residue Shifting","iterative");
			print(Total, "Total", "");
			std::cout<<std::endl<<std::endl;
		}

		void clearTimer() 
		{
			 ttMBasis.clear();
			 ttUpdateSerie.clear();
			 ttBasisMultiplication.clear();
			 ttPermutation.clear();   
			 ttTransformation.clear(); 
			 ttSigmaUp.clear();     
			 ttResidueUp.clear();     
			 ttSigmaSh.clear();        
			 ttResidueSh.clear(); 	 
			 Total.clear();
		}	
#endif 


	public:

		SigmaBasis(const Field &F, std::vector<Coefficient> &PowerSerie) : _F(F), _BMD(F), _MD(F), _Serie(PowerSerie), PM_domain(F)
		{
#ifdef  _BM_TIMING
			clearTimer();
#endif
		}

		void right_basis(std::vector<Coefficient>     &SigmaBase,
				size_t                           degree, 
				std::vector<size_t>             &defect) 
		{

			
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
				 std::vector<size_t>             &defect) 
		{

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
				      std::vector<size_t>              &defect2) 
		{
			
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
				       std::vector<size_t>              &defect2) 
		{
			
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
				      std::vector<size_t>                 &defect) 
		{ PadeApproximant(Approx, _Serie, degree, defect); }
	
		// function to compute the right denominator from Matrix Pade Approximant
		// compute Q(x) in S(x).Q(x) - R(x) = O(x^degree).
		
		void right_PadeMatrix (std::vector<Coefficient>            &Approx, 
				       size_t                               degree, 
				       std::vector<size_t>                 &defect) 
		{

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
					    std::vector<size_t>                  &defect) 
		{
			
			MultiPadeApproximant(Approx1, degree1, Approx2, degree2, _Serie, defect);
		}
		
		// function to compute two right denominators from Matrix Pade Approximant at two different degrees
		// compute Q1(x), Q2(x) such that  S(x).Q1(x) - R1(x) = O(x^degree1) and S(x).Q2(x) - R2(x) = O(x^degree2) 
		
		void multi_right_PadeMatrix (std::vector<Coefficient>            &Approx1, 			
					     size_t                               degree1,
					     std::vector<Coefficient>            &Approx2, 			
					     size_t                               degree2,
					     std::vector<size_t>                  &defect) 
		{

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
		// PowerSerie must have at least degree+1 element
#define MBASIS_THRESHOLD 16		
	 
		template <class Polynomial1, class Polynomial2>
		void PM_Basis(Polynomial1      &SigmaBase,
			      const Polynomial2     &PowerSerie, 
			      size_t               degree, 
			      std::vector<size_t> &defect) 
		{
						
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
				
				if (degree <= MBASIS_THRESHOLD) {
#ifdef _BM_TIMING				
					tMBasis.clear();tMBasis.start();
#endif 
					M_Basis(SigmaBase, PowerSerie, degree, defect);	 
#ifdef _BM_TIMING
					tMBasis.stop();	ttMBasis += tMBasis;
#endif			
				}
				
				else {
					size_t degree1,degree2,shift;
					shift=(degree & 0x1);
					degree1 = (degree >> 1) + shift;
					degree2 = degree - degree1;
														
					// Compute Sigma Base of half degree
					std::vector<Coefficient> Sigma1(degree1+1,ZeroSigma);
					
					std::vector<Coefficient> Serie1(degree1+1,ZeroSerie);
					for (size_t i=0;i< degree1+1;++i)
						Serie1[i] = PowerSerie[i];

					//Subvector<typename Polynomial2::iterator> Serie1(PowerSerie.begin(),PowerSerie.begin()+degree1);					
					PM_Basis(Sigma1, Serie1, degree1, defect);
					// !!! NEED TO RESIZE SIGMA1 
					// because MBasis remove all 0 matrix from leading coefficient of SigmaBase
					// while PM_Basis does not
					Sigma1.resize(degree1+1,ZeroSigma);
#ifdef _BM_TIMING			
					tUpdateSerie.clear();tUpdateSerie.start();
#endif
					// Compute Serie2 = x^(-degree1).Sigma.PowerSerie mod x^degree2
					// degree1 instead degree2 for using middle product computation
					std::vector<Coefficient> Serie2(degree1+1,ZeroSerie);	
				
					//if (2*Sigma1.size() !=  PowerSerie.size()+1)
					//std::cerr<<Sigma1.size()<<" "<<PowerSerie.size()<<" marche pas\n";
					
					
					//PM_domain.midproduct(Serie2,Sigma1,PowerSerie);
					ComputeNewSerie(Serie2,Sigma1,PowerSerie, degree1, degree2);
					Serie2.resize(degree2+1);
					
#ifdef _BM_TIMING				
					tUpdateSerie.stop();ttUpdateSerie += tUpdateSerie;
#endif					
					// Compute Sigma Base of half degree from updated Power Serie
					std::vector<Coefficient> Sigma2(degree2+1,ZeroSigma);
					
					PM_Basis(Sigma2, Serie2, degree2, defect);
						
#ifdef _BM_TIMING				
					tBasisMultiplication.clear();tBasisMultiplication.start();
#endif
					// Compute the whole Sigma Base: SigmaBase= Sigma1 x Sigma2	
					PM_domain.mul(SigmaBase,Sigma2,Sigma1);						
#ifdef _BM_TIMING				
					tBasisMultiplication.stop();ttBasisMultiplication += tBasisMultiplication;
#endif
 				}
			}
		}


		void print_multime() 
		{std::cout<<"multime: "<<PM_domain.multime<<std::endl;}

		// Computation of a minimal Sigma Base of a Power Serie up to length
		// algorithm is from Giorgi, Jeannerod and Villard  ISSAC'03		
		template <class Polynomial1, class Polynomial2>
		void M_Basis(Polynomial1        &SigmaBase,
			     Polynomial2        PowerSerie, 
			     size_t                 length, 
			     std::vector<size_t>   &defect) 
		{

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
			// I adjust the degree with the maximal difference between defects
			// this is just to be sure to catch degree increase according to elimination process
			size_t min_defect, max_defect;
			min_defect = max_defect = defect[0];
			for (size_t i=0;i<m;++i){
				if ( defect[i] > max_defect)
					max_defect = defect[i];
				if ( defect[i] < min_defect)
					min_defect = defect[i];
			}
			std::vector<size_t> degree(m,max_defect-min_defect);
		
			
			// Discrepancy
			Coefficient Discrepancy(m,n);
			Timer chrono;
#ifdef  _BM_TIMING		
                        int cptr=0;
#endif

			// Compute the minimal Sigma Base of the PowerSerie up to length
			for (size_t k=0; k< length; ++k) {

#ifdef  _BM_TIMING		
				chrono.start();
#endif
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
				
				// permute row degree
				for (size_t i=0;i<m;++i)
					std::swap(degree[i], degree[Perm1[i]]);	
	
		
#ifdef  _BM_TIMING
				chrono.stop();
				ttPermutation+=chrono;
				chrono.clear();
				chrono.start();
#endif
				// Apply Bperm1 to the current SigmaBase
				for (size_t i=0;i<SigmaBase.size();++i)
					_BMD.mulin_right(BPerm1,SigmaBase[i]);
		
#ifdef  _BM_TIMING
				chrono.stop();
				ttSigmaUp+=chrono;
				chrono.clear();
				chrono.start();
#endif				
				// Compute Discrepancy			
				_BMD.mul(Discrepancy,SigmaBase[0],PowerSerie[k]);
				for (size_t i=1;i<SigmaBase.size();i++){
					_BMD.axpyin(Discrepancy,SigmaBase[i],PowerSerie[k-i]);
				}
#ifdef  _BM_TIMING
				cptr+=SigmaBase.size();
				chrono.stop();
				ttResidueUp+=chrono;
				chrono.clear();
				chrono.start();
#endif

				//std::cout<<"MBasis: Discrepancy\n";  
				//Discrepancy.write(std::cout,_F);
			

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

#ifdef  _BM_TIMING
				chrono.stop();
				ttTransformation+=chrono;
				chrono.clear();
				chrono.start();
#endif
				// Update Sigma by L^(-1)
				// Sigma = L^(-1) . Sigma
				for (size_t i=0;i<SigmaBase.size();++i) 
					_BMD.mulin_right(invL,SigmaBase[i]);
			
#ifdef  _BM_TIMING
				chrono.stop();
				ttSigmaUp+=chrono;
				chrono.clear();
				chrono.start();
#endif
				//std::cout<<"BaseBis"<<k<<":=";
				//write_maple(F,SigmaBase);
				// Increase  degree and defect according to row choosen as pivot in LQUP
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
				if (SigmaBase.size()<= max_degree)
					{					
						SigmaBase.resize(size+1,Zero);					
						size++;
					}				
				// Mulitply by x the rows of Sigma involved as pivot in LQUP
				for (size_t i=0;i<n;++i){
					for (int j= (int) size-2;j>=0; --j){						
						for (size_t l=0;l<m;++l)
//BB: #warning Q[i] pour i>r ne veut rien dire...
							_F.assign(SigmaBase[j+1].refEntry(*(Qt.getPointer()+i),l),
								 SigmaBase[j].getEntry(*(Qt.getPointer()+i),l));			
					}
					for (size_t l=0;l<m;++l)
						_F.assign(SigmaBase[0].refEntry(*(Qt.getPointer()+i),l),zero);
				}
#ifdef  _BM_TIMING
				chrono.stop();
				ttSigmaSh+=chrono;
				chrono.clear();
				chrono.start();
#endif
				//write_maple("SS1",SigmaBase);			
			}
		}


		// Multiply a Power Serie by a Sigma Base.
		// only affect coefficients of the Power Serie between degree1 and degree1+degree2
		template<class Polynomial1, class Polynomial2,class Polynomial3>
		inline void ComputeNewSerie(Polynomial1          &NewSerie, 
					    const Polynomial2   &SigmaBase, 
					    const Polynomial3    &OldSerie,
					    size_t                 , //deg1
					    size_t                 )//deg2					
		{
			
			// degree1 >= degree2
			//size_t size = 2*degree1 + 1;
					
			const Coefficient ZeroSerie (OldSerie[0].rowdim(), OldSerie[0].coldim());
			//const Coefficient ZeroBase  (SigmaBase[0].rowdim(), SigmaBase[0].coldim());
			
			// Work on a copy of the old  Serie (increase size by one for computation of middle product)
			
			std::vector<Coefficient> Serie(OldSerie.size()+1, ZeroSerie);
			for (size_t i=0;i< OldSerie.size();++i)
				Serie[i] = OldSerie[i];			
			Serie[OldSerie.size()]=ZeroSerie;
			
			//  ** try to not use a Copy **
			// Work on a copy of the Sigma Base 
			//std::vector<Coefficient> Sigma(SigmaBase.size());
			//for (size_t i=0;i<SigmaBase.size();++i){
			//	Sigma[i] = SigmaBase[i];
			//}
		

			PM_domain.midproduct(NewSerie, SigmaBase, Serie);
		}
		

		void ComputeNewSerieClassic(std::vector<Coefficient>          &NewSerie, 
					    const std::vector<Coefficient>   &SigmaBase, 
					    const std::vector<Coefficient>    &OldSerie,
					    size_t                              degree1,
					    size_t                              degree2)
		{
						
			for (size_t i=0;i<degree2+1;++i)
				for (size_t j=0;j<degree1+1;++j)
					_BMD.axpyin(NewSerie[i],SigmaBase[j],OldSerie[degree1+i-j]);
		}


		


		// Multiply two Sigma Basis
		void MulSigmaBasis(std::vector<Coefficient> &C, 
				   std::vector<Coefficient> &A,
				   std::vector<Coefficient> &B)
		{PM_domain.mul(C,A,B);}
		
		
		void PadeApproximant (std::vector<Coefficient>            &Approx, 
				      const std::vector<Coefficient>  &PowerSerie, 
				      size_t                               length, 
				      std::vector<size_t>                 &defect) 
		{//pph
			

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
					   std::vector<size_t>                  &defect) 
		{
					   

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


		// Computation of a minimal Sigma Base of a Power Serie up to length
		// algorithm is from Giorgi, Jeannerod and Villard  ISSAC'03		
		void new_M_Basis(std::vector<Coefficient>     &SigmaBase,
				 std::vector<Coefficient>    &PowerSerie, 
				 size_t                           length, 
				 std::vector<size_t>             &defect) 
		{

			// Get the dimension of matrices inside 
			// the Matrix Power Serie
			size_t m,n;
			m = PowerSerie[0].rowdim();
			n = PowerSerie[0].coldim();

			// Set some useful constants
			const Coefficient Zeromm(m,m);
			const Coefficient Zeromn(m,n);
			Element one, zero;
			_F.init(one,1UL);
			_F.init(zero,0UL);
			
			// Reserve memory for the Sigma Base  
			SigmaBase.reserve(length+1);
			SigmaBase.resize(1);
			
			// set SigmaBase[0] to Identity
			Coefficient Identity(m,m);
			for (size_t i=0;i< m;++i)
				Identity.setEntry(i,i,one);
			SigmaBase[0]=Identity;
			
			// Define Truncated Residual
			std::vector<Coefficient>  Residual(length+1, Zeromn);
			
			// Set Truncated Residual to PowerSerie mod X^length
			for (size_t k=0;k<length; ++k)
				Residual[k] = PowerSerie[k];

			// Define Discrepancy 
			Coefficient Discrepancy(m,n);
			
			// Row degree of SigmaBase
			// Keep track on Sigma Base's row degree
			// I adjust the degree with the maximal difference between defects
			// this is just to be sure to catch degree increase according to elimination process
			int min_defect, max_defect;
			min_defect = max_defect = defect[0];
			for (size_t i=0;i<m;++i){
				if ( defect[i] > max_defect)
					max_defect = defect[i];
				if ( defect[i] < min_defect)
					min_defect = defect[i];
			}
			//std::vector<size_t> degree(m,max_defect-min_defect);
			std::vector<size_t> degree(m);
			//size_t max_diff= max_defect-min_defect;
			for (size_t i=0;i<m;++i)
				degree[i]= defect[i]-min_defect;


			std::vector<size_t> row_degree(m,0);

			//write_maple("PowerSerie",PowerSerie);		

#ifdef _BM_TIMING
			Timer chrono;
			double tSigmaUp, tResidueUp, tSigmaSh, tResidueSh, tLQUP, tPerm;
			tSigmaUp= tResidueUp= tSigmaSh= tResidueSh= tLQUP= tPerm =0.;
#endif

			std::vector<size_t> triv_column(m,0);
			std::vector<size_t> PermPivots(m);
			for (size_t i=0;i<m;++i)
				PermPivots[i]=i;

			//int optim=0;
			int cptr=0;
#ifdef DISCREPANCY_OPTIM
			double updis=0.;
#endif

			//bool UseTrivial=true;

			//std::cout<<"defect: ";
			//for (size_t i=0;i<m;++i)
			//	std::cout<<defect[i]<<", ";
			//std::cout<<"\n";

			// Compute the minimal Sigma Base of the PowerSerie up to length
			for (size_t k=0; k< length; ++k) {
		
				/*
				std::cout<<"row degree: ";
				for (size_t i=0;i<m;++i)
					std::cout<<degree[i]<<", ";
				std::cout<<"\n";
				write_maple("Sigma",SigmaBase);
				*/
				
#ifdef _BM_TIMING
				chrono.start();
#endif
				// Compute the number of trivial column in SigmaBase
				int nbr_triv=0;			       
				for (size_t i=0;i<m;i++) if (triv_column[i]==0) nbr_triv++;


#ifdef DISCREPANCY_OPTIM
				// Compute a Permutation to put all trivial columns of SigmaBase to right part of it
				std::vector<size_t> PermTrivial(m);
				bool PTrivial=true;
				if (nbr_triv > 0) {
					size_t idx_triv, idx_nontriv;
					idx_nontriv = 0; idx_triv = m-nbr_triv;
										
					for (size_t i=0;i<m;++i){
						if (triv_column[i]!=0){
							PermTrivial[i]=idx_nontriv;
							if (i!=idx_nontriv) PTrivial=false;
							idx_nontriv++;
						}
						else {						
							PermTrivial[i]=idx_triv;
							if (i!=idx_triv) PTrivial=false;
							idx_triv++;
						}
					}
				}

				BlasPermutation PPP (PermTrivial);
				TransposedBlasMatrix<BlasPermutation> PPPT(PPP);
				
				// set Discrepancy to Residual[k] 
				// -> can be optimized by directly using Residual[k]
				//Discrepancy = Residual[k];					


				BlasMatrix<Element> Db    (Discrepancy,m-nbr_triv,0,nbr_triv,n);
				// Compute Discrepancy using convolution
				if (nbr_triv > m){ 					
					if (k==0){
						_MD.copy(Discrepancy, PowerSerie[0]);
					}
					else {	
						BlasPermutation PPiv (PermPivots);
						TransposedBlasMatrix<BlasPermutation> PPivT(PPiv);
				
						if (!PTrivial){
							_BMD.mulin_left(SigmaBase[0], PPP);
							_BMD.mulin_right(PPPT, PowerSerie[k]);					
						}
						_BMD.mulin_right(PPivT,SigmaBase[0]);
						
						BlasMatrix<Element> SBl   (SigmaBase[0],0,0,m,m-nbr_triv);
						BlasMatrix<Element> PSt   (PowerSerie[k],0,0,m-nbr_triv,n);
						BlasMatrix<Element> PSb   (PowerSerie[k],m-nbr_triv,0,nbr_triv,n);
					
						Timer Dchrono;
						Dchrono.start();
						_BMD.mul (Discrepancy,SBl,PSt);
						if (k%2 == 0)
							_MD.addin(Db,PSb);							
						else
							_MD.subin(Db,PSb);																							
						Dchrono.stop();
						updis+= Dchrono.usertime();
						
						_BMD.mulin_right(PPiv,SigmaBase[0]);				
						if (!PTrivial){
							_BMD.mulin_right(PPP, PowerSerie[k]);
							_BMD.mulin_left(SigmaBase[0], PPPT);
						}
						
						for (size_t i=1;i<SigmaBase.size();i++){
							if (!PTrivial){
								_BMD.mulin_left(SigmaBase[i], PPP);
								_BMD.mulin_right(PPPT, PowerSerie[k-i]);
							}
							_BMD.mulin_right(PPivT,SigmaBase[i]);

							BlasMatrix<Element> SBli   (SigmaBase[i],0,0,m,m-nbr_triv);
							BlasMatrix<Element> PSti   (PowerSerie[k-i],0,0,m-nbr_triv,n);
							BlasMatrix<Element> PSbi   (PowerSerie[k-i],m-nbr_triv,0,nbr_triv,n);
							Dchrono.clear();
							Dchrono.start();					       
							_BMD.axpyin  (Discrepancy,SBli,PSti);
							Dchrono.stop();
							updis+= Dchrono.usertime();
							_BMD.mulin_right(PPiv,SigmaBase[i]);
							if (!PTrivial){
								_BMD.mulin_right(PPP, PowerSerie[k-i]);
								_BMD.mulin_left(SigmaBase[i], PPPT);											
							}
						}
						_BMD.mulin_right(PPiv,Discrepancy);
					}
				}
				else{
#endif			
					_BMD.mul(Discrepancy,SigmaBase[0],PowerSerie[k]);
					for (size_t i=1;i<SigmaBase.size();i++){
						_BMD.axpyin(Discrepancy,SigmaBase[i],PowerSerie[k-i]);
					}
#ifdef DISCREPANCY_OPTIM
				}
#endif
				

#ifdef _BM_TIMING
				chrono.stop();
				ttResidueUp+=chrono;
				chrono.clear();
				chrono.start();
#endif			
			
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

				// permute degree and row degree
				for (size_t i=0;i<m;++i)
					if ( i < Perm1[i]){
						std::swap(degree[i], degree[Perm1[i]]);	
						std::swap(row_degree[i], row_degree[Perm1[i]]);	
						if (k!=0)
							std::swap(PermPivots[i], PermPivots[Perm1[i]]);
					}

				if (k==0)
					for (size_t i=0;i<m;++i)
						PermPivots[i]=Perm1[i];									

				BlasPermutation BPerm1(Perm1);

				// Apply Bperm1 to the Discrepancy
				_BMD.mulin_right(BPerm1, Discrepancy);	
							
#ifdef _BM_TIMING		
				chrono.stop();
				ttPermutation+=chrono;
				chrono.clear();
				chrono.start();
#endif

				/* new version : use of columnReducedEchelon */							
				EchelonFormDomain<Field>  EFD(_F);			
				size_t rank = EFD.columnReducedEchelon(Discrepancy);
			
			
				// compute permutation such that all pivots are on the principal minor
				std::vector<size_t> perm(m);
				for (size_t i=0;i<m;++i)
					perm[i]=i;				
				size_t idx=0;
				for (size_t i=0;i<rank;++i){
					while(_F.isZero(Discrepancy.getEntry(idx,i))) idx++;
					perm[i]=idx;
					idx++;
				}
			
				BlasPermutation Qt(perm);						
				TransposedBlasMatrix<BlasPermutation> Q(Qt);

				// detect trivial permutation (Q=Qt=Identity)
				bool QisTrivial=true;
				size_t Qidx=0;
				while(QisTrivial && (Qidx <m)){
					if (perm[Qidx] != Qidx)
						QisTrivial=false;
					Qidx++;
				}				

				// put all pivots on the principal minor
				_BMD.mulin_right(Qt, Discrepancy);

				// Get the (m-r)*r left bottom submatrix of Reduced Echelon matrix 
				BlasMatrix<Element> G(Discrepancy, rank, 0,m-rank,rank);
#ifdef _BM_TIMING													
				chrono.stop();
				ttTransformation+=chrono;
				chrono.clear();
				chrono.start();
#endif
	
			
								    				
				// compute size of trivial part of SigmaBase
				size_t rsize, lsize;
				rsize=0;
				if (nbr_triv>rank)
					rsize=nbr_triv-rank;
				lsize=m-rsize;

				//std::cout<<"rsize: "<<rsize<<"\n";

				// compute maximal degree of first rank-th row of SigmaBase
				
				size_t maxs=0;  
				for(size_t d=0;d<rank;++d)
					maxs=std::max(maxs, degree[*(Qt.getPointer()+d)]);				
				maxs=std::min(maxs, SigmaBase.size()-1);
				
				//Discrepancy.write(std::cout,_F);

			
#ifndef OPTMIZED_SIGMA_UPDATE				
				// Compute a Permutation to put all trivial columns of SigmaBase to right part of it
				std::vector<size_t> PermTrivial(m);
				bool PTrivial=true;
				for(size_t i=0;i<m;++i)
					if (PermPivots[i]>i)
						std::swap(triv_column[i], triv_column[PermPivots[i]]);

				if (nbr_triv > rank) {
					size_t idx_triv, idx_nontriv;
					idx_nontriv = 0; idx_triv = m-nbr_triv;
										
					for (size_t i=0;i<m;++i){
						if (triv_column[i]!=0){
							PermTrivial[i]=idx_nontriv;
							if (i!=idx_nontriv) PTrivial=false;
							idx_nontriv++;
						}
						else {						
							PermTrivial[i]=idx_triv;
							if (i!=idx_triv) PTrivial=false;
							idx_triv++;
						}
					}
				}
				for(int i=m-1;i>=0;--i)
					if (PermPivots[i]>i)
						std::swap(triv_column[i], triv_column[PermPivots[i]]);
				
				// Modify Permutation of trivial columns to incorporate pivot columns
				if (nbr_triv>rank){
					for(size_t i=0;i<m;++i)
						std::swap(PermTrivial[i], PermTrivial[*(Qt.getPointer()+i)]);								
				}
			
				/*
				std::cout<<"PermTrivial: ";
				for (size_t i=0;i<m;++i)
					std::cout<<PermTrivial[i]<<", ";
				std::cout<<"\n";
				*/

				BlasPermutation PTr (PermTrivial);
				TransposedBlasMatrix<BlasPermutation> PTrT(PTr);
				//BlasPermutation PPP (PermTrivial);
				//TransposedBlasMatrix<BlasPermutation> PPPT(PPP);
				BlasPermutation PPiv (PermPivots);
				TransposedBlasMatrix<BlasPermutation> PPivT(PPiv);


				/*
				std::cout<<"MAXS: "<<maxs<<"\n\n ******************\n";
				
				std::cout<<"Perm1: ";
				for (size_t i=0;i<m;++i)
					std::cout<<Perm1[i]<<", ";
				std::cout<<"\n";
				
				write_maple("Sigma",SigmaBase);
				//Discrepancy.write(std::cout,_F);

			

				std::cout<<"trivial: ";
				for (size_t i=0;i<m;++i)
					std::cout<<triv_column[i]<<", ";
				std::cout<<"\n";
				
				SigmaBase[0].write(std::cout,_F);
				_BMD.mulin_right(BPerm1, SigmaBase[0]);			
				SigmaBase[0].write(std::cout,_F);			
				_BMD.mulin_left(SigmaBase[0],PPiv);
				SigmaBase[0].write(std::cout,_F);
				_BMD.mulin_left(SigmaBase[0], PTr);				
				SigmaBase[0].write(std::cout,_F);
				_BMD.mulin_left(SigmaBase[0], PTrT);
				_BMD.mulin_left(SigmaBase[0],PPivT);
				TransposedBlasMatrix<BlasPermutation> BPerm1T(BPerm1);
				_BMD.mulin_right(BPerm1T, SigmaBase[0]);
				*/

				// Update SigmaBase 
				for (size_t i=0;i<maxs+1;++i) {	
					//SigmaBase[0].write(std::cout,_F);
					// permute SigmaBase
					_BMD.mulin_right(BPerm1, SigmaBase[i]);
													
					if (!QisTrivial){
						_BMD.mulin_right(Qt, SigmaBase[i]);					
					}
																				
					if (nbr_triv > rank){
						_BMD.mulin_left(SigmaBase[i],PPivT);
						_BMD.mulin_left(SigmaBase[i], PTr);
					}
												
					
					// apply transformation to SigmaBase
				 	//BlasMatrix<Element>    S_top(SigmaBase[i], 0,0,rank,m);
				 	//BlasMatrix<Element> S_bottom(SigmaBase[i], rank,0,m-rank,m);	
					//_BMD.axmyin(S_bottom, G, S_top); 			
											
					BlasMatrix<Element> S_top_left    (SigmaBase[i], 0,0,rank,lsize);
					BlasMatrix<Element> S_bottom_left (SigmaBase[i], rank,0,m-rank,lsize);
					//if (i==0){						
					//	S_bottom_left.write(std::cout,_F);
					//}
					// deal with the left part of S_bottom
					_BMD.axmyin(S_bottom_left, G, S_top_left);
					//if (i==0){						
					//	S_bottom_left.write(std::cout,_F);
					//}
				
					
					// deal with the right part of S_bottom
					if (rsize > 0){	
						BlasMatrix<Element> S_bottom_right(SigmaBase[i],rank,lsize,m-rank, rsize);
						_MD.negin(S_bottom_right);
					}
					
					// undo the permutation on sigmaBase
					if (nbr_triv > rank){
						_BMD.mulin_left(SigmaBase[i], PTrT);
						_BMD.mulin_left(SigmaBase[i],PPiv);
					}

					if (!QisTrivial)
					 _BMD.mulin_right(Q, SigmaBase[i]);
					//if (i==0){						
					//	S_bottom_left.write(std::cout,_F);
					//}					
				} 		
				
				for (size_t i=maxs+1;i<SigmaBase.size();++i) {
					// permute SigmaBase
					_BMD.mulin_right(BPerm1, SigmaBase[i]);
					if (!QisTrivial)
						_BMD.mulin_right(Qt, SigmaBase[i]);					
					
					// apply transformation to SigmaBase
					BlasMatrix<Element> S_bottom(SigmaBase[i], rank,0,m-rank,m);				
				 	_MD.negin(S_bottom);

					// undo the permutation on sigmaBase
					if (!QisTrivial)
					 _BMD.mulin_right(Q, SigmaBase[i]);					
				} 				
		
#else		  
				/*********************************/
				/* OPTIMIZATION Update SigmaBase */
				/*********************************/
				
				// Permute SigmaBase
				for (size_t i=0;i<SigmaBase.size();++i) {
					_BMD.mulin_right(BPerm1, SigmaBase[i]);
					if (!QisTrivial)
						_BMD.mulin_right(Qt, SigmaBase[i]);
					if (nbr_triv > rank)
						_BMD.mulin_left(SigmaBase[i], PTr);
				}
		 						
				// split G into Gk of dimension r x r				
				size_t q = (m-rank) / rank;
				size_t q_last= (m-rank) % rank;
							
				for (size_t l=0;l<q;++l){					
					// get kth part of G
					BlasMatrix<Element> Gk(G,l*rank,0,rank,rank);

					// get maximal degree of kth part of sigma base
					size_t maxk=0;
					for(size_t d=(l+1)*rank;d<(l+2)*rank;++d)
						maxk=std::max(maxk, degree[*(Qt.getPointer()+d)]);
					maxk=std::min(maxk,SigmaBase.size()-1);

					for (size_t i=0;i<maxs+1;++i){
						BlasMatrix<Element> S_top_left (SigmaBase[i], 0,0,rank,lsize);
						BlasMatrix<Element> S_bottom   (SigmaBase[i], rank,0,m-rank,m);
						
						// deal with the left part of kth slice of S_bottom
						BlasMatrix<Element> Sbk_left (S_bottom,l*rank,0,rank, lsize);					
						_BMD.axmyin(Sbk_left, Gk, S_top_left);	
						
						// deal with the right part of kth slice of S_bottom
						if (rsize > 0){
							BlasMatrix<Element> Sbk_right(S_bottom,l*rank,lsize,rank, rsize);
							_MD.negin(Sbk_right);
						}
					} 
					for (size_t i=maxs+1;i<maxk+1;++i){
						BlasMatrix<Element> S_bottom(SigmaBase[i], rank,0,m-rank,m);
						BlasMatrix<Element> Sbk(S_bottom,l*rank,0,rank, m);
						_MD.negin(Sbk);
					} 
				}		
				if( q_last > 0) {
					// get last part of G
					BlasMatrix<Element> G_last(G,q*rank,0,q_last,rank);

					size_t maxk=0;
					for(size_t d=m-q_last;d<m;++d)
						maxk=std::max(maxk, degree[*(Qt.getPointer()+d)]);
					maxk=std::min(maxk,SigmaBase.size()-1);
					for (size_t i=0;i<maxs+1;++i){
						BlasMatrix<Element> S_top_left (SigmaBase[i], 0,0,rank,lsize);
					 	BlasMatrix<Element> S_bottom   (SigmaBase[i], rank,0,m-rank,m);
						
						// deal with the left part of kth slice of S_bottom
					 	BlasMatrix<Element> Sb_last_left (S_bottom,q*rank,0,q_last, lsize);					
					 	_BMD.axmyin(Sb_last_left, G_last, S_top_left);	
						
						// deal with the right part of kth slice of S_bottom
					 	if (rsize > 0){
							BlasMatrix<Element> Sb_last_right(S_bottom,q*rank,lsize,q_last, rsize);
						 	_MD.negin(Sb_last_right);
						}
						
					} 
					for (size_t i=maxs+1;i<maxk+1;++i){
						BlasMatrix<Element> S_bottom(SigmaBase[i], rank,0,m-rank,m);
						BlasMatrix<Element> Sb_last(S_bottom,q*rank,0,q_last, m);
						_MD.negin(Sb_last);
					}				
				} 
						
				// undo Permutation of sigma Base
				for (size_t i=0;i<SigmaBase.size();++i) {
					if (nbr_triv > rank)
						_BMD.mulin_left(SigmaBase[i], PTrT);
					if (!QisTrivial)
						_BMD.mulin_right(Q, SigmaBase[i]);				
				}						
#endif
				// END OF OPTIMIZATION
			
	
#ifdef _BM_TIMING
				chrono.stop();
				ttSigmaUp+=chrono;
				chrono.clear();
				chrono.start();
#endif			
				/*
				// Update  Residual (only monomials greater than k-1)
				for (size_t i=k;i<length;++i){
					_BMD.mulin_right(BPerm1,Residual[i]);
					//_BMD.mulin_right(invL,Residual[i]);

					// try optimization
					_BMD.mulin_right(Qt, Residual[i]);
					 BlasMatrix<Element>    R_top(Residual[i], 0,0,rank,n);
					 BlasMatrix<Element> R_bottom(Residual[i], rank,0,m-rank,n);
					 _BMD.axmyin(R_bottom, G, R_top);	
					 				
					 _BMD.mulin_right(Q, Residual[i]);
				}	
				*/
#ifdef _BM_TIMING
				chrono.stop();
				ttResidueUp+=chrono;
				chrono.clear();
				chrono.start();
#endif			
				//  Calculate the new degree of SigmaBase (looking only pivot's row)
				size_t max_degree=degree[*(Qt.getPointer())];
				for (size_t i=1;i<rank;++i) {
					if (degree[*(Qt.getPointer()+i)]>max_degree)
						max_degree=degree[*(Qt.getPointer()+i)];
				}					

				// Resize SigmaBase if necessary
				size_t size=SigmaBase.size(); 
				if (SigmaBase.size()<= max_degree+1)
					{					
						SigmaBase.resize(size+1,Zeromm);					
						size++;cptr++;
					}		
					
				// Mulitply by x the rows of SigmaBase involved as pivot 
				for (size_t i=0;i<rank;++i){
					for (int j= (int) size-2;j>=0; --j){						
						for (size_t l=0;l<m;++l)
							_F.assign(SigmaBase[j+1].refEntry(*(Qt.getPointer()+i),l), SigmaBase[j].getEntry(*(Qt.getPointer()+i),l));			
					}
					for (size_t l=0;l<m;++l)
						_F.assign(SigmaBase[0].refEntry(*(Qt.getPointer()+i),l),zero);
				}
#ifdef _BM_TIMING
				chrono.stop();
				ttSigmaSh+=chrono;
				chrono.clear();
				chrono.start();
#endif
				/*
				// Mulitply by x the rows of Residual involved as pivot 				
				for (size_t i=0;i<rank;++i){
					for (int j= (int) length-2;j>= (int) k; j--){
						for (size_t l=0;l<n;++l) 
							_F.assign(Residual[j+1].refEntry(*(Qt.getPointer()+i),l), Residual[j].getEntry(*(Qt.getPointer()+i),l));
					}					
				}
				*/
#ifdef _BM_TIMING
				chrono.stop();
				ttResidueSh+=chrono;
				chrono.clear();
#endif
				// Increase defect according to row index choosen as pivot 				
				for (size_t i=0;i<rank;++i){
					defect[*(Qt.getPointer()+i)]++;	
					degree[*(Qt.getPointer()+i)]++;
					row_degree[*(Qt.getPointer()+i)]++;
					triv_column[PermPivots[*(Qt.getPointer()+i)]]++;					
				}											
			}
			//std::cout<<"row degree: ";
			//for (size_t i=0;i<m;++i)
			//	std::cout<<degree[i]<<", ";
			//std::cout<<"\n";

			//write_maple("Sigma",SigmaBase);
			//std::cout<<"cpt:= "<<cptr<<"\n";
			//if (length > 2) 
			//	std::cout<<"permformed :"<<optim<<" s x s matrix multiplication\n";
			//std::cout<<"updsi time: "<<updis<<"s \n";
		}


		void new_PM_Basis(std::vector<Coefficient>     &SigmaBase,
				  std::vector<Coefficient>    &PowerSerie, 
				  size_t                           degree, 
				  std::vector<size_t>             &defect) 
		{
						
			size_t m,n;
			m = PowerSerie[0].rowdim();
			n = PowerSerie[0].coldim();
			Element one;
			_F.init(one,1UL);
			const Coefficient ZeroSigma(m,m);
			const Coefficient ZeroSerie(m,n);
			PowerSerie.resize(degree, ZeroSerie);	

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
					//write_maple("nPowerSerie", PowerSerie);
					//new_M_Basis(SigmaBase, PowerSerie, degree, defect);
					M_Basis(SigmaBase, PowerSerie, degree, defect);
					//write_maple("\nSigmaBase", SigmaBase);

#ifdef _BM_TIMING
					tMBasis.stop();
					ttMBasis += tMBasis;
#endif			
				}
				else {
					size_t degree1,degree2;
					degree1 = (degree >> 1) + (degree & 1);				
					degree2 = degree - degree1;									
					
					//write_maple("\n\nPowerSerie", PowerSerie);

					// Compute Sigma Base of half degree
					std::vector<Coefficient> Sigma1(degree1,ZeroSigma);
					

					std::vector<Coefficient> Serie1(degree1);
					for (size_t i=0;i< degree1;++i)
						Serie1[i] = PowerSerie[i];
					
					//write_maple("Serie1", Serie1);
				
					new_PM_Basis(Sigma1, Serie1, degree1, defect);
					//new_PM_Basis(Sigma1, PowerSerie, degree1, defect);

					size_t S1size= Sigma1.size();

#ifdef _BM_TIMING				
					tUpdateSerie.clear();
					tUpdateSerie.start();
#endif
					/*
					//write_maple("Sigma1", Sigma1);				
					Sigma1.resize(degree1+1, ZeroSigma);
					
					// Compute Serie2 = x^(-degree1).Sigma.PowerSerie mod x^degree2
					std::vector<Coefficient> Serie2(degree1+1,ZeroSerie);										
								
					// Work on a copy of the old  Serie (increase size by one for computation of middle product)					
					std::vector<Coefficient> Serie(2*degree1+1,ZeroSerie);					
					for (size_t i=0;i< PowerSerie.size();++i)
						Serie[i] = PowerSerie[i];
					
					PM_domain.midproduct(Serie2, Sigma1, Serie);
					//ClassicMulDomain<Field, std::vector<Coefficient> > CM_domain(_F);
					//CM_domain.midproduct(Serie2, Sigma1, Serie); 
					Serie2.resize(degree2, ZeroSerie);
					*/
					std::vector<Coefficient> Serie2(degree2,ZeroSerie);
					UpdateSerie(Serie2, Sigma1, PowerSerie, degree1, degree2);
					
					
#ifdef _BM_TIMING				
					tUpdateSerie.stop();
					ttUpdateSerie += tUpdateSerie;
#endif
					//write_maple("Serie2", Serie2);
					
					// Compute Sigma Base of half degree from updated Power Serie					
					std::vector<Coefficient> Sigma2(degree2,ZeroSigma);
					new_PM_Basis(Sigma2, Serie2, degree2, defect);

					//std::cout<<"starting multiplication "<<Sigma2.size()<<" x "<<Sigma1.size()<<"...\n"; 

					//write_maple("Sigma2", Sigma2);
							
					// Compute the whole Sigma Base through the product 
					// of the Sigma Basis Sigma1 x Sigma2						
				
										
#ifdef _BM_TIMING				
					tBasisMultiplication.clear();
					tBasisMultiplication.start();
#endif	
					// Remove leading Zero coefficient of Sigma1 and Sigma2
					/*
					size_t idx1,idx2;
					idx1=Sigma1.size();
					idx2=Sigma2.size();
					while( _MD.isZero(Sigma1[idx1-1]) && idx1 >0) {idx1--;}
					while( _MD.isZero(Sigma2[idx2-1]) && idx2 >0) {idx2--;}					
				
					//std::cout<<"zero removed: "<<Sigma1.size()-idx1+Sigma2.size()-idx2<<"\n";
					// resize Sigma1 ad Sigma2
					Sigma1.resize(idx1);
					Sigma2.resize(idx2);
					*/
					
					//resize SigmaBase
					SigmaBase.resize(Sigma1.size()+Sigma2.size()-1, ZeroSigma);
					
					PM_domain.mul(SigmaBase,Sigma2,Sigma1);	

					// Remove leading Zero coefficient of SigmaBase
					size_t idx;
					idx=SigmaBase.size();
					while( _MD.isZero(SigmaBase[idx-1]) && idx >0) {idx--;}
					SigmaBase.resize(idx, ZeroSigma);
					
					
					
#ifdef _BM_TIMING				
					tBasisMultiplication.stop();
					ttBasisMultiplication += tBasisMultiplication;
#endif
					//write_maple("SigmaBase", SigmaBase);
				}
			}
		}


		// Multiply a Power Serie by a Sigma Base.
		// only affect coefficients of the Power Serie between degree1 and degree1+degree2-1
		void UpdateSerie(std::vector<Coefficient>                &NewSerie, 
				 std::vector<Coefficient>               &SigmaBase, 
				 const std::vector<Coefficient>          &OldSerie,
				 size_t                                    degree1,
				 size_t                                    degree2)
		{

			size_t m,n;
			m = OldSerie[0].rowdim();
			n = OldSerie[0].coldim();
			const Coefficient ZeroSigma(m,m);
			const Coefficient ZeroSerie(m,n);
			size_t Ssize = SigmaBase.size();

			if (SigmaBase.size() < 5){
			
				// do the calculation by hand
				for (size_t j=degree1;j<degree1+degree2;++j){
					_BMD.mul(NewSerie[j-degree1], SigmaBase[0], OldSerie[j]);
					for (size_t i=1;i<Ssize; ++i)				
						_BMD.axpyin(NewSerie[j-degree1], SigmaBase[i], OldSerie[j-i]);
				}						
			}
			else{
				//std::cout<<"Sigma size: "<<Ssize<<" -> "<<degree1<<" -> "<<degree2<<"\n";
				// resize to fit the requirement of middle product algorithm
				SigmaBase.resize(degree1+1, ZeroSigma);
				NewSerie.resize (degree1+1, ZeroSerie);				
				std::vector<Coefficient> Serie(2*degree1+1,ZeroSerie);					
				for (size_t i=0;i< OldSerie.size();++i)
					Serie[i] = OldSerie[i];
				
				// call middle product
				PM_domain.midproduct(NewSerie, SigmaBase, Serie);

				// resize results and entries
				NewSerie.resize(degree2, ZeroSerie);
				SigmaBase.resize(Ssize, ZeroSigma);
			}
		}


		void write_maple(const char* name, const Coefficient & C) 
		{
			size_t m,n;
			m = C.rowdim();
			n = C.coldim();
			std::cout<<"Fx:=proc(P) local i; return eval(sum(x^(i-1)*P[i],i=1..nops(P))); end proc;";

			std::cout<<name<<":=Fx(["; 
			std::cout<<"Matrix([";
			for (size_t i=0;i<m-1;++i){
				std::cout<<"[";
				for (size_t j=0;j<n-1;++j)
					_F.write(std::cout,C.getEntry(i,j))<<",";
				_F.write(std::cout,C.getEntry(i,n-1))<<"] , ";
			}
			std::cout<<"[";
			for (size_t j=0;j<n-1;++j)
				_F.write(std::cout,C.getEntry(m-1,j))<<",";				
			_F.write(std::cout, C.getEntry(m-1,n-1))<<"]])]); ";

		}

		void write_maple(const char* name, const std::vector<Coefficient> & P) 
		{
			size_t m,n;
			m = P[0].rowdim();
			n = P[0].coldim();
			std::cout<<"Fx:=proc(P) local i; return eval(sum(x^(i-1)*P[i],i=1..nops(P))); end proc;";

			std::cout<<name<<":=Fx(["; 
			for (size_t k=0;k<P.size()-1;++k){
				std::cout<<"Matrix([";
				for (size_t i=0;i<m-1;++i){
					std::cout<<"[";
					for (size_t j=0;j<n-1;++j)
						_F.write(std::cout,P[k].getEntry(i,j))<<",";
					_F.write(std::cout, P[k].getEntry(i,n-1))<<"] , ";
				}
				std::cout<<"[";
				for (size_t j=0;j<n-1;++j)
					_F.write(std::cout,P[k].getEntry(m-1,j))<<",";				
				_F.write(std::cout, P[k].getEntry(m-1,n-1))<<"]]) , ";	
			}
			
			std::cout<<"Matrix([";
			for (size_t i=0;i<m-1;++i){
				std::cout<<"[";
				for (size_t j=0;j<n-1;++j)
					_F.write(std::cout,P[P.size()-1].getEntry(i,j))<<",";
				_F.write(std::cout, P[P.size()-1].getEntry(i,n-1))<<"] , ";
			}
			std::cout<<"[";
			for (size_t j=0;j<n-1;++j)
				_F.write(std::cout,P[P.size()-1].getEntry(m-1,j))<<",";				
			_F.write(std::cout, P[P.size()-1].getEntry(m-1,n-1))<<"]])]); \n\n";	
		}


	}; // end of class SigmaBasis

} // end of namespace LinBox

#undef OPTMIZED_SIGMA_UPDATE

#endif //__LINBOX_sigma_basis_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

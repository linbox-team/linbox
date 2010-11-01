/* linbox/algorithms/
 * Copyright (C) 2005  Pascal Giorgi
 *
 * Written by Pascal Giorgi <pgiorgi@uwaterloo.ca>
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


#ifndef __LINBOX_matpoly_mult_H
#define __LINBOX_matpoly_mult_H

#include <linbox/randiter/random-fftprime.h>
#include <linbox/algorithms/blas-domain.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/util/error.h>
#include <linbox/util/debug.h>
#include <linbox/util/timer.h>
#include <vector>
#ifdef __LINBOX_HAVE_OPENMP
#include <omp.h>
#endif
//#define FFT_TIMING

namespace LinBox 
{


#define FFT_DEG_THRESHOLD   64
#define KARA_DEG_THRESHOLD  1
#ifndef FFT_PRIME_SEED
	// random seed
#define FFT_PRIME_SEED 0
#endif


	template <class Field>
	class KaratsubaMulDomain;

	template <class _Field>
	class FFTMulDomain;

	template <class Field>
	class ClassicMulDomain;
	
	template <class Field>
	class PolynomialMatrixDomain {
	protected:
		KaratsubaMulDomain<Field>     _kara;
		FFTMulDomain<Field>            _fft;
		ClassicMulDomain<Field>    _classic;
		
	public:
		Timer multime;
		
		PolynomialMatrixDomain ( const Field &F) : _kara(F), _fft(F), _classic(F) {multime.clear();}

		template< class Polynomial1, class Polynomial2, class Polynomial3>
		void mul (Polynomial1 &a, const Polynomial2 &b, const Polynomial3 &c) {
			size_t d = b.size()+c.size();
			linbox_check(a.size() >= (b.size()+c.size()-1));
			//Timer mul;
			//mul.start();
			if (d > FFT_DEG_THRESHOLD)
				_fft.mul(a,b,c);
			else
				if ( d > KARA_DEG_THRESHOLD)
					_kara.mul(a,b,c);		
				else
					_classic.mul(a,b,c);
			/*
			mul.stop();multime+=mul;
			std::cout.width(30);
			std::cout<<"mul "<<b.size()<<"x"<<c.size()<<" : ";		
			std::cout.precision(3);
			std::cout<<mul.usertime()<<std::endl;
			*/
		}
		
		template< class Polynomial1, class Polynomial2, class Polynomial3>
		void midproduct (Polynomial1 &a, const Polynomial2 &b, const Polynomial3 &c) {
			linbox_check( 2*a.size() == c.size()+1);
			linbox_check( 2*b.size() == c.size()+1);			

			size_t d = b.size()+c.size();
			//std::cout<<"midp "<<a.size()<<" = "<<b.size()<<" x "<<c.size()<<"...\n"; 
			//Timer mul;
			//mul.start();
			if (d > FFT_DEG_THRESHOLD)
				_fft.midproduct(a,b,c);
			else
				if ( d > KARA_DEG_THRESHOLD)
					_kara.midproduct(a,b,c);		
				else
					_classic.midproduct(a,b,c);
			//std::cout<<"done\n";
			/*
			mul.stop();multime+=mul;
			std::cout.width(30);
			std::cout<<"mul "<<b.size()<<"x"<<c.size()<<" : ";		
			std::cout.precision(3);
			std::cout<<mul.usertime()<<std::endl;	
			*/
		}

		template< class Polynomial1, class Polynomial2, class Polynomial3>
		void midproductgen (Polynomial1 &a, const Polynomial2 &b, const Polynomial3 &c) {
			linbox_check( a.size()+b.size() == c.size()+1);
			
			size_t d = b.size()+c.size();
			//std::cout<<"midp "<<a.size()<<" = "<<b.size()<<" x "<<c.size()<<"...\n"; 
			//Timer mul;
			//mul.start();
			if (d > FFT_DEG_THRESHOLD)
				_fft.midproductgen(a,b,c);
			else
				if ( d > KARA_DEG_THRESHOLD)
					_kara.midproductgen(a,b,c);		
				else
					_classic.midproductgen(a,b,c);
			//std::cout<<"done\n";
			/*
			mul.stop();multime+=mul;
			std::cout.width(30);
			std::cout<<"mul "<<b.size()<<"x"<<c.size()<<" : ";		
			std::cout.precision(3);
			std::cout<<mul.usertime()<<std::endl;	
			*/
		}
	};

	
	template <class Field>
	class ClassicMulDomain {
	private:
		Field                       _F;
		BlasMatrixDomain<Field>   _BMD;
		MatrixDomain<Field>        _MD;
		
	public:
		
		ClassicMulDomain(const Field &F) : _F(F), _BMD(F), _MD(F) {}
		
		template< class Polynomial1, class Polynomial2, class Polynomial3>
		void mul(Polynomial1 &a, const Polynomial2 &b, const Polynomial3 &c) {
			
			linbox_check(a.size() >= (b.size()+c.size()-1));
			for (size_t i=0;i<b.size();++i){
				for (size_t j=0;j<c.size();++j)
					_BMD.axpyin(a[i+j],b[i],c[j]);
			}
		}

		template< class Polynomial1, class Polynomial2, class Polynomial3>
		void midproduct (Polynomial1 &a, const Polynomial2 &b, const Polynomial3 &c) {
			linbox_check( 2*a.size() == c.size()+1);
			linbox_check( 2*b.size() == c.size()+1);
			
			for (size_t i=0;i<b.size();++i){
				for (size_t j=0;j<c.size();++j){					
					if ((i+j<2*a.size()-1) && (i+j>=a.size()-1)){
						_BMD.axpyin(a[i+j - a.size()+1],b[i],c[j]);
					}
				}
			}			
		}

		template< class Polynomial1, class Polynomial2, class Polynomial3>
		void midproductgen (Polynomial1 &a, const Polynomial2 &b, const Polynomial3 &c) {
			linbox_check(a.size()+b.size() == c.size()+1);
		
			for (size_t i=0;i<b.size();++i){
				for (size_t j=0;j<c.size();++j){					
					if ((i+j>=a.size()-1) && (i+j<=c.size()-1)){
						_BMD.axpyin(a[i+j - a.size()+1],b[i],c[j]);
					}
				}
			}			
		}		
	};
	
	
	
	template<class _Field>
	class KaratsubaMulDomain {
	public:
		typedef _Field                         Field;
	private:
		Field                       _F;
		BlasMatrixDomain<Field>   _BMD;
		MatrixDomain<Field>        _MD;
		size_t                    _mul;
	public:
		
		KaratsubaMulDomain(const Field &F) : _F(F), _BMD(F), _MD(F) {_mul=0;}
		
		template< class Polynomial1, class Polynomial2, class Polynomial3>
		void mul(Polynomial1 &a, const Polynomial2 &b, const Polynomial3 &c) {			
			linbox_check(a.size() >= (b.size()+c.size()-1));
			Karatsuba_mul(a, 0, b, 0, b.size(), c, 0, c.size());
		}
		
		template< class Polynomial1, class Polynomial2, class Polynomial3>
		void midproduct(Polynomial1 &a, const Polynomial2 &b, const Polynomial3 &c) {
			linbox_check( 2*a.size() == c.size()+1);
			linbox_check( 2*b.size() == c.size()+1);
			midproduct_Karatsuba(a, b, c);	
		}
		
		template< class Polynomial1, class Polynomial2, class Polynomial3>
		void midproductgen(Polynomial1 &a, const Polynomial2 &b, const Polynomial3 &c) {
			linbox_check(a.size()+b.size() == c.size()+1);
			midproduct_Karatsubagen(a, b, c);	
		}

	protected:
		
		template< class Polynomial1, class Polynomial2, class Polynomial3>		
		void Karatsuba_mul(Polynomial1 &C, size_t shiftC,
				   const Polynomial2 &A, size_t shiftA, size_t degA,
				   const Polynomial3 &B, size_t shiftB, size_t degB){
			
			typedef  typename Polynomial1::value_type Coefficient;			       
			const Coefficient ZeroC(C[0].rowdim(), C[0].coldim());
			const Coefficient ZeroA(A[0].rowdim(), A[0].coldim());
			const Coefficient ZeroB(B[0].rowdim(), B[0].coldim());
						
			if ((degA == 1) || (degB == 1)) {
				
				if ((degA == 1) && (degB == 1))
					{_BMD.mul(C[shiftC],A[shiftA],B[shiftB]); _mul++;}
				else 
					if (degA == 1) 
						for (size_t i=0;i< degB;++i)
							{_BMD.mul(C[shiftC+i],A[shiftA],B[shiftB+i]);_mul++;}
					else 
						for (size_t i=0;i< degA;++i)
							{_BMD.mul(C[shiftC+i],A[shiftA+i],B[shiftB]);_mul++;}
			}
			else {
				size_t degA_low, degA_high, degB_low, degB_high, half_degA, half_degB, degSplit;
				half_degA= (degA & 1) + (degA >>1);
				half_degB= (degB & 1) + (degB >>1);
				degSplit= (half_degA > half_degB) ? half_degA : half_degB;
				
				degB_low = (degB < degSplit) ? degB : degSplit;
				degA_low = (degA < degSplit) ? degA : degSplit;
				degA_high= degA - degA_low;
				degB_high= degB - degB_low;
				
			
				// multiply low degrees
				Karatsuba_mul(C, shiftC, A, shiftA, degA_low, B, shiftB, degB_low);   
				
				// multiply high degrees (only if they are both different from zero)
				if ((degA_high !=0) && (degB_high != 0)) 	
					Karatsuba_mul(C, shiftC+(degSplit << 1), A, shiftA+degSplit, degA_high, B, shiftB+degSplit, degB_high);
								
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
				Karatsuba_mul(C_tmp, 0, A_tmp, 0, degA_low, B_tmp, 0, degB_low);
				
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


		template< class Polynomial1, class Polynomial2, class Polynomial3>
		void midproduct_Karatsuba(Polynomial1 &C, const Polynomial2 &A, const Polynomial3 &B){
		
			if (A.size() == 1){
				_BMD.mul(C[0],A[0],B[0]);
			}
			else {				
				size_t k0= A.size()>>1;
				size_t k1= A.size()-k0;
			
				typedef  typename Polynomial1::value_type Coefficient;			       
				const Coefficient ZeroC(C[0].rowdim(), C[0].coldim());
				const Coefficient ZeroA(A[0].rowdim(), A[0].coldim());
				const Coefficient ZeroB(B[0].rowdim(), B[0].coldim());
			
				std::vector<Coefficient> alpha(k1,ZeroC), beta(k1,ZeroC), gamma(k0,ZeroC);
				std::vector<Coefficient> A_low(k0, ZeroA), A_high(k1,ZeroA);
				std::vector<Coefficient> B1(2*k1-1,ZeroB), B2(2*k1-1,ZeroB);

				for (size_t i=0;i<k0;++i)
					A_low[i] = A[i];

				for (size_t i=k0;i<A.size();++i)
					A_high[i-k0] = A[i];

				for (size_t i=0;i<2*k1-1;++i){
					B1[i] = B[i];
					B2[i] = B[i+k1];
					_MD.addin(B1[i],B2[i]);
				}		
				midproduct_Karatsuba(alpha, A_high, B1);			

				if (k0 == k1) {
					for (size_t i=0;i<k1;++i)
						_MD.subin(A_high[i],A_low[i]);
					midproduct_Karatsuba(beta, A_high, B2);					
				}
				else {
					for (size_t i=1;i<k1;++i)
						_MD.subin(A_high[i],A_low[i-1]);
					midproduct_Karatsuba(beta, A_high, B2);
				}				
				
				std::vector<Coefficient> B3(2*k0-1,ZeroB);
				for (size_t i=0;i<2*k0-1;++i)
					_MD.add(B3[i],B[i+2*k1],B[i+k1]);
				
				midproduct_Karatsuba(gamma, A_low, B3);

				for (size_t i=0;i<k1;++i)
					_MD.sub(C[i],alpha[i],beta[i]);
				
				for (size_t i=0;i<k0;++i){
					C[k1+i]=gamma[i];
					_MD.addin(C[k1+i],beta[i]);
				}
			}
		}

		template< class Polynomial1, class Polynomial2, class Polynomial3>
		void midproduct_Karatsubagen(Polynomial1 &C, const Polynomial2 &A, const Polynomial3 &B){
		
			if (A.size() == 1){
				for (size_t i=0;i<B.size();i++)
					_BMD.mul(C[i],A[0],B[i]);
			}
			else {	
				size_t dA=A.size();
				size_t dB=B.size();
				size_t Ak0= dA>>1;
				size_t Ak1= dA-Ak0;
				size_t Bk0= (dB+1-dA)>>1;
				size_t Bk1= dB+1-dA-Bk0;
				std::cout<<C.size()<<" "<<A.size()<<" "<<B.size()<<"("<<Bk0<<","<<Bk1<<")\n";
				
				typedef  typename Polynomial1::value_type Coefficient;			       
				const Coefficient ZeroC(C[0].rowdim(), C[0].coldim());
				const Coefficient ZeroA(A[0].rowdim(), A[0].coldim());
				const Coefficient ZeroB(B[0].rowdim(), B[0].coldim());
			
				std::vector<Coefficient> A_low(Ak0,ZeroA), A_high(Ak1,ZeroA);
				std::vector<Coefficient> alpha(Bk1,ZeroC), beta(Bk1,ZeroC), gamma(Bk0,ZeroC);		      	
				std::vector<Coefficient> B1(Ak1+Bk1-1,ZeroB), B2(Ak1+Bk1-1,ZeroB), B3(Ak0+Bk0-1, ZeroB);;
				

				for (size_t i=0;i<Ak0;++i)
					A_low[i] = A[i];

				for (size_t i=Ak0;i<dA;++i)
					A_high[i-Ak0] = A[i];

				for (size_t i=0;i<Ak1+Bk1-1;++i){
					B1[i] = B[i];
					B2[i] = B[i+Ak1];
					_MD.addin(B1[i],B2[i]);
				}		
				midproduct_Karatsubagen(alpha, A_high, B1);			

				if (Ak0 == Ak1) {
					for (size_t i=0;i<Ak1;++i)
						_MD.subin(A_high[i],A_low[i]);
					midproduct_Karatsubagen(beta, A_high, B2);					
				}
				else {
					for (size_t i=1;i<Ak1;++i)
						_MD.subin(A_high[i],A_low[i-1]);
					midproduct_Karatsubagen(beta, A_high, B2);
				}				
				
				
				if (Bk0>0) {
					for (size_t i=0;i<Ak0+Bk0-1;++i)
						_MD.add(B3[i],B[i+Ak0+Bk1],B[i+Ak0]);				
					midproduct_Karatsubagen(gamma, A_low, B3);
				}
			

				for (size_t i=0;i<Bk1;++i)
					_MD.sub(C[i],alpha[i],beta[i]);
				
				for (size_t i=0;i<Bk0;++i){
					C[Bk1+i]=gamma[i];
					_MD.addin(C[Bk1+i],beta[i]);
				}
			}
		}

	}; // end of class KaratsubaMulDomain<Field, Matrix>

	template <class _Field>
	class SpecialFFTMulDomain;

	// FFT domain for every prime
	template <class _Field>
	class FFTMulDomain {
	public:
		typedef _Field                              Field;
		typedef typename Field::Element           Element;
		typedef SpecialFFTMulDomain<Field>  FFTDomainBase;

	private:
		Field                _F; 
		integer              _p;
		size_t         _fftsize;
			
	public:

		FFTMulDomain (const Field &F) :  _F(F){
			
			_F.characteristic(_p);

			_fftsize=0;
			//check if field is based on fft prime
			size_t p = _p;
			if (p&1){
				p-=1;
				do { p=p>>1; _fftsize++;} while(!(p&0x0001));
			}								
		}
		
		template< class Polynomial1, class Polynomial2, class Polynomial3>
		void mul(Polynomial1 &a, const Polynomial2 &b, const Polynomial3 &c) {
			linbox_check(a.size() >= (b.size()+c.size()-1));
			size_t deg     = b.size()+c.size()-1;
			size_t lpts = 0; 
			size_t pts =1; while (pts < deg) { pts= pts<<1; ++lpts; }

			typedef  typename Polynomial1::value_type Coefficient;			       
			const Coefficient ZeroC(c[0].rowdim(), c[0].coldim());
			const Coefficient ZeroA(a[0].rowdim(), a[0].coldim());
			const Coefficient ZeroB(b[0].rowdim(), b[0].coldim());
			
			// check if fft prime and good enough
			if (lpts < _fftsize){
				FFTDomainBase fftdomain(_F);
				fftdomain.mul(a, b, c);
			}
			else {				
				// computation done using CRT and few fft primes
				
				// get number of bits of feasible fft prime
				int k= b[0].coldim();
				size_t n=k;
				size_t ln=0;
				while ( k>0) {k>>=1; ln++;}
				
				// taking primes greater than current prime
				size_t bit = std::max((53-ln)>>1, _p.bitsize());

				// get number of necessary primes				
				integer ibound = n * _p * _p * std::max(b.size(), c.size());
				integer primesprod=1; size_t nbrprimes=1;
				RandomFFTPrime fftprime(bit, FFT_PRIME_SEED);
				std::vector<integer> lprimes(10); lprimes.resize(nbrprimes);
				lprimes[0] = fftprime.randomPrime();				
				primesprod = lprimes[0];
				while (primesprod < ibound) {
					++nbrprimes;
					lprimes.resize(nbrprimes);				
					do {lprimes[nbrprimes-1] = fftprime.randomPrime();} while (primesprod % lprimes[nbrprimes-1] == 0);					
					primesprod *= lprimes[nbrprimes-1]; 
				}
#ifdef FFT_TIMING
				std::cout<<"num of primes "<<nbrprimes<<"\n";
#endif
				// allocate fftprime fields
				Field * f_i = new Field[nbrprimes];

				// allocate polynomial matrices for multimodular results
				std::vector<Coefficient> * a_i = new std::vector<Coefficient>[nbrprimes];
				
				// set fftprimes, fftdomains, polynomial matrix results
				for (size_t i=0; i< nbrprimes; ++i){
					f_i[i] = Field(lprimes[i]);
					FFTDomainBase fftdomain(f_i[i]);					
					a_i[i] = std::vector<Coefficient>(a.size(), ZeroA);
					// does not work if original field representation is not Fp seen as integers mod p
					fftdomain.mul(a_i[i], b ,c);				
				}
				
				LinBox::Timer chrono;
				chrono.start();
				// reconstruct the solution modulo the original prime
				if (nbrprimes < 2) {
					for (size_t k=0;k<a.size();++k)
						for (size_t i=0;i<a[0].rowdim();++i)
							for (size_t j=0;j<a[0].coldim();++j){
								_F.init(a[k].refEntry(i,j), a_i[0][k].getEntry(i,j));
							}
				}
				else {
					integer * crt = new integer[nbrprimes];
					Element * crt_inv = new Element[nbrprimes];
					Element tmp;
					for (size_t i=0;i<nbrprimes; ++i){
						crt[i]=primesprod/lprimes[i];
						f_i[i].init(tmp,crt[i]);
						f_i[i].inv(crt_inv[i], tmp);
					}
					
					integer res,acc;
					for (size_t k=0;k<deg;++k)
						for (size_t i=0;i<a[0].rowdim();++i)
							for (size_t j=0;j<a[0].coldim();++j){
								acc= integer(0);
								for (size_t l=0;l<nbrprimes; ++l){  							
									f_i[l].mul(tmp, a_i[l][k].getEntry(i,j), crt_inv[l]);  
									res= f_i[l].convert(res,tmp);
									acc+= res*crt[l];
									if (acc > primesprod)
										acc-= primesprod;
								}
								_F.init(a[k].refEntry(i,j), acc);
							}
#ifdef FFT_TIMING
					chrono.stop();std::cout<<"reconstruction time: "<<chrono<<"\n";
#endif
					delete [] crt;
					delete [] crt_inv;
				}
				delete [] f_i;
				delete [] a_i;
			}			
		}
		
		template< class Polynomial1, class Polynomial2, class Polynomial3>
		void midproduct(Polynomial1 &a, const Polynomial2 &b, const Polynomial3 &c) {
			linbox_check(2*a.size() == c.size()+1 );
			linbox_check(2*b.size() == c.size()+1 );
			linbox_check(b[0].coldim() == c[0].rowdim());

			size_t m = b[0].rowdim();
			size_t k = b[0].coldim();
			size_t n = c[0].coldim();

			typedef  typename Polynomial1::value_type Coefficient;			       
			const Coefficient ZeroA(m,n), ZeroB(m,k), ZeroC(k,n);

			size_t deg  = c.size()+1;
			size_t lpts = 0; 
			size_t pts =1; while (pts < deg) { pts= pts<<1; ++lpts; }
			
			// check if fft prime and good enough
			if (lpts < _fftsize){
				FFTDomainBase fftdomain(_F);
				fftdomain.midproduct(a, b, c);
			}
			else {
				// computation done using CRT and few fft primes
				
				// get number of bits of feasible fft prime
				int k= b[0].coldim();
				size_t n=k;
				size_t ln=0;
				while ( k>0) {k>>=1; ln++;}
				
				// taking primes greater than current prime
				size_t bit = std::max((53-ln)>>1, _p.bitsize());

				// get number of necessary primes				
				integer ibound = n * _p * _p * std::max(b.size(), c.size());
				integer primesprod=1; size_t nbrprimes=1;
				RandomFFTPrime fftprime(bit, FFT_PRIME_SEED);
				std::vector<integer> lprimes(10); lprimes.resize(nbrprimes);
				lprimes[0] = fftprime.randomPrime();				
				primesprod = lprimes[0];
				while (primesprod < ibound) {
					++nbrprimes;
					lprimes.resize(nbrprimes);				
					do {lprimes[nbrprimes-1] = fftprime.randomPrime();} while (primesprod % lprimes[nbrprimes-1] == 0);					
					primesprod *= lprimes[nbrprimes-1]; 
				}
				
#ifdef FFT_TIMING			
				std::cout<<"num of primes "<<nbrprimes<<"\n";
#endif				

				// allocate fftprime fields
				Field * f_i = new Field[nbrprimes];

				// allocate polynomial matrices for multimodular results
				std::vector<Coefficient> * a_i = new std::vector<Coefficient>[nbrprimes];
				
				// set fftprimes, fftdomains, polynomial matrix results
				for (size_t i=0; i< nbrprimes; ++i){
					f_i[i] = Field(lprimes[i]);
					FFTDomainBase fftdomain(f_i[i]);					
					a_i[i] = std::vector<Coefficient>(a.size(), ZeroA);
					// does not work if original field representation is not Fp seen as integers mod p
					fftdomain.midproduct(a_i[i], b ,c);				
				}
				
				LinBox::Timer chrono;
				chrono.start();
				// reconstruct the solution modulo the original prime
				if (nbrprimes < 2) {
					for (size_t k=0;k<a.size();++k)
						for (size_t i=0;i<a[0].rowdim();++i)
							for (size_t j=0;j<a[0].coldim();++j){
								_F.init(a[k].refEntry(i,j), a_i[0][k].getEntry(i,j));
							}
				}
				else {
					integer * crt = new integer[nbrprimes];
					Element * crt_inv = new Element[nbrprimes];
					Element tmp;
					for (size_t i=0;i<nbrprimes; ++i){
						crt[i]=primesprod/lprimes[i];
						f_i[i].init(tmp,crt[i]);
						f_i[i].inv(crt_inv[i], tmp);
					}
					
					integer res,acc;
					for (size_t k=0;k<a.size();++k)
						for (size_t i=0;i<a[0].rowdim();++i)
							for (size_t j=0;j<a[0].coldim();++j){
								acc= integer(0);
								for (size_t l=0;l<nbrprimes; ++l){  							
									f_i[l].mul(tmp, a_i[l][k].getEntry(i,j), crt_inv[l]);  
									res= f_i[l].convert(res,tmp);
									acc+= res*crt[l];
									if (acc > primesprod)
										acc-= primesprod;
								}
								_F.init(a[k].refEntry(i,j), acc);
							}
					delete [] crt;
					delete [] crt_inv;
				}
#ifdef FFT_TIMING			
				chrono.stop();std::cout<<"reconstruction time: "<<chrono<<"\n";
#endif	
				delete [] f_i;
				delete [] a_i;
			}		
		}
	};

	
	// FFT Domain when prime is a FFT prime
	template <class _Field>
	class SpecialFFTMulDomain {
	public:
		typedef _Field                                    Field;
		typedef typename Field::Element                 Element;
		
	private:
		Field                      _F;
		integer                    _p;
		long                      _pl;
		MatrixDomain<Field>       _MD;
		BlasMatrixDomain<Field>  _BMD;
		double fftadd, fftmul, fftcopy;
		mutable             long _gen;
	public:

		SpecialFFTMulDomain(const Field &F) : _F(F), _MD(F), _BMD(F) {
			F.characteristic(_p);
			_pl = _p;
			fftadd= fftmul = fftcopy=0.;
			
			// find a pseudo primitive element of the multiplicative group _p -1						
			long m = _pl - 1;
			long k = 0;srand(time(NULL));
			while ((m & 1) == 0) {
				m = m >> 1;
				k++;
			}
			long long y,z,j;
			for (;;) {
				_gen = rand() % _pl; if (_gen <= 0) continue;
				
				z = 1; for (long i=0;i<m;++i) z = z*_gen % _pl;
				if (z == 1) continue;
				
				//_gen = z;
				j = 0;
				do {
					y = z;
					z = y*y % _pl;
					j++;
				} while (j != k && z != 1);				
				if (j == k) 
					break;
			}
		}
		
		template< class Polynomial1, class Polynomial2, class Polynomial3>
		void mul(Polynomial1 &a, const Polynomial2 &b, const Polynomial3 &c) {
#ifdef FFT_TIMING			
			Timer chrono;
			chrono.start();
#endif
			linbox_check(a.size() >=  b.size()+c.size()-1 );
			linbox_check(b[0].coldim() == c[0].rowdim());

			size_t m = b[0].rowdim();
			size_t k = b[0].coldim();
			size_t n = c[0].coldim();
			typedef  typename Polynomial1::value_type Coefficient;	
			const Coefficient ZeroA(m,n), ZeroB(m,k), ZeroC(k,n);

			size_t deg     = b.size()+c.size()-1;
			size_t lpts = 0; 
			size_t pts =1; while (pts < deg) { pts= pts<<1; ++lpts; }

#ifdef FFT_TIMING			
			std::cout<<"FFT: points "<<pts<<"\n";
#endif			
			if (_p%pts != 1) {
				std::cout<<"Error the prime is not a FFTPrime or it has too small power of 2\n";
				throw LinboxError("LinBox ERROR: bad FFT Prime\n");
			}
			
			long w ;
			// find a pseudo nth primitive root of unity
			for (;;) {						
				// compute the nth primitive root
				w=  (long) ::powmod(_gen, _pl>>lpts, _p);
				if ((w !=1) && (w != _pl-1))
					break;

				// find a pseudo primitive element of the multiplicative group _p-1
				long mm = _pl - 1;
				long kk = 0;srand(time(NULL));
				while ((mm & 1) == 0) {
					mm = mm >> 1;
					kk++;
				}
				long yy,zz,jj;
				for (;;) {
					_gen = rand() % _pl; if (_gen <= 0) continue;
					
					zz = 1; for (long i=0;i<mm;++i) zz = zz*_gen % _pl;
					if (zz == 1) continue;
					
					jj = 0;
					do {
						yy = zz;
						zz = yy*yy % _pl;
						jj++;
					} while (jj != kk && zz != 1);				
					if (jj == kk) 
						break;
				}
		 	}

			long inv_w;

			// compute w^(-1) mod p using extended euclidean algorithm 
			long x_int, y_int, q, tx, ty, temp;
			x_int = (long) _p;
			y_int = (long) w;
			tx = 0; ty = 1;
			while (y_int != 0) {
				q = x_int / y_int; // integer quotient
				temp = y_int; y_int = x_int - q * y_int;
				x_int = temp;
				temp = ty; ty = tx - q * ty;
				tx = temp;
			}			
			if (tx < 0) tx += (long) _p;			
			inv_w = tx;

			
			Element _w, _inv_w; 
			_F.init(_w, w);
			_F.init(_inv_w, inv_w); 
			std::vector<Element> pow_w(pts);
			std::vector<Element> pow_inv_w(pts);
			
			//std::cout<<"w: "<<w<<"\n w^-1: "<<inv_w<<"\n";
			//std::cout<<"degree: "<<pts<<"\n";

			// compute power of w and w^(-1)
			_F.init(pow_w[0],1);
			_F.init(pow_inv_w[0],1);		
			for (size_t i=1;i<pts;++i){
				_F.mul(pow_w[i], pow_w[i-1], _w);
				_F.mul(pow_inv_w[i], pow_inv_w[i-1], _inv_w);	
			}			 

			// compute reverse bit ordering
			size_t revbit[pts];
			for (long i = 0, j = 0; i < static_cast<long>(pts); i++, j = RevInc(j, lpts))
				revbit[i]=j;
			
			// set the data
			std::vector<Coefficient> fft_a(pts, ZeroA), fft_b(pts, ZeroB), fft_c(pts,ZeroC);
			for (size_t i=0;i<b.size();++i)
				fft_b[i]=b[i];
			for (size_t i=b.size();i<pts;++i)
				fft_b[i]=ZeroB;			
			for (size_t i=0;i<c.size();++i)
				fft_c[i]=c[i];
			for (size_t i=c.size();i<pts;++i)
				fft_c[i]=ZeroC;

			
			
#ifdef FFT_TIMING
			chrono.stop();
			std::cout<<"FFT: init                       : "<<chrono.usertime()<<"\n";
			chrono.clear();
			chrono.start();
#endif
			// compute the DFT of b and c using FFT (parallel if __LINBOX_HAVE_OPENMP)
			launch_FFT(fft_b, pts, pow_w);
			launch_FFT(fft_c, pts, pow_w);	
			
#ifdef FFT_TIMING
			chrono.stop();		       
			std::cout<<"FFT: DFT of inputs              : "<<chrono
				 <<" = "<<fftcopy<<" (copy), "<<fftadd<<" (add), "<<fftmul<<" (mul)\n";
			fftcopy=fftadd=fftmul=0.;
			chrono.clear();
			chrono.start();
#endif
			// do the multiplication componentwise
#ifdef __LINBOX_HAVE_OPENMP
#ifdef FFT_TIMING
			//LinBox::Timer chrono_mul[omp_get_max_threads()], chrono_mul_t[omp_get_max_threads()];
#endif
//std::cerr << "mul|| th: " << omp_get_max_threads() << ", pts: " << pts << ", n: " << fft_a[0].rowdim() << std::endl;
		
#pragma omp parallel for shared(fft_a,fft_b,fft_c) private(i) schedule(dynamic)
#endif			
			for (long i=0;i<static_cast<long>(pts);++i)
				// #ifdef FFT_TIMING
				// #ifdef __LINBOX_HAVE_OPENMP
				// 				{chrono_mul[omp_get_thread_num()].start();
				// #endif
				// #endif 
				_BMD.mul(fft_a[i], fft_b[i], fft_c[i]);							  
#ifdef FFT_TIMING
			chrono.stop();
#ifdef __LINBOX_HAVE_OPENMP
			//chrono_mul[omp_get_thread_num()].stop();chrono_mul_t[omp_get_thread_num()]+=chrono_mul[omp_get_thread_num()];}	
			//for (size_t i=0;i<omp_get_max_threads();i++)
			//std::cout<<"FFT: componentwise mul thread["<<i<<"] -> "<<chrono_mul_t[i]<<std::endl;
#endif						
			std::cout<<"FFT: componentwise mul total      : "<<chrono<<"\n";
			chrono.clear();
			chrono.start();	
#endif
			
			Element swapping;
			
			
			// reorder the term in the FFT according to reverse bit ordering
			// #ifdef __LINBOX_HAVE_OPENMP
			// #pragma omp parallel for shared(fft_a,revbit) private(Element) schedule(runtime)
			// #endif						
			for (size_t i=0; i< pts; ++i){
				if (revbit[i]>i){
					typename Coefficient::RawIterator it_a1=fft_a[i].rawBegin();
					typename Coefficient::RawIterator it_a2=fft_a[revbit[i]].rawBegin();
					for (; it_a1 != fft_a[i].rawEnd(); ++it_a1, ++it_a2){
						_F.assign(swapping,*it_a1);
						_F.assign(*it_a1, *it_a2);
						_F.assign(*it_a2,swapping);
					}					
				}
			}
			
			
#ifdef FFT_TIMING
			chrono.stop();
			std::cout<<"FFT: reverse bit ordering       : "<<chrono<<"\n";
			chrono.clear();
			chrono.start();	
#endif
	

			// compute the DFT inverse of fft_a
			launch_FFT(fft_a, pts, pow_inv_w);			
			//iterative_FFT(fft_a, pts, lpts, pow_inv_w);			
	
#ifdef FFT_TIMING	
			chrono.stop();
			std::cout<<"FFT: DFT inverse                : "<<chrono
				 <<" = "<<fftcopy<<" (copy), "<<fftadd<<" (add), "<<fftmul<<" (mul)\n";
			chrono.clear();
			chrono.start();	
#endif			

			// set the result according to bitreverse ordering and multiply by 1/pts
			Element inv_pts;
			_F.init(inv_pts, pts);
			_F.invin(inv_pts);
			
			// #ifdef __LINBOX_HAVE_OPENMP
			// #pragma omp parallel for shared(a,fft_a,revbit,inv_pts) schedule(stati
			// #endif
			
			for (long i=0; i< static_cast<long>(deg); ++i){
				//a[i] = fft_a[revbit[i]];
				//_MD.mulin(a[i], inv_pts);
				_MD.mul(a[i],fft_a[revbit[i]],inv_pts);
			}
#ifdef FFT_TIMING
			chrono.stop();
			std::cout<<"FFT: order and scale the result : "<<chrono.usertime()<<"\n\n";					
#endif
				
		}
			
		// middle product: a[0..n-1] = (b.c)[n..2n-1]
		template< class Polynomial1, class Polynomial2, class Polynomial3>
		void midproduct (Polynomial1 &a, const Polynomial2 &b, const Polynomial3 &c) {
			
			linbox_check(2*a.size() == c.size()+1 );
			linbox_check(2*b.size() == c.size()+1 );
			linbox_check(b[0].coldim() == c[0].rowdim());
			
			size_t m = b[0].rowdim();
			size_t k = b[0].coldim();
			size_t n = c[0].coldim();
			typedef  typename Polynomial1::value_type Coefficient;	
			const Coefficient ZeroA(m,n), ZeroB(m,k), ZeroC(k,n);

			size_t deg  = c.size()+1;
			size_t lpts = 0; 
			size_t pts =1; while (pts < deg) { pts= pts<<1; ++lpts; }
						
			if (_p%pts != 1) {
				std::cout<<"Error the prime is not a FFTPrime or it has too small power of 2\n";
				throw LinboxError("LinBox ERROR: bad FFT Prime\n");
			}
			
			long w ;
			// find a pseudo nth primitive root of unity
			for (;;) {
						
				// compute the nth primitive root
				w=  (long) ::powmod(_gen, _pl>>lpts, _p);
				//std::cout<<w<<" : "<<_gen<<"\n"<<(_pl>>lpts)<<"\n";

				if ((w !=1) && (w != _pl-1))
					break;

				// find a pseudo primitive element of the multiplicative group _p-1
				long mm = _pl - 1;
				long kk = 0;srand(time(NULL));
				while ((mm & 1) == 0) {
					mm = mm >> 1;
					kk++;
				}
				long yy,zz,jj;
				for (;;) {
					_gen = rand() % _pl; if (_gen <= 0) continue;
					
					zz = 1; for (long i=0;i<mm;++i) zz = zz*_gen % _pl;
					if (zz == 1) continue;
					
					jj = 0;
					do {
						yy = zz;
						zz = yy*yy % _pl;
						jj++;
					} while (jj != kk && zz != 1);				
					if (jj == kk) 
						break;
				}
			}

			long inv_w;

			// compute w^(-1) mod p using extended euclidean algorithm 
			long x_int, y_int, q, tx, ty, temp;
			x_int = (long) _p;
			y_int = (long) w;
			tx = 0; ty = 1;
			while (y_int != 0) {
				q = x_int / y_int; // integer quotient
				temp = y_int; y_int = x_int - q * y_int;
				x_int = temp;
				temp = ty; ty = tx - q * ty;
				tx = temp;
			}			
			if (tx < 0) tx += (long) _p;			
			inv_w = tx;
	
			
			Element _w, _inv_w;
			_F.init(_w,w);
			_F.init(_inv_w, inv_w); 
			std::vector<Element> pow_w(pts);
			std::vector<Element> pow_inv_w(pts);
			
			// compute power of w and w^(-1)
			_F.init(pow_w[0],1);
			_F.init(pow_inv_w[0],1);
			for (size_t i=1;i<pts;++i){
				_F.mul(pow_w[i], pow_w[i-1], _w);
				_F.mul(pow_inv_w[i], pow_inv_w[i-1], _inv_w);				
			}

			// compute reverse bit ordering
			size_t revbit[pts];
			for (long i = 0, j = 0; i < static_cast<long>(pts); i++, j = RevInc(j, lpts))
				revbit[i]=j;
		
			// set the data
			std::vector<Coefficient> fft_a(pts, ZeroA), fft_b(pts, ZeroB), fft_c(pts, ZeroC);
			for (size_t i=0;i<b.size();++i)
				fft_b[i]=b[b.size()-i-1];// reverse b
			for (size_t i=b.size();i<pts;++i)
				fft_b[i]=ZeroB;;			
			for (size_t i=0;i<c.size();++i)
				fft_c[i]=c[i];
			for (size_t i=c.size();i<pts;++i)
				fft_c[i]=ZeroC;


			// compute the DFT of b and DFT^-1 of c (parallel if __LINBOX_HAVE_OPENMP)
			launch_FFT(fft_b, pts, pow_w);							       	
			launch_FFT(fft_c, pts, pow_inv_w);

			// do the multiplication componentwise

#ifdef __LINBOX_HAVE_OPENMP
//std::cerr << "midproduct|| th: " << omp_get_max_threads() << ", pts: " << pts << ", n: " << fft_a[0].rowdim() << std::endl;
#pragma omp parallel for shared(fft_a,fft_b,fft_c) schedule(dynamic)
#endif
			for (long i=0;i<static_cast<long>(pts);++i)
				_BMD.mul(fft_a[i], fft_b[i], fft_c[i]);
						
			Element swapping;
			// reorder the term in the FFT according to reverse bit ordering
// #ifdef __LINBOX_HAVE_OPENMP
// #pragma omp parallel for shared(fft_a,revbit,inv_pts) schedule(runtime)
// #endif
			for (size_t i=0; i< pts; ++i){
				if (revbit[i]>i){
					typename Coefficient::RawIterator it_a1=fft_a[i].rawBegin();
					typename Coefficient::RawIterator it_a2=fft_a[revbit[i]].rawBegin();
					for (; it_a1 != fft_a[i].rawEnd(); ++it_a1, ++it_a2){
						_F.assign(swapping,*it_a1);
						_F.assign(*it_a1, *it_a2);
						_F.assign(*it_a2,swapping);
					}					
				}
			}

			// compute the DFT of fft_a (parallel if __LINBOX_HAVE_OPENMP)
			launch_FFT(fft_a, pts, pow_w);			
			
			
			// set the result according to bitreverse ordering and multiply by 1/pts
			Element inv_pts;
			_F.init(inv_pts, pts);
			_F.invin(inv_pts);
// #ifdef __LINBOX_HAVE_OPENMP
// #pragma omp parallel for shared(fft_a,revbit,inv_pts) schedule(static)
// #endif
				for (long i=0; i< static_cast<long>(a.size()); ++i){
				a[i] = fft_a[revbit[i]];
				_MD.mulin(a[i], inv_pts);
				}			
		}


// unbalanced middle product: a[0..n-1] = (b.c)[n..n+m-1]
		template< class Polynomial1, class Polynomial2, class Polynomial3>
		void midproductgen (Polynomial1 &a, const Polynomial2 &b, const Polynomial3 &c) {
			
			linbox_check(a.size()+b.size() == c.size()+1 );
			linbox_check(b[0].coldim() == c[0].rowdim());
			
			size_t m = b[0].rowdim();
			size_t k = b[0].coldim();
			size_t n = c[0].coldim();
			typedef  typename Polynomial1::value_type Coefficient;	
			const Coefficient ZeroA(m,n), ZeroB(m,k), ZeroC(k,n);

			size_t deg  = c.size()+1;
			size_t lpts = 0; 
			size_t pts =1; while (pts < deg) { pts= pts<<1; ++lpts; }
						
			if (_p%pts != 1) {
				std::cout<<"Error the prime is not a FFTPrime or it has too small power of 2\n";
				throw LinboxError("LinBox ERROR: bad FFT Prime\n");
			}
			
			long w ;
			// find a pseudo nth primitive root of unity
			for (;;) {
						
				// compute the nth primitive root
				w=  (long) ::powmod(_gen, _pl>>lpts, _p);
				//std::cout<<w<<" : "<<_gen<<"\n"<<(_pl>>lpts)<<"\n";

				if ((w !=1) && (w != _pl-1))
					break;

				// find a pseudo primitive element of the multiplicative group _p-1
				long mm = _pl - 1;
				long kk = 0;srand(time(NULL));
				while ((mm & 1) == 0) {
					mm = mm >> 1;
					kk++;
				}
				long yy,zz,jj;
				for (;;) {
					_gen = rand() % _pl; if (_gen <= 0) continue;
					
					zz = 1; for (long i=0;i<mm;++i) zz = zz*_gen % _pl;
					if (zz == 1) continue;
					
					jj = 0;
					do {
						yy = zz;
						zz = yy*yy % _pl;
						jj++;
					} while (jj != kk && zz != 1);				
					if (jj == kk) 
						break;
				}
			}

			long inv_w;

			// compute w^(-1) mod p using extended euclidean algorithm 
			long x_int, y_int, q, tx, ty, temp;
			x_int = (long) _p;
			y_int = (long) w;
			tx = 0; ty = 1;
			while (y_int != 0) {
				q = x_int / y_int; // integer quotient
				temp = y_int; y_int = x_int - q * y_int;
				x_int = temp;
				temp = ty; ty = tx - q * ty;
				tx = temp;
			}			
			if (tx < 0) tx += (long) _p;			
			inv_w = tx;
	
			
			Element _w, _inv_w;
			_F.init(_w,w);
			_F.init(_inv_w, inv_w); 
			std::vector<Element> pow_w(pts);
			std::vector<Element> pow_inv_w(pts);
			
			// compute power of w and w^(-1)
			_F.init(pow_w[0],1);
			_F.init(pow_inv_w[0],1);
			for (size_t i=1;i<pts;++i){
				_F.mul(pow_w[i], pow_w[i-1], _w);
				_F.mul(pow_inv_w[i], pow_inv_w[i-1], _inv_w);				
			}

			// compute reverse bit ordering
			size_t revbit[pts];
			for (long i = 0, j = 0; i < static_cast<long>(pts); i++, j = RevInc(j, lpts))
				revbit[i]=j;
		
			// set the data
			std::vector<Coefficient> fft_a(pts, ZeroA), fft_b(pts, ZeroB), fft_c(pts, ZeroC);
			for (size_t i=0;i<b.size();++i)
				fft_b[i]=b[b.size()-i-1];// reverse b
			for (size_t i=b.size();i<pts;++i)
				fft_b[i]=ZeroB;;			
			for (size_t i=0;i<c.size();++i)
				fft_c[i]=c[i];
			for (size_t i=c.size();i<pts;++i)
				fft_c[i]=ZeroC;


			// compute the DFT of b and DFT^-1 of c (parallel if __LINBOX_HAVE_OPENMP)
			launch_FFT(fft_b, pts, pow_w);							       	
			launch_FFT(fft_c, pts, pow_inv_w);

			// do the multiplication componentwise

#ifdef __LINBOX_HAVE_OPENMP
//std::cerr << "midproductgen|| th: " << omp_get_max_threads() << ", pts: " << pts << ", n: " << fft_a[0].rowdim() << std::endl;

#pragma omp parallel for shared(fft_a,fft_b,fft_c) schedule(dynamic)
#endif
			for (long i=0;i<pts;++i)
				_BMD.mul(fft_a[i], fft_b[i], fft_c[i]);
						
			Element swapping;
			// reorder the term in the FFT according to reverse bit ordering
// #ifdef __LINBOX_HAVE_OPENMP
// #pragma omp parallel for shared(fft_a,revbit,inv_pts) schedule(runtime)
// #endif
			for (size_t i=0; i< pts; ++i){
				if (revbit[i]>i){
					typename Coefficient::RawIterator it_a1=fft_a[i].rawBegin();
					typename Coefficient::RawIterator it_a2=fft_a[revbit[i]].rawBegin();
					for (; it_a1 != fft_a[i].rawEnd(); ++it_a1, ++it_a2){
						_F.assign(swapping,*it_a1);
						_F.assign(*it_a1, *it_a2);
						_F.assign(*it_a2,swapping);
					}					
				}
			}

			// compute the DFT of fft_a (parallel if __LINBOX_HAVE_OPENMP)
			launch_FFT(fft_a, pts, pow_w);			
			
			
			// set the result according to bitreverse ordering and multiply by 1/pts
			Element inv_pts;
			_F.init(inv_pts, pts);
			_F.invin(inv_pts);
// #ifdef __LINBOX_HAVE_OPENMP
// #pragma omp parallel for shared(fft_a,revbit,inv_pts) schedule(static)
// #endif
				for (long i=0; i< a.size(); ++i){
				a[i] = fft_a[revbit[i]];
				_MD.mulin(a[i], inv_pts);
				}			
		}

		//protected:

		inline long RevInc(long a, long k)
		{
			long j, m;
			
			j = k; 
			m = 1L << (k-1);
			
			while (j && (m & a)) {
				a ^= m;
				m >>= 1;
				j--;
			}
			if (j) a ^= m;
			return a;
		}

		template<class Coeff>
		inline void Butterfly (Coeff &A, Coeff &B, const Element &alpha) {
			typename Coeff::RawIterator it_a= A.rawBegin();
			typename Coeff::RawIterator it_b= B.rawBegin();
			Element tmp;
			for (; it_a != A.rawEnd(); ++it_a, ++it_b){
				_F.assign(tmp,*it_a);
				_F.addin(*it_a, *it_b);
				_F.sub(*it_b, tmp, *it_b);
				_F.mulin(*it_b, alpha);
			}
		}


		template<class Coefficient>
		void my(Coefficient &A, const Coefficient &B) {
			size_t n2 = A.rowdim()*A.coldim();
			Element       *aptr = A.getPointer();
			const Element *bptr = B.getPointer();
			for (size_t i=0;i<n2;++i){
				aptr[i]+= bptr[i];
			}
		}
		
		template<class Coefficient>
		void myAddSub(Coefficient &A, Coefficient &B) {
			size_t n2 = A.rowdim()*A.coldim();
			Element *aptr = A.getPointer();
			Element *bptr = B.getPointer();
			Element tmp;
			for (size_t i=0;i<n2;++i){
				tmp = aptr[i];
				_F.addin(aptr[i],bptr[i]);
				_F.sub(bptr[i], tmp, bptr[i]);
			}
		}

		
		template <class Polynomial>
		void launch_FFT (Polynomial &fft, size_t pts, const std::vector<Element> &pow_w){
#ifdef __LINBOX_HAVE_OPENMP
			// do blocking (by row) on the matrix coefficient to perform FFT in parallel on each block 
			size_t m,n;
			m=fft[0].rowdim();
			n=fft[0].coldim();			
			long nump=omp_get_max_threads();
			long nb_bsize=m%nump;
			long bsize   =m/nump+(nb_bsize?1:0);
			long lbsize  =m/nump;

// std::cerr << "launch_FFT|| th: " << omp_get_max_threads() << ", nump: " << nump << ", n: " << fft.size() << std::endl;
// #pragma omp parallel for shared(fft) schedule(dynamic)					
#pragma omp parallel for  schedule(dynamic)					
			for (int i=0;i<nump;i++){
#ifdef FFT_TIMING
				LinBox::Timer chrono1,chrono2;
				chrono1.start();		
#endif
				std::vector<DenseSubmatrix<Element> > block(fft.size());
				int row_idx,row_size;
				if (i>=nb_bsize) {
					row_size=lbsize;
					row_idx = nb_bsize*bsize+  (i-nb_bsize)*lbsize;
				} 
				else {
					row_size= bsize;
					row_idx = i*bsize;
				}
				
				for (size_t j=0;j<fft.size();j++)
					block[j]=DenseSubmatrix<Element>(fft[j],row_idx,0,row_size,n);

#ifdef FFT_TIMING
				chrono1.stop();
				chrono2.start();
#endif
				FFT(block,pts,pow_w);				

#ifdef FFT_TIMING
				chrono2.stop();
				std::cout<<"thread["<<omp_get_thread_num()<<"]: "<<chrono2<<"( "<<chrono1<<" )"<<std::endl;
#endif
			}
#else
			// call directly FFT code
			FFT(fft,pts,pow_w);
#endif	
		}

		template <class Polynomial>
		void FFT (Polynomial &fft, size_t n, const std::vector<Element> &pow_w, size_t idx_w=1, size_t shift=0){
			
			if (n != 1){

				size_t n2= n>>1;
				//size_t mn= fft[0].rowdim()* fft[0].coldim();
				//Coefficient tmp(fft[0].rowdim(), fft[0].coldim());
			
				//_MD.copy(tmp, fft[shift]);
				//_MD.addin(fft[shift],fft[shift+n2]);			
				//_MD.sub(fft[shift+n2], tmp, fft[shift+n2]);
				//myAddSub(fft[shift],fft[shift+n2]);
				Element one;_F.init(one,integer(1));
				Butterfly(fft[shift],fft[shift+n2],one);

				for (size_t i=1; i< n2; ++i){		
					Butterfly(fft[shift+i],fft[shift+i+n2],pow_w[idx_w*i]); 
					
					//_MD.copy(tmp, fft[shift+i]);
					//_MD.addin(fft[shift+i],fft[shift+i+n2]);					
					//_MD.sub(fft[shift+i+n2], tmp, fft[shift+i+n2]);
					//myAddSub(fft[shift+i],fft[shift+i+n2]);
					
					//_MD.mulin(fft[shift+i+n2],  pow_w[idx_w*i]);
					//FFLAS::fscal(_F, mn, pow_w[idx_w*i], fft[shift+i+n2].getPointer(), 1);
				}				
				FFT(fft, n2, pow_w, idx_w<<1, shift);				
				FFT(fft, n2, pow_w, idx_w<<1, shift+n2);						
			}
		}
			
		// fft entries are already in bit reverse order
		template< class Polynomial>
		void iterative_FFT (Polynomial &fft, size_t n, size_t ln, const std::vector<Element> &pow_w){

			typedef  typename Polynomial::value_type Coefficient;
			Coefficient tmp(fft[0].rowdim(), fft[0].coldim());
			
			if (ln == 0)
				return;
			if (ln == 1){
				_MD.copy(tmp, fft[0]);
				_MD.addin(fft[0], fft[1]);
				_MD.sub(fft[1], tmp, fft[1]);
				return;
			}

			// bottom level s = 1
			for (size_t i=0; i< n; i+=2) {
				_MD.copy(tmp, fft[i]);
				_MD.addin(fft[i], fft[i+1]);
				_MD.sub(fft[i+1], tmp, fft[i+1]);
			}
								
			// others levels s = 2..ln-1
			for (size_t s=2; s< ln; ++s){
				size_t m  = 1<<s;
				size_t m2 = 1<<(s-1);
				//size_t m4 = 1<<(s-2);		

				size_t w=(ln-s)<<1; 

				for (size_t i=0;i<n; i+=m){
					
					Coefficient *t, *t1,  *tt, *tt1, *u, *u1, *uu, *uu1;
					
					t  = &fft[i+m2];
					u  = &fft[i];
					t1 = &fft[i+1+m2]; _MD.mulin(*t1, pow_w[w]);
					u1 = &fft[i+1];

					for (size_t j=0; j<m2-2; j+=2){
						tt  = &fft[i+j+2+m2];_MD.mulin(*tt, pow_w[(j+2)*w]);
						uu  = &fft[i+j+2];
						tt1 = &fft[i+j+3+m2];_MD.mulin(*tt1, pow_w[(j+3)*w]);
						uu1 = &fft[i+j+3];
						
						_MD.copy(tmp, *u);
						_MD.addin(*u, *t);
						_MD.sub(*t,tmp,*t);
						
						_MD.copy(tmp,*u1);
						_MD.addin(*u1, *t1);
						_MD.sub(*t1,tmp,*t1);

						u  = uu;
						t  = tt;
						t1 = tt1;
						u1 = uu1;
					}

					_MD.copy(tmp, *u);
					_MD.addin(*u, *t);
					_MD.sub(*t,tmp,*t);
					
					_MD.copy(tmp,*u1);
					_MD.addin(*u1, *t1);
					_MD.sub(*t1,tmp,*t1);					
				}				
			}


			// top level s = ln
			//size_t m  = n;
			size_t m2 = 1<<(ln-1);
			//size_t m4 = 1<<(ln-2);
				
			
			_MD.copy(tmp, fft[0]);
			_MD.addin(fft[0], fft[m2]);
			_MD.sub(fft[m2],tmp,fft[m2]);
			
			_MD.copy(tmp,fft[1]);
			_MD.mulin(fft[m2+1], pow_w[1]);
			_MD.addin(fft[1], fft[m2+1]);
			_MD.sub(fft[m2+1],tmp,fft[m2+1]);

			for (size_t j=0; j< m2; j+=2){
				
				_MD.copy(tmp, fft[j]);
				_MD.mulin(fft[j+m2], pow_w[j>>1]);
				_MD.addin(fft[j], fft[j+m2]);
				_MD.sub(fft[j+m2],tmp,fft[j+m2]);
				
				_MD.copy(tmp,fft[j+1]);
				_MD.mulin(fft[j+m2+1], pow_w[j>>1]);
				_MD.addin(fft[j+1], fft[j+m2+1]);
				_MD.sub(fft[j+m2+1],tmp,fft[j+m2+1]);
			}
		}


	}; // end of class special FFT mul domain

} // end of namespace LinBox

#endif //__LINBOX_matpoly_mult_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

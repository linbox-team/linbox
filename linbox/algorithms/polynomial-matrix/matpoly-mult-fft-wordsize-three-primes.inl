/*
 * Copyright (C) 2015  Pascal Giorgi
 *
 * Written by Pascal Giorgi   <pascal.giorgi@lirmm.fr>
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
#ifndef __LINBOX_matpoly_mult_ftt_wordsize_three_primes_INL
#define __LINBOX_matpoly_mult_ftt_wordsize_three_primes_INL

#include "givaro/modular.h"
#include "fflas-ffpack/fflas-ffpack.h"
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/algorithms/polynomial-matrix/matpoly-mult-fft-wordsize-fast.inl"
#include "linbox/randiter/random-fftprime.h"
#include <cmath>

namespace LinBox {

	/***********************************************************************************
	 **** Polynomial Matrix Multiplication over Zp[x] with p (FFLAS prime) ***
	 *********************************x**************************************************/
	template<class Field>
	class PolynomialMatrixThreePrimesFFTMulDomain {
	public:
		// Polynomial matrix stored as a matrix of polynomial
		typedef PolynomialMatrix<Field,PMType::polfirst> MatrixP;
		// Polynomial matrix stored as a polynomial of matrix
		typedef PolynomialMatrix<Field,PMType::matfirst> PMatrix;
		//typedef Givaro::Modular<double>                ModField;
		typedef Field ModField;
	       		
	private:
		const Field              *_field;  // Read only
		uint64_t                      _p;
	  
	public:
		inline const Field & field() const { return *_field; }
	  
		PolynomialMatrixThreePrimesFFTMulDomain(const Field &F)
			: _field(&F), _p(field().cardinality())
		{
			if (integer(_p).bitsize()>29) {
				std::cout<<"MatPoly MUL FFT 3-primes: error initial prime has more than 29 bits exiting.."<<std::endl;
				throw LinboxError("LinBox ERROR: too large FFT Prime (more than 29 bits \n");
			}
		}

		template<typename Matrix1, typename Matrix2, typename Matrix3>
		void mul (Matrix1 &c, const Matrix2 &a, const Matrix3 &b, size_t max_rowdeg=0) const {
			linbox_check(a.coldim()==b.rowdim());
			// deg is the max rowdegree of the product
			size_t deg  = (max_rowdeg?max_rowdeg:a.size()+b.size()-2); //size_t deg  = a.size()+b.size()-1;
			c.resize(deg+1);
			size_t lpts = 0;
			size_t pts  = 1; while (pts <= deg) { pts= pts<<1; ++lpts; }
			// padd the input a and b to 2^lpts (convert to MatrixP representation)
			MatrixP a2(field(),a.rowdim(),a.coldim(),pts);
			MatrixP b2(field(),b.rowdim(),b.coldim(),pts);
			a2.copy(a,0,a.degree());
			b2.copy(b,0,b.degree());
			MatrixP c2(field(),c.rowdim(),c.coldim(),pts);
			integer bound=integer(_p-1)*integer(_p-1)
				*integer((uint64_t)a.coldim())*integer((uint64_t)std::min(a.size(),b.size()));
			mul_fft (lpts,c2, a2, b2, bound);
			c.copy(c2,0,deg);
		}

		void mul (MatrixP &c, const MatrixP &a, const MatrixP &b, size_t max_rowdeg=0) const {
			linbox_check(a.coldim()==b.rowdim());
			// deg is the max rowdegree of the product
			size_t deg  = (max_rowdeg?max_rowdeg:a.size()+b.size()-2); //size_t deg  = a.size()+b.size()-1;
			size_t lpts = 0;
			size_t pts  = 1; while (pts <= deg) { pts= pts<<1; ++lpts; }
			// padd the input a and b to 2^lpts
			MatrixP a2(field(),a.rowdim(),a.coldim(),pts);
			MatrixP b2(field(),b.rowdim(),b.coldim(),pts);
			a2.copy(a,0,a.degree());
			b2.copy(b,0,b.degree());
			// resize c to 2^lpts
			c.resize(pts);
			integer bound=integer(_p-1)*integer(_p-1)
				*integer((uint64_t)a.coldim())*integer((uint64_t)std::min(a.size(),b.size()));

			mul_fft (lpts,c, a2, b2, bound);
			c.resize(deg+1);
		}
		
		// a,b and c must have size: 2^lpts
		void mul_fft (size_t lpts, MatrixP &c, MatrixP &a, MatrixP &b, const integer& bound) const {
			size_t pts=c.size();			
			if ((_p-1) % pts == 0){
				PolynomialMatrixFFTPrimeMulDomain<ModField> fftprime_domain (field());
				fftprime_domain.mul_fft(lpts,c,a,b);
                		return;
			}			
			//std::cout<<"a:="<<a<<std::endl;
			//std::cout<<"b:="<<b<<std::endl;			
			linbox_check(a.coldim() == b.rowdim());
			size_t m = a.rowdim();
			size_t k = a.coldim();
			size_t n = b.coldim();
			uint64_t prime_max=maxFFTPrimeValue(k,pts); // CAREFUL: only for Modular<double>;
			std::vector<integer> bas;
			if (!RandomFFTPrime::generatePrimes (bas, prime_max, bound, lpts)){
				std::cout<<"COULD NOT FIND ENOUGH FFT PRIME in MatPoly 3-Primes FFT MUL exiting..."<<std::endl;
				std::cout<<"Parameters : ("<<m<<"x"<<k<<") ("<<k<<"x"<<n<<") with "<<pts<<" points with prime= "<<_p<<std::endl; 
				throw LinboxError("LinBox ERROR: not enough FFT Prime\n");
			}
			//std::cout<<"MUL FFT 3-PRIME found enough prime"<<std::endl;
			size_t num_primes = bas.size();
			std::vector<double> basis(num_primes);
			std::copy(bas.begin(),bas.end(),basis.begin());
	    
			std::vector<MatrixP*> c_i (num_primes);
			std::vector<ModField> f(num_primes,ModField(2));
			for (size_t l=0;l<num_primes;l++)
				f[l]=ModField(basis[l]);
	    
			for (size_t l=0;l<num_primes;l++){
				PolynomialMatrixFFTPrimeMulDomain<ModField> fftdomain (f[l]);
				MatrixP ai(f[l],m,k,pts);
				MatrixP bi(f[l],k,n,pts);
				if (basis[l]> _p) {
					//FFLAS::fassign(f[l],m*k*pts,a.getPointer(),1,ai.getPointer(),1);
					//FFLAS::fassign(f[l],k*n*pts,b.getPointer(),1,bi.getPointer(),1);
					// fassign is buggy (size < 2^31) with double
					std::copy(a.getPointer(),a.getPointer()+m*k*pts,ai.getPointer());
					std::copy(b.getPointer(),b.getPointer()+k*n*pts,bi.getPointer());
				}
				else {
					FFLAS::finit(f[l],m*k*pts,a.getPointer(),1,ai.getPointer(),1);
					FFLAS::finit(f[l],k*n*pts,b.getPointer(),1,bi.getPointer(),1);
				
				}
				c_i[l] = new MatrixP(f[l], m, n, pts);
 				fftdomain.mul_fft(lpts, *c_i[l], ai, bi);				
				//std::cout<<"pi:="<<(uint64_t)basis[l]<<std::endl;
				//std::cout<<"ci:="<<*c_i[l]<<std::endl;
			}

			// reconstruct the result with MRS
			typename Field::Element alpha,tmp;
			typename Field::Element beta=field().one;
			FFLAS::freduce(field(),m*n*pts,c_i[0]->getPointer(),1,c.getPointer(),1);
			//for(size_t k=0;k<m*n*pts;k++) field().reduce(c.getPointer()[k],c_i[0]->getPointer()[k]);

			for (size_t i=1;i<num_primes;i++){
				for(size_t j=0;j<i;j++){
					f[i].init(alpha,basis[j]);
					f[i].invin(alpha);
					FFLAS::fsubin (f[i],m*n*pts,c_i[j]->getPointer(),1,c_i[i]->getPointer(),1);
					//for(size_t k=0;k<m*n*pts;k++) {f[i].init(tmp,c_i[j]->getPointer()[k]); f[i].subin(c_i[i]->getPointer()[k], tmp); }
					FFLAS::fscalin(f[i],m*n*pts,alpha,c_i[i]->getPointer(),1);
					//for(size_t k=0;k<m*n*pts;k++) f[i].mulin(c_i[i]->getPointer()[k], alpha);					
				}
				field().init(tmp,basis[i-1]);
				field().mulin(beta,tmp);
				//field().mulin(beta,basis[i-1]);
				FFLAS::faxpy(field(),m*n*pts,beta,c_i[i]->getPointer(),1,c.getPointer(),1);
				//for(size_t k=0;k<m*n*pts;k++) field().axpyin(c.getPointer()[k], c_i[i]->getPointer()[k],beta);
			}
			
			//std::cout<<"c:="<<c<<std::endl;
			//#ifdef CHECK_MATPOL_MUL
			//std::cerr<<"(3 - prime CRT ) - "<<_p<<" - ";
			//check_mul(c,a,b,c.size());
			//#endif

			for (size_t i=0;i<num_primes;i++)
				delete c_i[i];
		
			
		}

		// compute  c= (a*b x^(-n0-1)) mod x^n1
		// by defaut: n0=c.size() and n1=2*c.size()-1;
		template<typename Matrix1, typename Matrix2, typename Matrix3>
		void midproduct (Matrix1 &c, const Matrix2 &a, const Matrix3 &b,
				 bool smallLeft=true, size_t n0=0,size_t n1=0) const {
			linbox_check(a.coldim()==b.rowdim());
			size_t hdeg = (n0==0?c.size():n0);
			size_t deg  = (n1==0?2*hdeg-1:n1);
			linbox_check(c.size()>=deg-hdeg);
			if (smallLeft){
				linbox_check(b.size()<hdeg+deg);
			}
			else
				linbox_check(a.size()<hdeg+deg);

			size_t lpts = 0;
			size_t pts  = 1; while (pts < deg) { pts= pts<<1; ++lpts; }
			// padd the input a and b to 2^lpts (use MatrixP representation)
			MatrixP a2(field(),a.rowdim(),a.coldim(),pts);
			MatrixP b2(field(),b.rowdim(),b.coldim(),pts);
			MatrixP c2(field(),c.rowdim(),c.coldim(),pts);
			a2.copy(a,0,a.size()-1);
			b2.copy(b,0,b.size()-1);

			// reverse the element of the smallest polynomial according to h(x^-1)*x^(hdeg)
			if (smallLeft)
				for (size_t j=0;j<a2.rowdim()*a2.coldim();j++)
					for (size_t i=0;i<hdeg/2;i++)
						std::swap(a2.ref(j,i),a2.ref(j,hdeg-1-i));
			else
				for (size_t j=0;j<b2.rowdim()*b2.coldim();j++)
					for (size_t i=0;i<hdeg/2;i++)
						std::swap(b2.ref(j,i),b2.ref(j,hdeg-1-i));
			integer bound=integer(_p-1)*integer(_p-1)
				*integer((uint64_t)a.coldim())*integer((uint64_t)std::min(a.size(),b.size()));
			
			midproduct_fft (lpts,c2, a2, b2, bound, smallLeft);
			c.copy(c2,0,c.size()-1);
		}

		
		// a,b and c must have size: 2^lpts
		// -> a must have been already reversed according to the midproduct algorithm
		void midproduct_fft (size_t lpts, MatrixP &c, MatrixP &a, MatrixP &b,
				     const integer& bound, bool smallLeft=true) const {
			size_t pts=c.size();			
			if ((_p-1) % pts == 0){
				//std::cerr<<"3-prime FFT midp switching to FFTPrime  "<<std::endl;
				PolynomialMatrixFFTPrimeMulDomain<ModField> fftprime_domain (field());
				fftprime_domain.midproduct_fft(lpts,c,a,b,smallLeft);
				return;
			}
			size_t m = a.rowdim();
			size_t k = a.coldim();
			size_t n = b.coldim();

			// compute bit size of feasible prime for FFLAS
			// size_t _k=k,lk=0;
			// while ( _k ) {_k>>=1; ++lk;}
			// size_t prime_bitsize= (53-lk)>>1;

			// compute max prime value for FFLAS
			//uint64_t prime_max= std::min(uint64_t(std::sqrt( (1ULL<<53) / k)+1), uint64_t(Givaro::Modular<double>::maxCardinality())) 
			uint64_t prime_max=maxFFTPrimeValue(k,pts); // CAREFUL: only for Modular<double>;
			
			std::vector<integer> bas;
			if (!RandomFFTPrime::generatePrimes (bas, prime_max, bound, lpts)){
				std::cout<<"COULD NOT FIND ENOUGH FFT PRIME  in MatPoly 3-Primes FFT MIDP exiting..."<<std::endl;
				std::cout<<"Parameters : ("<<m<<"x"<<k<<") ("<<k<<"x"<<n<<") with "<<pts<<" points with prime= "<<_p<<std::endl; 
				throw LinboxError("LinBox ERROR: not enough FFT Prime\n");
			}
			size_t num_primes = bas.size();
 
			std::vector<double> basis(num_primes);
			std::copy(bas.begin(),bas.end(),basis.begin());
	    
			std::vector<MatrixP*> c_i (num_primes);
			std::vector<ModField> f(num_primes,ModField(2));
			for (size_t l=0;l<num_primes;l++)
				f[l]=ModField(basis[l]);
	    
			for (size_t l=0;l<num_primes;l++){
				//std::cerr<<"3-prime FFT midp over "; f[l].write(std::cerr)<<std::endl;
				PolynomialMatrixFFTPrimeMulDomain<ModField> fftdomain (f[l]);
				MatrixP ai(f[l],m,k,pts);
				MatrixP bi(f[l],k,n,pts);
				if (basis[l]> _p) {
					//FFLAS::fassign(f[l],m*k*pts,a.getPointer(),1,ai.getPointer(),1);
					//FFLAS::fassign(f[l],k*n*pts,b.getPointer(),1,bi.getPointer(),1);
					// fassign is buggy (size < 2^31) with double
					std::copy(a.getPointer(),a.getPointer()+m*k*pts,ai.getPointer());
					std::copy(b.getPointer(),b.getPointer()+k*n*pts,bi.getPointer());				
				}
				else {
					FFLAS::finit(f[l],m*k*pts,a.getPointer(),1,ai.getPointer(),1);
					FFLAS::finit(f[l],k*n*pts,b.getPointer(),1,bi.getPointer(),1);
				
				}			       
				c_i[l] = new MatrixP(f[l], m, n, pts);
				fftdomain.midproduct_fft(lpts, *c_i[l], ai, bi,smallLeft);				
				//std::cout<<"pi:="<<(uint64_t)basis[l]<<std::endl;
				//std::cout<<"ci:="<<*c_i[l]<<std::endl;
			}
	    
			// reconstruct the result with MRS
			typename Field::Element alpha,tmp;
			typename Field::Element beta=field().one;
			FFLAS::freduce(field(),m*n*pts,c_i[0]->getPointer(),1,c.getPointer(),1);
			for (size_t i=1;i<num_primes;i++){
				for(size_t j=0;j<i;j++){
					f[i].init(alpha,basis[j]);
					f[i].invin(alpha);
					FFLAS::fsubin (f[i],m*n*pts,c_i[j]->getPointer(),1,c_i[i]->getPointer(),1);
					FFLAS::fscalin(f[i],m*n*pts,alpha,c_i[i]->getPointer(),1);
				}
 				field().init(tmp,basis[i-1]);
				field().mulin(beta,tmp);
				//field().mulin(beta,basis[i-1]);
				FFLAS::faxpy(field(),m*n*pts,beta,c_i[i]->getPointer(),1,c.getPointer(),1);
			}

			//std::cout<<"c:="<<c<<std::endl;
			
			for (size_t i=0;i<num_primes;i++)
				delete c_i[i];
		
		}
		
	};
} // end of namespace LinBox

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

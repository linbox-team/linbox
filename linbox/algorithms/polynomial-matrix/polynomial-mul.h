/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/*
 * Copyright (C) 2016  Pascal Giorgi
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
#include <givaro/zring.h>
#include "linbox/ring/modular.h"
#include "linbox/randiter/random-fftprime.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/util/error.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-fft-transform.h"
#include <fflas-ffpack/field/rns-double.h>

namespace LinBox {
	
	template<typename T>
	class PolynomialFFTMulDomain;


	
	template<class Field>
	class PolynomialFFTPrimeMulDomain {
	private:
		const Field              *_field;  // Read only
		uint64_t                      _p;
		FieldSimd<Field>      _fieldSIMD;
	
	public:
		typedef typename Field::Element       Element;
		typedef std::vector<Element, AlignedAllocator<Element, Alignment::CACHE_LINE>> Polynomial; 
		
		inline const Field & field() const { return *_field; }
	
		PolynomialFFTPrimeMulDomain(const Field &F)
			: _field(&F), _p(field().cardinality()), _fieldSIMD(F) {}

		template<typename Poly>
		void mul (Poly& c,  const Poly& a,  const Poly& b){
			FFT_PROFILE_START(1);
			size_t deg  = a.size()+b.size()-1;
			size_t lpts = 0;
			size_t pts  = 1; while (pts < deg) { pts= pts<<1; ++lpts; }
			// padd the input a and b to 2^lpts
			Polynomial a2(pts), b2(pts), c2(pts);
			std::copy(a.begin(),a.end(),a2.begin());
			std::copy(b.begin(),b.end(),b2.begin());
			FFT_PROFILING(1,"padding/copying input");
			mul_fft (lpts,c2.data(), a2.data(), b2.data());
			std::copy(c2.begin(),c2.begin()+deg,c.begin());
			FFT_PROFILING(1,"copying output");		
		}
		
		// a,b,c are of size 2^lpts and their entries will be modified
		void mul_fft (size_t lpts, Element* c, Element* a, Element* b){
			FFT_PROFILE_START(1);
			size_t pts=1<<lpts;		
			if ((_p-1) % pts != 0) {
				std::cout<<"Error the prime is not a FFTPrime or it has too small power of 2\n";
				std::cout<<"prime="<<_p<<std::endl;
				std::cout<<"nbr points="<<pts<<std::endl;
				throw LinboxError("LinBox ERROR: bad FFT Prime\n");
			}
			FFT_transform<Field> FFTer (field(), lpts);
			FFT_transform<Field> FFTinv (field(), lpts, FFTer.getInvRoot());
			FFT_PROFILING(1,"FFT init");

			FFTer.FFT_DIF(a);
			FFTer.FFT_DIF(b);
			FFT_PROFILING(1,"direct FFT_DIF");

#ifdef __FFLASFFPACK_USE_SIMD
			//cout<<"OK SIMD USED\n";
			using simd = Simd<Element>;
			using vect_t = typename simd::vect_t;
			vect_t av,bv,cv;
			size_t j = 0;
			for( ; j < ROUND_DOWN(pts, simd::vect_size) ; j+=simd::vect_size){
				av = simd::load(a+j);
				bv = simd::load(b+j);
				cv = simd::mul(av,bv);
				_fieldSIMD.mod(cv);
				simd::store(c+j, cv);
			}
			for( ; j < pts ; ++j){
				field().mul(c[j],a[j],b[j]);
			}						
#else			
			for(size_t i=0;i<pts;i++)
				field().mul(c[i],a[i],b[i]);
#endif

			FFT_PROFILING(1,"Pointwise mult");
#if DEBUG_FFT
			cout<<"w="<<FFTer.getRoot()<<endl;;
			cout<<"A(w)= "; std::copy(a, a+pts, std::ostream_iterator<Element>(std::cout, " "));cout<<endl;
			cout<<"B(w)= "; std::copy(b, b+pts, std::ostream_iterator<Element>(std::cout, " "));cout<<endl;
			cout<<"C(w)= "; std::copy(c, c+pts, std::ostream_iterator<Element>(std::cout, " "));cout<<endl;
#endif 
			FFTinv.FFT_DIT(c);
			FFT_PROFILING(1,"inverse FFT_DIT");
			typename Field::Element inv_pts;
			field().init(inv_pts, pts);
			field().invin(inv_pts);

			FFLAS::fscalin(field(),pts, inv_pts, c,1);
			FFT_PROFILING(1,"scaling the result");
			//for(size_t i=0;i<pts;i++)
			//	field().mulin(c[i],inv_pts);
		}
	};


	template<class Field>
	class PolynomialThreePrimesFFTMulDomain {
	public:
		typedef Field ModField;
		typedef typename Field::Element       Element;
		typedef std::vector<Element, AlignedAllocator<Element, Alignment::CACHE_LINE>> Polynomial; 

	private:
		const Field              *_field;  // Read only
		uint64_t                      _p;
	  
	public:
		inline const Field & field() const { return *_field; }
	  
		PolynomialThreePrimesFFTMulDomain(const Field &F)
			: _field(&F), _p(field().cardinality())
		{
			if (integer(_p).bitsize()>29) {
				std::cout<<"MatPoly MUL FFT 3-primes: error initial prime has more than 29 bits exiting.."<<std::endl;
				throw LinboxError("LinBox ERROR: too large FFT Prime (more than 29 bits \n");
			}
		}
		template<typename Poly>
		void mul (Poly& c, const Poly& a, const Poly& b){
			FFT_PROFILE_START(1);
			size_t deg  = a.size()+b.size()-1;
			size_t lpts = 0;
			size_t pts  = 1; while (pts < deg) { pts= pts<<1; ++lpts; }
			// padd the input a and b to 2^lpts
			Polynomial a2(pts), b2(pts), c2(pts);
			std::copy(a.begin(),a.end(),a2.begin());
			std::copy(b.begin(),b.end(),b2.begin());			
			FFT_PROFILING(1,"padding/copying input");
			mul_fft (lpts,c2.data(), a2.data(), b2.data());
			std::copy(c2.begin(),c2.begin()+deg,c.begin());
			FFT_PROFILING(1,"copying output");
		}

		void mul_fft (size_t lpts, Element* c, Element* a,  Element* b) {
			FFT_PROFILE_START(1);
			size_t pts=(uint64_t(1)<<lpts);
			if ((_p-1) % pts == 0){
				PolynomialFFTPrimeMulDomain<ModField> fftprime_domain (field());
				fftprime_domain.mul_fft(lpts,c,a,b);
				return;
			}			
			
			uint64_t prime_max= uint64_t(Givaro::Modular<double>::maxCardinality());
			RandomFFTPrime RdFFT(prime_max);
			integer bound=_p*_p*integer((uint64_t)pts);
			std::vector<integer> bas;
			if (!RdFFT.generatePrimes(lpts,bound,bas)){
				std::cout<<"COULD NOT FIND ENOUGH FFT PRIME in Polynomial 3-Primes FFT MUL exiting..."<<std::endl;
				std::cout<<"Parameters : ("<<pts<<" points with prime= "<<_p<<std::endl; 
				throw LinboxError("LinBox ERROR: not enough FFT Prime\n");
			}
			size_t num_primes = bas.size();
			std::vector<double> basis(num_primes);
			std::copy(bas.begin(),bas.end(),basis.begin());

#ifdef FFT_PROFILER
			std::cout << "number of FFT primes :" << num_primes << std::endl;
			std::cout << "max prime            : "<<prime_max<<" ("<<integer(prime_max).bitsize()<<")"<<std::endl;
			std::cout << "bitsize of the output: "<<bound.bitsize();
			std::cout <<" +++++++++++++++++++++++++++++++"<<std::endl;	
#endif
			Polynomial c_i(num_primes*pts);
			Polynomial ai(pts),bi(pts);
			FFT_PROFILING(1,"3Primes: init");
			for (size_t l=0;l<num_primes;l++){
				ModField F(basis[l]);
				PolynomialFFTPrimeMulDomain<ModField> fftdomain (F);
				//FFLAS::fassign(F,pts,a,1,ai.data(),1);
				//FFLAS::fassign(F,pts,b,1,bi.data(),1);
				//fftdomain.mul_fft(lpts, c_i.data()+l*pts, ai.data(), bi.data());
				fftdomain.mul_fft(lpts, c_i.data()+l*pts, a, b);				
			}
			FFT_PROFILING(1,"3Primes: red+mul");
	    
			// reconstruct the result with MRS
			typename Field::Element alpha;
			typename Field::Element beta=field().one;
			FFLAS::freduce(field(),pts,c_i.data(),1,c,1);
			for (size_t i=1;i<num_primes;i++){
				ModField F(basis[i]);
				for(size_t j=0;j<i;j++){
					F.init(alpha,basis[j]);
					F.invin(alpha);
					FFLAS::fsubin (F,pts,c_i.data()+j*pts,1,c_i.data()+i*pts,1);
					FFLAS::fscalin(F,pts,alpha,c_i.data()+i*pts,1);
				}
				field().mulin(beta,basis[i-1]);
				FFLAS::faxpy(field(),pts,beta,c_i.data()+i*pts,1,c,1);
			}
			FFT_PROFILING(1,"3Primes: reconstruct");

		}

	};
	
        template <class T1, class T2>
        class PolynomialFFTMulDomain<Givaro::Modular<T1,T2> > {
	public:
		typedef Givaro::Modular<T1,T2> Field;
        private:
		
                const Field            *_field;  // Read only
                uint64_t                    _p;
        public:
                inline const Field & field() const { return *_field; }

                PolynomialFFTMulDomain (const Field& F) : _field(&F), _p(F.cardinality()) {}

                template<typename Poly>
                void mul (Poly &c, const Poly &a,  const Poly &b) {
                        uint64_t pts= 1<<(integer((uint64_t)a.size()+b.size()-1).bitsize());
                        if ( _p< (1<<26)  &&  ((_p-1) % pts)==0){
				PolynomialFFTPrimeMulDomain<Field> MulDom(field());
                                MulDom.mul(c,a,b);
                        }
                        else {
				if (_p< 536870912ULL){
					PolynomialThreePrimesFFTMulDomain<Field> MulDom(field());
					MulDom.mul(c,a,b);
				}
				else {
					std::cout<<"NOT YET IMPLEMENTED"<<std::endl;
				}
			}
		}
	};
		
		
	
		
	
	template<>
	class PolynomialFFTMulDomain <Givaro::ZRing<integer> > {
	public:
		typedef Givaro::ZRing<integer>       Field;
		typedef Givaro::Modular<double>      ModField;
		typedef std::vector<integer>       Polynomial;
		typedef std::vector<double, AlignedAllocator<double, Alignment::CACHE_LINE>> PolynomialDbl; 

	private:
		const Field     *_field;
		integer           _maxnorm;
	

		size_t logmax(const Polynomial& A) const {
			size_t mm=A[0].bitsize();
			for(size_t k=1;k<A.size();k++)
				mm=std::max(mm,A[k].bitsize());
			return mm;
		}

	public:
		void getFFTPrime(uint64_t prime_max, size_t lpts, integer bound, std::vector<integer> &bas){

			RandomFFTPrime RdFFT(prime_max);
			size_t nbp=0;
			if (!RdFFT.generatePrimes(lpts,bound,bas)){
				integer MM=1;
				for(std::vector<integer>::size_type i=0;i<bas.size();i++)
					MM*=bas[i];
				RandomPrimeIter Rd(integer(prime_max).bitsize());
				integer tmp;
				do {
					do {Rd.random(tmp);}
					while (MM%tmp==0 || tmp>prime_max);
					bas.push_back(tmp);
					nbp++;
					MM*=tmp;
				} while (MM<bound);	
			}
#ifdef VERBOSE_FFT
			std::cout<<"MatPoly Multiprecision FFT : using "<<bas.size()-nbp<<" FFT primes and "<<nbp<<" normal primes "<<std::endl;
#endif
			for(auto i: bas)
				if (i>prime_max) std::cout<<"ERROR\n";
		}

		PolynomialFFTMulDomain (const Field& F) : _field(&F) {}

		template<typename Poly>
		void mul (Poly &c, const Poly &a, const Poly &b){
			FFT_PROFILE_START(2);
			integer maxA,maxB;
				maxA=maxB=_maxnorm;			
				if (_maxnorm==0){
					maxA=1;maxA<<=uint64_t(logmax(a));
					maxB=1;maxB<<=uint64_t(logmax(b));
				}
				integer bound=2*maxA*maxB*uint64_t(std::min(a.size(),b.size()));
				FFT_PROFILING(2,"max norm computation");
				
				mul_crtla(c,a,b,maxA,maxB,bound);				
			}

		template<typename Poly>
		void mul_crtla(Poly &c, const Poly &a, const Poly &b,
			       const integer& maxA, const integer& maxB, const integer& bound) {
			FFT_PROFILE_START(2);
			size_t s= a.size()+b.size()-1;
			c.resize(s);
			size_t lpts=0;
			size_t pts  = 1; while (pts < s) { pts= pts<<1; ++lpts; }

			// compute max prime value for FFLAS      
			uint64_t prime_max= uint64_t(Givaro::Modular<double>::maxCardinality());
			std::vector<integer> bas;
			getFFTPrime(prime_max,lpts,bound,bas);
			std::vector<double> basis(bas.size());
			std::copy(bas.begin(),bas.end(),basis.begin());
			FFPACK::rns_double RNS(basis);
			size_t num_primes = RNS._size;
#ifdef FFT_PROFILER
			std::cout << "number of FFT primes :" << num_primes << std::endl;
			std::cout << "max prime            : "<<prime_max<<" ("<<integer(prime_max).bitsize()<<")"<<std::endl;
			std::cout << "bitsize of the output: "<<bound.bitsize()
				  <<"( "<< RNS._M.bitsize()<<" )"<<std::endl;
			std::cout <<" +++++++++++++++++++++++++++++++"<<std::endl;	
#endif
			FFT_PROFILING(2,"init of CRT approach");
			// reduce t_a and t_b modulo each FFT primes
			//double* a_mod= new double[pts*num_primes];
			//double* b_mod= new double[pts*num_primes];
			//double* c_mod= new double[pts*num_primes];
			PolynomialDbl amod(pts*num_primes,0);
			PolynomialDbl bmod(pts*num_primes,0);
			PolynomialDbl cmod(pts*num_primes,0);
			double* a_mod=amod.data();
			double* b_mod=bmod.data();
			double* c_mod=cmod.data();
			
			RNS.init(1, a.size(), a_mod, pts, a.data(), a.size(), maxA);
			RNS.init(1, b.size(), b_mod, pts, b.data(), b.size(), maxB);
			FFT_PROFILING(2,"reduction mod pi of input polynomial");

			// FFT_PROFILE_START(2);
			// auto sp=SPLITTER();
			// PARFOR1D(l,num_primes,sp,
			for (size_t l=0;l<num_primes;l++)
				{
					//FFT_PROFILE_START;
					ModField f(RNS._basis[l]);		
					PolynomialThreePrimesFFTMulDomain<ModField> fftdomain (f);       
					fftdomain.mul_fft(lpts, c_mod+l*pts,a_mod+l*pts,b_mod+l*pts);	
					//FFT_PROFILE_GET(tMul);
				}
#ifdef DEBUG_POLMUL	
			std::cout<<"*** POL MUL RNS ***"<<std::endl;
			for(size_t i=0;i<num_primes;i++)
				{std::copy(c_mod+i*pts, c_mod+(i+1)*pts, std::ostream_iterator<integer>(std::cout, ", "));std::cout<<std::endl;}
			std::cout<<"*******************"<<std::endl;
#endif
			//)
			FFT_PROFILING(2,"FFTprime mult");
			//delete[] a_mod;
			//delete[] b_mod;
	
			// reconstruct the result in C
			RNS.convert(1,s,0,c.data(), s, c_mod, pts);
			//delete[] c_mod;
			FFT_PROFILING(2,"k prime reconstruction");
		}



	};
	
} // end of namespace


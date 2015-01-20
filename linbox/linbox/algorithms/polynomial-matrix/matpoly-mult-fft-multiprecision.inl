// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/*
 * Copyright (C) 2013  Pascal Giorgi
 *                     Romain Lebreton
 *
 * Written by Pascal Giorgi   <pascal.giorgi@lirmm.fr>
 *            Romain Lebreton <lebreton@lirmm.fr>
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
#ifndef __LINBOX_matpoly_mult_ftt_multiprecision_INL
#define __LINBOX_matpoly_mult_ftt_multiprecision_INL

#include "linbox/field/unparametric.h"
#include "linbox/field/modular.h"
#include "linbox/randiter/random-fftprime.h"
#include <fflas-ffpack/field/rns-double.h>

namespace LinBox{

	/***************************************************
	 **** Polynomial Matrix Multiplication over Z[x] ***
	 ***************************************************/
	template<>
	class PolynomialMatrixFFTMulDomain<Givaro::UnparametricRing<integer> > {
	public:
		typedef Givaro::UnparametricRing<integer>       IntField;
		typedef Givaro::Modular<int32_t> ModField;
		typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,ModField> MatrixP_F; // Polynomial matrix stored as a polynomial of matrix
		typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,IntField> MatrixP_I; // Polynomial matrix stored as a polynomial of matrix

	private:
		const IntField     *_field;
		integer           _maxnorm;

		template<typename PMatrix1>
		size_t logmax(const PMatrix1 A) const {
			size_t mm=A.get(0,0,0).bitsize();
			for(size_t k=0;k<A.size();k++)
				for (size_t i=0;i<A.rowdim()*A.coldim();i++)
					mm=std::max(mm,A.get(i,k).bitsize());
			return mm;
		}


	public:
		inline const IntField & field() const { return *_field; }


		PolynomialMatrixFFTMulDomain (const IntField &F, const integer maxnorm=0) :
			_field(&F) {}

		template<typename PMatrix1, typename PMatrix2, typename PMatrix3>
		void mul (PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b) {
			//compute a bound on the entry of the input matrix a and b
			FFT_PROFILE_START;
			integer maxA,maxB;
			maxA=maxB=_maxnorm;
			if (_maxnorm==0){
				maxA=1;maxA<<=(logmax(a));
				maxB=1;maxB<<=(logmax(b));
			}
			integer bound=2*maxA*maxB*a.coldim()*(std::min(a.size(),b.size())-1);
			FFT_PROFILING(2,"max norm computation");

			mul_crtla(c,a,b,maxA,maxB,bound);
		}

		template<typename PMatrix1, typename PMatrix2, typename PMatrix3>
		void midproduct (PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b,
				 bool smallLeft=true, size_t n0=0, size_t n1=0) {
			//compute a bound on the entry of the input matrix a and b
			FFT_PROFILE_START;
			integer maxA,maxB;
			maxA=maxB=_maxnorm;
			if (_maxnorm==0){
				maxA=1;maxA<<=(logmax(a));
				maxB=1;maxB<<=(logmax(b));
			}
			integer bound=2*maxA*maxB*a.coldim();
			if (smallLeft)
				bound*= a.size();
			else
				bound*= b.size();

			FFT_PROFILING(2,"max norm computation");

			midproduct_crtla(c,a,b,maxA,maxB,bound,smallLeft, n0,n1);
		}


		template< typename PMatrix1,typename PMatrix2, typename PMatrix3>
		void mul_crtla(PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b,
			       const integer& maxA, const integer& maxB, const integer& bound) {
			// (convert to MatrixP representation)
			MatrixP_I a2(field(),a.rowdim(),a.coldim(),a.size());
			MatrixP_I b2(field(),b.rowdim(),b.coldim(),b.size());
			a2.copy(a,0,a.size()-1);
			b2.copy(b,0,b.size()-1);
			MatrixP_I c2(field(),c.rowdim(),c.coldim(),c.size());
			mul_crtla(c2,a2,b2,maxA,maxB,bound);
			c.copy(c2,0,c.size()-1);
		}


		void mul_crtla(MatrixP_I &c, const MatrixP_I &a, const MatrixP_I &b,
			       const integer& maxA, const integer& maxB, const integer& bound){
			FFT_PROFILE_START;
			linbox_check(a.coldim() == b.rowdim());
			size_t m = a.rowdim();
			size_t k = a.coldim();
			size_t n = b.coldim();
			size_t s= a.size()+b.size()-1;
			size_t lpts=0;
			size_t pts  = 1; while (pts < s) { pts= pts<<1; ++lpts; }
			size_t _k=k,lk=0;
			// compute bit size of feasible prime for FFLAS
			while ( _k ) {_k>>=1; ++lk;}
			size_t prime_bitsize= (53-lk)>>1;
			RandomFFTPrime RdFFT(prime_bitsize);
			std::vector<integer> bas=RdFFT.generatePrimes(bound);
			std::vector<double> basis(bas.size());
			std::copy(bas.begin(),bas.end(),basis.begin());
			FFPACK::rns_double RNS(basis);
			size_t num_primes = RNS._size;
#ifdef FFT_PROFILER
			double tMul=0.,tCopy=0;;
			if (FFT_PROF_LEVEL<3){
				cout << "number of FFT primes :" << num_primes << endl;
				cout << "feasible prime bitsize : "<<prime_bitsize<<endl;
				cout << "bitsize of the output: "<<bound.bitsize()
				     <<"( "<< RNS._M.bitsize()<<" )"<<endl;
				cout <<" +++++++++++++++++++++++++++++++"<<endl;
			}
#endif
			FFT_PROFILING(2,"init of CRT approach");
			// reduce t_a and t_b modulo each FFT primes
			size_t n_ta=m*k*a.size(), n_tb=k*n*b.size();
			double* t_a_mod= new double[n_ta*num_primes];
			double* t_b_mod= new double[n_tb*num_primes];
			RNS.init(1, n_ta, t_a_mod, n_ta, a.getPointer(), n_ta, maxA);
			RNS.init(1, n_tb, t_b_mod, n_tb, b.getPointer(), n_tb, maxB);
			FFT_PROFILING(2,"reduction mod pi of input matrices");

			vector<MatrixP_F*> c_i (num_primes);
			vector<ModField> f(num_primes,ModField(1));
			for (size_t l=0;l<num_primes;l++)
				f[l]=ModField(RNS._basis[l]);
#ifdef HAVE_OPENMP
#pragma omp parallel for shared(c_i,f, t_a_mod,t_b_mod) schedule(dynamic)
#endif
			for (size_t l=0;l<num_primes;l++){
				FFT_PROFILE_START;
				MatrixP_F a_i (f[l], m, k, pts);
				MatrixP_F b_i (f[l], k, n, pts);
				c_i[l] = new MatrixP_F(f[l], m, n, pts);
				// copy reduced data
				for (size_t i=0;i<m*k;i++)
					for (size_t j=0;j<a.size();j++)
						a_i.ref(i,j)=t_a_mod[l*n_ta+j+i*a.size()];
				for (size_t i=0;i<k*n;i++)
					for (size_t j=0;j<b.size();j++)
						b_i.ref(i,j)=t_b_mod[l*n_tb+j+i*b.size()];
				FFT_PROFILE_GET(tCopy);
				PolynomialMatrixFFTPrimeMulDomain fftdomain (f[l]);
				fftdomain.mul_fft(lpts, *c_i[l], a_i, b_i);
				FFT_PROFILE_GET(tMul);
			}
			delete[] t_a_mod;
			delete[] t_b_mod;
			FFT_PROFILE(2,"copying linear reduced matrix",tCopy);
			FFT_PROFILE(2,"FFTprime multiplication",tMul);
			FFT_PROFILE_START;

			if (num_primes < 2) {
				FFT_PROFILE_START;
				c.copy(*c_i[0],0,s-1);
			} else {
				FFT_PROFILE_START;
				// construct contiguous storage for c_i
				double *t_c_mod;
				size_t n_tc=m*n*s;
				t_c_mod = new double[n_tc*num_primes];
				for (size_t l=0;l<num_primes;l++)
					for (size_t i=0;i<m*n;i++)
						for (size_t j=0;j<s;j++)
							t_c_mod[l*n_tc + (j+i*s)]= c_i[l]->get(i,j);
				FFT_PROFILING(2,"linearization of results mod pi");

				// reconstruct the result in C
				RNS.convert(1,n_tc,0,c.getWritePointer(),n_tc, t_c_mod, n_tc);

				delete[] t_c_mod;

			}
			FFT_PROFILING(2,"k prime reconstruction");
		}


		template< typename PMatrix1,typename PMatrix2, typename PMatrix3>
		void midproduct_crtla(PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b,
				      const integer& maxA, const integer& maxB, const integer& bound,
				      bool smallLeft=true, size_t n0=0, size_t n1=0) {
			// (convert to MatrixP representation)
			MatrixP_I a2(field(),a.rowdim(),a.coldim(),a.size());
			MatrixP_I b2(field(),b.rowdim(),b.coldim(),b.size());
			a2.copy(a,0,a.size()-1);
			b2.copy(b,0,b.size()-1);
			MatrixP_I c2(field(),c.rowdim(),c.coldim(),c.size());
			midproduct_crtla(c2,a2,b2,maxA,maxB,bound,smallLeft,n0,n1);
			c.copy(c2,0,c2.size()-1);
		}

		void midproduct_crtla(MatrixP_I &c, const MatrixP_I &a, const MatrixP_I &b,
				      const integer& maxA, const integer& maxB, const integer& bound,
				      bool smallLeft=true, size_t n0=0, size_t n1=0) {
			FFT_PROFILE_START;
			linbox_check(a.coldim() == b.rowdim());
			size_t m = a.rowdim();
			size_t k = a.coldim();
			size_t n = b.coldim();
			size_t hdeg = (n0==0?c.size():n0);
			size_t deg  = (n1==0?2*hdeg:n1);
			linbox_check(c.size()>=deg-hdeg);

			if (smallLeft){
				linbox_check(b.size()<hdeg+deg);
			}
			else
				linbox_check(a.size()<hdeg+deg);

			//linbox_check(2*c.size()-1 == b.size());
			//size_t deg= b.size()+1;
			//size_t hdeg= deg/2;
			size_t lpts=0;
			size_t pts  = 1; while (pts < deg) { pts= pts<<1; ++lpts; }
			size_t _k=k,lk=0;
			// compute bit size of feasible prime for FFLAS
			while ( _k ) {_k>>=1; ++lk;}
			size_t prime_bitsize= (53-lk)>>1;

			FFPACK::rns_double RNS(bound, prime_bitsize);
			size_t num_primes = RNS._size;
#ifdef FFT_PROFILER
			double tMul=0.,tCopy=0;;
			if (FFT_PROF_LEVEL<3){
				cout << "number of FFT primes :" << num_primes << endl;
				cout << "feasible prime bitsize : "<<prime_bitsize<<endl;
				cout << "bitsize of the output: "<<bound.bitsize()
				     <<"( "<< RNS._M.bitsize()<<" )"<<endl;
				cout <<" +++++++++++++++++++++++++++++++"<<endl;
			}
#endif
			FFT_PROFILING(2,"init of CRT approach");
			// reduce t_a and t_b modulo each FFT primes
			size_t n_ta=m*k*a.size(), n_tb=k*n*b.size();
			double* t_a_mod= new double[n_ta*num_primes];
			double* t_b_mod= new double[n_tb*num_primes];
			RNS.init(1, n_ta, t_a_mod, n_ta, a.getPointer(), n_ta, maxA);
			RNS.init(1, n_tb, t_b_mod, n_tb, b.getPointer(), n_tb, maxB);
			FFT_PROFILING(2,"reduction mod pi of input matrices");

			vector<MatrixP_F> c_i (num_primes);

			for (size_t l=0;l<num_primes;l++){
				FFT_PROFILE_START;
				ModField f(RNS._basis[l]);
				MatrixP_F a_i (f, m, k, pts);
				MatrixP_F b_i (f, k, n, pts);
				c_i[l] = MatrixP_F(f, m, n, pts);
				// copy reduced data and reversed when necessary according to midproduct algo
				for (size_t i=0;i<m*k;i++)
					for (size_t j=0;j<a.size();j++)
						if (smallLeft)
							a_i.ref(i,hdeg-1-j)=t_a_mod[l*n_ta+j+i*a.size()];
						else
							a_i.ref(i,j)=t_a_mod[l*n_ta+j+i*a.size()];
				for (size_t i=0;i<k*n;i++)
					for (size_t j=0;j<b.size();j++)
						if (smallLeft)
							b_i.ref(i,j)=t_b_mod[l*n_tb+j+i*b.size()];
						else
							b_i.ref(i,hdeg-1-j)=t_b_mod[l*n_tb+j+i*b.size()];
				FFT_PROFILE_GET(tCopy);
				PolynomialMatrixFFTPrimeMulDomain fftdomain (f);
				fftdomain.midproduct_fft(lpts, c_i[l], a_i, b_i, smallLeft);
				FFT_PROFILE_GET(tMul);
			}
			delete[] t_a_mod;
			delete[] t_b_mod;
			FFT_PROFILE(2,"copying linear reduced matrix",tCopy);
			FFT_PROFILE(2,"FFTprime multiplication",tMul);
			FFT_PROFILE_START;

			if (num_primes < 2) {
				FFT_PROFILE_START;
				c.copy(c_i[0],0,c.size()-1);
			} else {
				FFT_PROFILE_START;
				// construct contiguous storage for c_i
				double *t_c_mod;
				size_t n_tc=m*n*c.size()-1;
				t_c_mod = new double[n_tc*num_primes];
				for (size_t l=0;l<num_primes;l++)
					for (size_t i=0;i<m*n;i++)
						for (size_t j=0;j<c.size();j++)
							t_c_mod[l*n_tc + (j+i*c.size())]= c_i[l].get(i,j);
				FFT_PROFILING(2,"linearization of results mod pi");

				// reconstruct the result in C
				RNS.convert(1,n_tc,0,c.getWritePointer(),n_tc, t_c_mod, n_tc);
				delete[] t_c_mod;


				FFT_PROFILING(2,"k prime reconstruction");
			}
		}
	};


	/***************************************************************************
	 **** Polynomial Matrix Multiplication over Fp[x], with p multiprecision ***
	 ***************************************************************************/
	template <>
	class PolynomialMatrixFFTMulDomain<Givaro::Modular<integer> > {
	public:
		typedef Givaro::Modular<integer>              Field;
		typedef typename Field::Element     Element;
		typedef Givaro::UnparametricRing<integer>  IntField;
		// Polynomial matrix stored as a polynomial of matrix
		typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,Field> MatrixP_F;
		// Polynomial matrix stored as a polynomial of matrix
		typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,IntField> MatrixP_I;

	private:
		const Field            *_field;  // Read only
		integer                     _p;

	public:
		inline const Field & field() const { return *_field; }

		PolynomialMatrixFFTMulDomain(const Field &F) : _field(&F) {
			field().cardinality(_p);
		}

		template<typename Matrix1, typename Matrix2, typename Matrix3>
		void mul (Matrix1 &c, const Matrix2 &a, const Matrix3 &b) {
			MatrixP_F a2(field(),a.rowdim(),a.coldim(),a.size());
			MatrixP_F b2(field(),b.rowdim(),b.coldim(),b.size());
			a2.copy(a,0,a.size()-1);
			b2.copy(b,0,b.size()-1);
			MatrixP_F c2(field(),c.rowdim(),c.coldim(),c.size());
			mul(c2,a2,b2);
			c.copy(c2,0,c.size()-1);
		}

		void mul (MatrixP_F &c, const MatrixP_F &a, const MatrixP_F &b) {
			IntField Z;
			PolynomialMatrixFFTMulDomain<IntField> Zmul(Z,_p);
			const MatrixP_I* a2 = reinterpret_cast<const MatrixP_I*>(&a);
			const MatrixP_I* b2 = reinterpret_cast<const MatrixP_I*>(&b);
			MatrixP_I* c2       = reinterpret_cast<MatrixP_I*>(&c);
			Zmul.mul(*c2,*a2,*b2);
			// reduce the result mod p
			for (size_t i=0;i<c2->rowdim()*c2->coldim();i++)
				for (size_t j=0;j<c2->size();j++)
					c2->ref(i,j)%=_p;

		}

		template<typename Matrix1, typename Matrix2, typename Matrix3>
		void midproduct (Matrix1 &c, const Matrix2 &a, const Matrix3 &b,
				 bool smallLeft=true, size_t n0=0, size_t n1=0) {

			MatrixP_F a2(field(),a.rowdim(),a.coldim(),a.size());
			MatrixP_F b2(field(),b.rowdim(),b.coldim(),b.size());
			a2.copy(a,0,a.size()-1);
			b2.copy(b,0,b.size()-1);
			MatrixP_F c2(field(),c.rowdim(),c.coldim(),c.size());
			midproduct(c2,a2,b2,smallLeft,n0,n1);
			c.copy(c2,0,c.size()-1);
		}

		void midproduct (MatrixP_F &c, const MatrixP_F &a, const MatrixP_F &b,
				 bool smallLeft=true, size_t n0=0, size_t n1=0) {
			IntField Z;
			PolynomialMatrixFFTMulDomain<IntField> Zmul(Z,_p);
			const MatrixP_I* a2 = reinterpret_cast<const MatrixP_I*>(&a);
			const MatrixP_I* b2 = reinterpret_cast<const MatrixP_I*>(&b);
			MatrixP_I* c2       = reinterpret_cast<MatrixP_I*>(&c);
			Zmul.midproduct(*c2,*a2,*b2,smallLeft,n0,n1);
			// reduce the result mod p
			for (size_t i=0;i<c2->rowdim()*c2->coldim();i++)
				for (size_t j=0;j<c2->size();j++)
					c2->ref(i,j)%=_p;
		}
	};

} // end of namespace LinBox

#endif // __LINBOX_matpoly_mult_ftt_multiprecision_INL

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
#ifndef __LINBOX_matpoly_mult_ftt_wordsize_fast_INL
#define __LINBOX_matpoly_mult_ftt_wordsize_fast_INL

#include "givaro/modular.h"
#include "fflas-ffpack/fflas-ffpack.h"
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/algorithms/polynomial-matrix/fft.h"

namespace LinBox {

	/***********************************************************************************
	 **** Polynomial Matrix Multiplication over Zp[x] with p (FFTPrime, FFLAS prime) ***
	 ***********************************************************************************/
	template<class _Field>
	class PolynomialMatrixFFTPrimeMulDomain {

	public:
        typedef _Field Field;
		// Polynomial matrix stored as a matrix of polynomial
		typedef PolynomialMatrix<Field,PMType::polfirst> MatrixP;
		// Polynomial matrix stored as a polynomial of matrix
		typedef PolynomialMatrix<Field,PMType::matfirst> PMatrix;
        
	private:
		const Field              *_field;  // Read only
		uint64_t                      _p;
		BlasMatrixDomain<Field>     _BMD;

	public:
		inline const Field & field() const { return *_field; }

		PolynomialMatrixFFTPrimeMulDomain(const Field &F)
			: _field(&F), _p(field().cardinality()),  _BMD(F){}

		template<typename Matrix1, typename Matrix2, typename Matrix3>
        void mul (Matrix1 &c, const Matrix2 &a, const Matrix3 &b,
                                                    size_t max_rowdeg=0) const {
			FFT_PROFILE_START(1);
            linbox_check(a.coldim()==b.rowdim());
            size_t deg  = (max_rowdeg ? max_rowdeg : a.size()+b.size()-2);
            size_t lpts = 0;
            size_t pts  = 1; while (pts <= deg) { pts= pts<<1; ++lpts; }

            /* padd the input a and b to 2^lpts and convert to MatrixP
             * representation if necessary.
             */
            MatrixP a2(field(),a.rowdim(),a.coldim(),pts);
            MatrixP b2(field(),b.rowdim(),b.coldim(),pts);
            a2.copy(a,0,a.size()-1);
            b2.copy(b,0,b.size()-1);
            MatrixP c2(field(),c.rowdim(),c.coldim(),pts);
            FFT_PROFILING(1,"padding input/output to 2^k");
            mul_fft (lpts, c2, a2, b2);
            c.copy(c2,0,deg);
            FFT_PROFILING(1,"copying back the result");
        }

		// a,b and c must have size: 2^lpts
		// -> use TFT to circumvent the padding issue
		void mul_fft (size_t lpts, MatrixP &c, MatrixP &a, MatrixP &b) const {
			FFT_PROFILE_START(1);
			size_t m = a.rowdim();
			size_t k = a.coldim();
			size_t n = b.coldim();
			size_t pts=c.size();
			//std::cout<<"mul : 2^"<<lpts<<std::endl;

#ifdef CHECK_MATPOL_MUL
			MatrixP a_copy(field(),a.rowdim(),a.coldim(),pts);
			MatrixP b_copy(field(),b.rowdim(),b.coldim(),pts);
			a_copy.copy(a);
			b_copy.copy(b);
#endif

#ifdef FFT_PROFILER
			Timer totalTime;
			totalTime.start();
			if (FFT_PROF_LEVEL==1) std::cout<<"FFT: points "<<pts<<"\n";
#endif

			if ((_p-1) % pts != 0) {
				std::cout<<"Error the prime is not a FFTPrime or it has too small power of 2\n";
				std::cout<<"prime="<<_p<<std::endl;
				std::cout<<"nbr points="<<pts<<std::endl;
				throw LinboxError("LinBox ERROR: bad FFT Prime\n");
			}
            FFT<Field> FFTer(field(), lpts);
            FFT<Field> FFTinv (field(), lpts, FFTer.invroot());
            
			FFT_PROFILING(1,"init");

            // std::cout<<"FFT prime: "<<_p<<std::endl;
			// std::cout<<"FFT Root: "<<FFTer.getRoot()<<std::endl;
			// std::cout<<"FFT InvRoot: "<<FFTer.getInvRoot()<<std::endl;
			// std::cout<<a<<std::endl;
			// std::cout<<b<<std::endl;
			
			// FFT transformation on the input matrices
			for (size_t i = 0; i < m * k; i++)
				FFTer.FFT_direct(&(a.ref(i,0)));
			for (size_t i = 0; i < k * n; i++)
				FFTer.FFT_direct(&(b.ref(i,0)));
			FFT_PROFILING(1,"direct FFT_DIF");
			
			// std::cout<<"DIF:  w="<<FFTer._w<<std::endl;
			// std::cout<<a<<std::endl;
			// std::cout<<b<<std::endl;
			
			
			// convert the matrix representation to matfirst (with double coefficient)
			PMatrix vm_c (field(), m, n, pts);
#ifdef TRY1
			BlasMatrix<Field> vm_a(field(),m,k);
			BlasMatrix<Field> vm_b(field(),k,n);
			FFT_PROFILING(1,"creation of Matfirst");

			// Pointwise multiplication
			for (size_t i = 0; i < pts; ++i){
				a.setMatrix(vm_a,i);
				b.setMatrix(vm_b,i);
				_BMD.mul(vm_c[i], vm_a, vm_b);
			}
			FFT_PROFILING(1,"Pointwise mult");
			
#else
			PMatrix vm_a (field(), m, k, pts);
			PMatrix vm_b (field(), k, n, pts);
			FFT_PROFILING(1,"creation of Matfirst");
			vm_a.copy(a);
			vm_b.copy(b);
			FFT_PROFILING(1,"Polfirst to Matfirst");

			// Pointwise multiplication
			for (size_t i = 0; i < pts; ++i){
                auto vm_c_i = vm_c[i];
				_BMD.mul(vm_c_i, vm_a[i], vm_b[i]);
                vm_c.setMatrix(vm_c_i,i); // normally does nothing
            }
			FFT_PROFILING(1,"Pointwise mult");
#endif			
			// Transformation into matrix of polynomials (with int32_t coefficient)
			c.copy(vm_c);
			FFT_PROFILING(1,"Matfirst to Polfirst");

			//std::cout<<"pointwise:"<<std::endl;
			//std::cout<<c<<std::endl;			
			
			// Inverse FFT on the output matrix
			for (size_t i = 0; i < m * n; i++)
				FFTinv.FFT_inverse(&(c.ref(i,0)));
			FFT_PROFILING(1,"inverse FFT_DIT");

			// std::cout<<"DIT:"<<std::endl;
			// std::cout<<c<<std::endl;

			// Divide by pts = 2^lpts
			typename Field::Element inv_pts;
			field().init(inv_pts, pts);
			field().invin(inv_pts);
			// for (size_t i = 0; i < m * n; i++) 
			// 	for (size_t j = 0; j < pts; j++)
			// 		field().mulin(c.ref(i,j), inv_pts);
			FFLAS::fscalin(field(),c.rowdim()*c.coldim()*c.size(), inv_pts,  c.getPointer(),1);

			// std::cout<<"SCALIN:"<<std::endl;
			// std::cout<<c<<std::endl;

			FFT_PROFILING(1,"scaling the result");
#ifdef FFT_PROFILER
			totalTime.stop();
			//std::cout<<"FFT(1): total time : "<<totalTime<<std::endl;
#endif

#ifdef CHECK_MATPOL_MUL
			std::cerr<<"(Fourier prime) - "<<_p<<" - ";
			check_mul(c,a_copy,b_copy,c.size());
#endif
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

			midproduct_fft (lpts,c2, a2, b2, smallLeft);
			c.copy(c2,0,c.size()-1);
		}

		// a,b and c must have size: 2^lpts
		// -> a must have been already reversed according to the midproduct algorithm
		void midproduct_fft (size_t lpts, MatrixP &c, MatrixP &a, MatrixP &b,
				     bool smallLeft=true) const {
			FFT_PROFILE_START(1);
			size_t m = a.rowdim();
			size_t k = a.coldim();
			size_t n = b.coldim();
			size_t pts=c.size();
			//cout<<"mid : "<<pts<<endl;
#ifdef FFT_PROFILER
			if (FFT_PROF_LEVEL==1) std::cout<<"FFT: points "<<pts<<"\n";
#endif
			if ((_p-1) % pts != 0) {
				std::cout<<"Error the prime is not a FFTPrime or it has too small power of 2\n";
				std::cout<<"prime="<<_p<<std::endl;
				std::cout<<"nbr points="<<pts<<std::endl;
				throw LinboxError("LinBox ERROR: bad FFT Prime\n");
			}
			FFT<Field> FFTer (field(), lpts);
			FFT<Field> FFTinv(field(), lpts, FFTer.invroot());
			FFT_PROFILING(1,"init");

			// FFT transformation on the input matrices
			if (smallLeft){
				for (size_t i = 0; i < m * k; i++)
					FFTer.FFT_direct(&(a(i)[0]));
				for (size_t i = 0; i < k * n; i++)
					FFTinv.FFT_direct(&(b(i)[0]));
			}
			else {
				for (size_t i = 0; i < m * k; i++)
					FFTinv.FFT_direct(&(a(i)[0]));
				for (size_t i = 0; i < k * n; i++)
					FFTer.FFT_direct(&(b(i)[0]));
			}
			FFT_PROFILING(1,"direct FFT_DIF");

			// convert the matrix representation to matfirst (with double coefficient)
			PMatrix vm_c (field(), m, n, pts);
			PMatrix vm_a (field(), m, k, pts);
			PMatrix vm_b (field(), k, n, pts);
			FFT_PROFILING(1,"creation of Matfirst");
			vm_a.copy(a);
			vm_b.copy(b);
			FFT_PROFILING(1,"Polfirst to Matfirst");

			// Pointwise multiplication
			for (size_t i = 0; i < pts; ++i){
                auto vm_c_i = vm_c[i];
				_BMD.mul(vm_c_i, vm_a[i], vm_b[i]);
                vm_c.setMatrix(vm_c_i,i); // normally does nothing
            }
			FFT_PROFILING(1,"pointwise mult");

			// Transformation into matrix of polynomials (with int32_t coefficient)
			c.copy(vm_c);
			FFT_PROFILING(1,"Matfirst to Polfirst");

			// Inverse FFT on the output matrix
			for (size_t i = 0; i < m * n; i++)
				FFTer.FFT_inverse(&(c(i)[0]));
			FFT_PROFILING(1,"inverse FFT_DIT");

			// Divide by pts = 2^ltps
			typename Field::Element inv_pts;
			field().init(inv_pts, pts);
			field().invin(inv_pts);
			// for (size_t i = 0; i < m * n; i++)
			// 	for (size_t j = 0; j < pts; j++)
			// 		field().mulin(c.ref(i,j), inv_pts);
			FFLAS::fscalin(field(),c.rowdim()*c.coldim()*c.size(), inv_pts,  c.getPointer(),1);
			FFT_PROFILING(1,"scaling the result");
		}
	}; // end of class special FFT mul domain



}//end of namespace LinBox

#endif // __LINBOX_matpoly_mult_ftt_wordsize_fast_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

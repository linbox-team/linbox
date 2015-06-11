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

#ifndef __LINBOX_matpoly_mult_ftt_H
#define __LINBOX_matpoly_mult_ftt_H

#include "linbox/util/error.h"
#include "linbox/util/debug.h"
#include "linbox/util/timer.h"

#include "linbox/integer.h"
#include <givaro/zring.h>
#include "linbox/field/modular.h"
#include "givaro/givtimer.h"

#ifdef FFT_PROFILER
#include <iostream>
#ifndef FFT_PROF_LEVEL
#define  FFT_PROF_LEVEL 1
#endif

//size_t  FFT_PROF_LEVEL=1;
Givaro::Timer mychrono;
#define FFT_PROF_MSG_SIZE 35
#define FFT_PROFILE_START  mychrono.clear();mychrono.start();

#define FFT_PROFILING(lvl,msg)						\
if (lvl>=FFT_PROF_LEVEL) {					\
  mychrono.stop();std::cout<<"FFT: ";					\
	std::cout.width(FFT_PROF_MSG_SIZE);std::cout<<std::left<<msg<<" : "; \
	std::cout.precision(6);std::cout<<mychrono.usertime()<<" s (real: "<<mychrono.realtime()<<")"<<std::endl; \
	mychrono.clear();mychrono.start();				\
}
#ifdef HAVE_OPENMP								
#define FFT_PROFILE_GET(x)			\
  mychrono.stop();(x)+=mychrono.realtime();mychrono.clear();mychrono.start();
#else
#define FFT_PROFILE_GET(x)					\
  mychrono.stop();(x)+=mychrono.usertime();mychrono.clear();mychrono.start();
#endif
#define FFT_PROFILE(lvl,msg,x)						\
if ((lvl)>=FFT_PROF_LEVEL) {					\
  std::cout<<"FFT: ";						   \
  std::cout.width(FFT_PROF_MSG_SIZE);std::cout<<std::left<<msg<<" : ";	\
  std::cout.precision(6);std::cout<<x<<" s"<<std::endl;			\
}
#else
#define FFT_PROFILE_START
#define FFT_PROFILING(lvl,msg)
#define FFT_PROFILE_GET(x)
#define FFT_PROFILE(lvl,msg,x)
#endif // FFT_PROFILER


#ifndef FFT_DEG_THRESHOLD   
#define FFT_DEG_THRESHOLD   4
#endif

namespace LinBox
{
	// generic handler for multiplication using FFT
	template <class Field>
	class PolynomialMatrixFFTMulDomain {
	public:
		inline const Field & field() const;

		PolynomialMatrixFFTMulDomain (const Field& F);

		template<typename Matrix1, typename Matrix2, typename Matrix3>
		void mul (Matrix1 &c, const Matrix2 &a, const Matrix3 &b);

		template<typename Matrix1, typename Matrix2, typename Matrix3>
		void midproduct (Matrix1 &c, const Matrix2 &a, const Matrix3 &b, bool smallLeft=true, size_t n0=0,size_t n1=0);
	};
		
	
	//class PolynomialMatrixFFTPrimeMulDomain ;                         // Mul in Zp[x] with p <2^32, (fflas, fourier)
		
	// template <class T>
	// class PolynomialMatrixFFTMulDomain<Givaro::Modular<T> > ;        // Mul in Zp[x] with p^2 storable in type T

	// template<>
	// class PolynomialMatrixFFTMulDomain<Givaro::ZRing<integer> >;  // Mul in Z[x]

	// template <>
	// class PolynomialMatrixFFTMulDomain<Givaro::Modular<integer> > ;           // Mul in Zp[x] with p multiprecision

} // end of namespace LinBox

#include "linbox/algorithms/polynomial-matrix/matpoly-mult-fft-wordsize-fast.inl"
#include "linbox/algorithms/polynomial-matrix/matpoly-mult-fft-wordsize-three-primes.inl"
#include "linbox/algorithms/polynomial-matrix/matpoly-mult-fft-multiprecision.inl"
#include "linbox/algorithms/polynomial-matrix/matpoly-mult-fft-wordsize.inl"

#endif // __LINBOX_matpoly_mult_ftt_H

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
#include "linbox/field/unparametric.h"
#include "linbox/field/modular.h"

#ifdef FFT_PROFILER	
#include <iostream>
using namespace std;
size_t  FFT_PROF_LEVEL=1; 
myTimer chrono;				
#define FFT_PROF_MSG_SIZE 35						
#define FFT_PROFILE_START  chrono.start();
 
#define FFT_PROFILING(lvl,msg)						\
	if (lvl>=FFT_PROF_LEVEL) {					\
		chrono.stop();cout<<"FFT: ";				\
		cout.width(FFT_PROF_MSG_SIZE);cout<<std::left<<msg<<" : "; \
		cout.precision(6);cout<<chrono.gettime()<<" s"<<endl;	\
		chrono.start();						\
	}
#define FFT_PROFILE_GET(x)					\
	chrono.stop();x+=chrono.gettime();chrono.start();

#define FFT_PROFILE(lvl,msg,x)						\
	if (lvl>=FFT_PROF_LEVEL) {					\
		cout<<"FFT: ";						\
		cout.width(FFT_PROF_MSG_SIZE);cout<<std::left<<msg<<" : "; \
		cout.precision(6);cout<<x<<" s"<<endl;			\
	}
#else
#define FFT_PROFILE_START			
#define FFT_PROFILING(lvl,msg)
#define FFT_PROFILE_GET(x)
#define FFT_PROFILE(lvl,msg,x)						
#endif						
	

namespace LinBox
{
		
	template <class _Field>
	class PolynomialMatrixFFTMulDomain ;                        // generic handler for multiplication using FFT
	template <>
	class PolynomialMatrixFFTMulDomain<Modular<uint32_t> > ;   // Mul in Zp[x] with p <2^32
	class PolynomialMatrixFFTPrimeMulDomain ;                  // Mul in Zp[x] with p FFTPrime and FFLAS
	
	//template<>
	//class PolynomialMatrixFFTMulDomain<FFPACK::UnparametricField<integer> >;           // Mul in Z[x]
	//template <>
	//class PolynomialMatrixFFTMulDomain<FFPACK::Modular<integer> > ;    // Mul in Zp[x] with p multiprecision

} // end of namespace LinBox

#include "linbox/algorithms/polynomial-matrix/matpoly-mult-fft-wordsize.inl"
//#include "linbox/algorithms/polynomial-matrix/matpoly-mult-fft-multiprecision.inl"


#endif

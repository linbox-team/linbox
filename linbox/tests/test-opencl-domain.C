/* tests/test-blas-domain.C
 * Copyright (C) 2011 Matthew Wezowicz
 *
 * Written by Matthew Wezowicz <mwezz@udel.edu>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file  tests/test-blas-domain.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */

#include <string>
#include <iostream>

#include "linbox/field/modular.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/algorithms/opencl-domain.h"
#include "linbox/randiter/nonzero.h"
#include "linbox/util/commentator.h"

#include "test-common.h"

#ifndef _OPENCL_TIMER_
#define _OPENCL_TIMER_

#include <stack>
#include <vector>
#include <sys/time.h>

class OpenCLTimer{

private:
	std::stack<struct timeval, std::vector<struct timeval> > timerStack;

public:
	OpenCLTimer(){}

	~OpenCLTimer(){}

	void tic(){
		struct timeval tempTimeVal;
		gettimeofday(&tempTimeVal, NULL);

		timerStack.push(tempTimeVal);
	}

	double toc(){
		double returnValue;
		struct timeval now;
		struct timeval then;

		gettimeofday(&now, NULL);

		then = timerStack.top();
		timerStack.pop();

		returnValue = ((double)(now.tv_sec - then.tv_sec) +
		              ((double)(now.tv_usec - then.tv_usec) * 1e-6));

		return returnValue;
	}

};

#endif

using namespace LinBox;

const int maxpretty = 35;

const char* pretty(std::string a) {
	std::string blank;
	blank = a;
	int msgsize= maxpretty - (int)blank.size();
	std::string dot(".");
	for(int i=0;i<msgsize ;++i){
		blank += dot;
	}
	return blank.c_str();
}

template <class Field>
static bool testMul(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Field> Matrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing mul"),"testMul",(size_t)iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator().progress(i);
		Matrix A(F,n,n);
		Matrix B(F,n,n);
		Matrix C_b(F,n,n);
		Matrix C_o(F,n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
			}
		}

		OpenCLTimer timer;
		timer.tic();
		BMD.mul(C_b,A,B);
		commentator().report() << "Blas mul() execution time: "
		                       << timer.toc() << std::endl;
		timer.tic();
		OMD.mul(C_o,A,B);
		commentator().report() << "OpenCL mul() execution time: "
		                       << timer.toc() << std::endl;

		if(!OMD.areEqual(C_b,C_o)){
			ret = false;
		}
	}

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMul");

	return ret;
}

template <class Field>
static bool testMulinLeft(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Field> Matrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing mulin_left"),"testMulinLeft",(size_t)iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator().progress(i);
		Matrix A_b(F,n,n);
		Matrix A_o(F,n,n);
		Matrix B(F,n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A_b.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
			}
		}

		A_o = A_b;

		OpenCLTimer timer;
		timer.tic();
		BMD.mulin_left(A_b,B);
		commentator().report() << "Blas mulin_left() execution time: "
		                       << timer.toc() << std::endl;
		timer.tic();
		OMD.mulin_left(A_o,B);
		commentator().report() << "OpenCL mulin_left() execution time: "
		                       << timer.toc() << std::endl;

		if(!OMD.areEqual(A_b,A_o)){
			ret = false;
		}
	}

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMulinLeft");

	return ret;
}

template <class Field>
static bool testMulinRight(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Field> Matrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing mulin_right"),"testMulinRight",(size_t)iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator().progress(i);
		Matrix A(F,n,n);
		Matrix B_b(F,n,n);
		Matrix B_o(F,n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B_b.setEntry(k,j,G.random(tmp));
			}
		}

		B_o = B_b;

		OpenCLTimer timer;
		timer.tic();
		BMD.mulin_right(A,B_b);
		commentator().report() << "Blas mulin_right() execution time: "
		                       << timer.toc() << std::endl;
		timer.tic();
		OMD.mulin_right(A,B_o);
		commentator().report() << "OpenCL mulin_right() execution time: "
		                       << timer.toc() << std::endl;

		if(!OMD.areEqual(B_b,B_o)){
			ret = false;
		}
	}

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMulinRight");

	return ret;
}

template <class Field>
static bool testAxpy(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Field> Matrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing axpy"),"testAxpy",(size_t)iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator().progress(i);
		Matrix A(F,n,n);
		Matrix B(F,n,n);
		Matrix C(F,n,n);
		Matrix D_b(F,n,n);
		Matrix D_o(F,n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
				C.setEntry(k,j,G.random(tmp));
			}
		}

		OpenCLTimer timer;
		timer.tic();
		BMD.axpy(D_b,A,B,C);
		commentator().report() << "Blas axpy() execution time: "
		                       << timer.toc() << std::endl;
		timer.tic();
		OMD.axpy(D_o,A,B,C);
		commentator().report() << "OpenCL axpy() execution time: "
		                       << timer.toc() << std::endl;

		if(!OMD.areEqual(D_b,D_o)){
			ret = false;
		}
	}

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testAxpy");

	return ret;
}

template <class Field>
static bool testAxpyin(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Field> Matrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing axpyin"),"testAxpyin",(size_t)iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator().progress(i);
		Matrix A(F,n,n);
		Matrix B(F,n,n);
		Matrix C_b(F,n,n);
		Matrix C_o(F,n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
				C_b.setEntry(k,j,G.random(tmp));
			}
		}

		C_o = C_b;

		OpenCLTimer timer;
		timer.tic();
		BMD.axpyin(C_b,A,B);
		commentator().report() << "Blas axpyin() execution time: "
		                       << timer.toc() << std::endl;
		timer.tic();
		OMD.axpyin(C_o,A,B);
		commentator().report() << "OpenCL axpyin() execution time: "
		                       << timer.toc() << std::endl;

		if(!OMD.areEqual(C_b,C_o)){
			ret = false;
		}
	}

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testAxpyin");

	return ret;
}

template <class Field>
static bool testMaxpy(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Field> Matrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing maxpy"),"testMaxpy",(size_t)iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator().progress(i);
		Matrix A(F,n,n);
		Matrix B(F,n,n);
		Matrix C(F,n,n);
		Matrix D_b(F,n,n);
		Matrix D_o(F,n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
				C.setEntry(k,j,G.random(tmp));
			}
		}

		OpenCLTimer timer;
		timer.tic();
		BMD.maxpy(D_b,A,B,C);
		commentator().report() << "Blas maxpy() execution time: "
		                       << timer.toc() << std::endl;
		timer.tic();
		OMD.maxpy(D_o,A,B,C);
		commentator().report() << "OpenCL maxpy() execution time: "
		                       << timer.toc() << std::endl;

		if(!OMD.areEqual(D_b,D_o)){
			ret = false;
		}
	}

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMaxpy");

	return ret;
}

template <class Field>
static bool testMaxpyin(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Field> Matrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing maxpyin"),"testMaxpyin",(size_t)iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator().progress(i);
		Matrix A(F,n,n);
		Matrix B(F,n,n);
		Matrix C_b(F,n,n);
		Matrix C_o(F,n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
				C_b.setEntry(k,j,G.random(tmp));
			}
		}

		C_o = C_b;

		OpenCLTimer timer;
		timer.tic();
		BMD.maxpyin(C_b,A,B);
		commentator().report() << "Blas maxpyin() execution time: "
		                       << timer.toc() << std::endl;
		timer.tic();
		OMD.maxpyin(C_o,A,B);
		commentator().report() << "OpenCL maxpyin() execution time: "
		                       << timer.toc() << std::endl;

		if(!OMD.areEqual(C_b,C_o)){
			ret = false;
		}
	}

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMaxpyin");

	return ret;
}

template <class Field>
static bool testAxmy(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Field> Matrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing axmy"),"testAxmy",(size_t)iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator().progress(i);
		Matrix A(F,n,n);
		Matrix B(F,n,n);
		Matrix C(F,n,n);
		Matrix D_b(F,n,n);
		Matrix D_o(F,n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
				C.setEntry(k,j,G.random(tmp));
			}
		}

		OpenCLTimer timer;
		timer.tic();
		BMD.axmy(D_b,A,B,C);
		commentator().report() << "Blas axmy() execution time: "
		                       << timer.toc() << std::endl;
		timer.tic();
		OMD.axmy(D_o,A,B,C);
		commentator().report() << "OpenCL axmy() execution time: "
		                       << timer.toc() << std::endl;

		if(!OMD.areEqual(D_b,D_o)){
			ret = false;
		}
	}

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testAxmy");

	return ret;
}

template <class Field>
static bool testAxmyin(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Field> Matrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing axmyin"),"testAxmyin",(size_t)iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator().progress(i);
		Matrix A(F,n,n);
		Matrix B(F,n,n);
		Matrix C_b(F,n,n);
		Matrix C_o(F,n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
				C_b.setEntry(k,j,G.random(tmp));
			}
		}

		C_o = C_b;

		OpenCLTimer timer;
		timer.tic();
		BMD.axmyin(C_b,A,B);
		commentator().report() << "Blas axmyin() execution time: "
		                       << timer.toc() << std::endl;
		timer.tic();
		OMD.axmyin(C_o,A,B);
		commentator().report() << "OpenCL axmyin() execution time: "
		                       << timer.toc() << std::endl;

		if(!OMD.areEqual(C_b,C_o)){
			ret = false;
		}
	}

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testAxmyin");

	return ret;
}

template <class Field>
static bool testMuladd(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Field> Matrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing muladd"),"testMuladd",(size_t)iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator().progress(i);
		Matrix A(F,n,n);
		Matrix B(F,n,n);
		Matrix C(F,n,n);
		Matrix D_b(F,n,n);
		Matrix D_o(F,n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
				C.setEntry(k,j,G.random(tmp));
			}
		}

		OpenCLTimer timer;
		timer.tic();
		BMD.muladd(D_b,1.0,C,2.0,A,B);
		commentator().report() << "Blas muladd() execution time: "
		                       << timer.toc() << std::endl;
		timer.tic();
		OMD.muladd(D_o,1.0,C,2.0,A,B);
		commentator().report() << "OpenCL muladd() execution time: "
		                       << timer.toc() << std::endl;

		if(!OMD.areEqual(D_b,D_o)){
			ret = false;
		}

		/*
		BMD.muladd(D_b,-1.0,C,2.0,A,B);
		OMD.muladd(D_o,-1.0,C,2.0,A,B);

		if(!OMD.areEqual(D_b,D_o)){
			ret = false;
		}

		BMD.muladd(D_b,1.0,C,-2.0,A,B);
		OMD.muladd(D_o,1.0,C,-2.0,A,B);

		if(!OMD.areEqual(D_b,D_o)){
			ret = false;
		}

		BMD.muladd(D_b,-1.0,C,-2.0,A,B);
		OMD.muladd(D_o,-1.0,C,-2.0,A,B);

		if(!OMD.areEqual(D_b,D_o)){
			ret = false;
		}
		*/
	}

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMuladd");

	return ret;
}

template <class Field>
static bool testMuladdin(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Field> Matrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing muladdin"),"testMuladdin",(size_t)iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator().progress(i);
		Matrix A(F,n,n);
		Matrix B(F,n,n);
		Matrix C_b(F,n,n);
		Matrix C_o(F,n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
				C_b.setEntry(k,j,G.random(tmp));
			}
		}

		C_o = C_b;

		OpenCLTimer timer;
		timer.tic();
		BMD.muladdin(1.0,C_b,2.0,A,B);
		commentator().report() << "Blas muladdin() execution time: "
		                       << timer.toc() << std::endl;
		timer.tic();
		OMD.muladdin(1.0,C_o,2.0,A,B);
		commentator().report() << "OpenCL muladdin() execution time: "
		                       << timer.toc() << std::endl;

		if(!OMD.areEqual(C_b,C_o)){
			ret = false;
		}
	}

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMuladdin");

	return ret;
}

template <class Field>
static bool testMulscale(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Field> Matrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing mulscale"),"testMulscale",(size_t)iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator().progress(i);
		Matrix A(F,n,n);
		Matrix B(F,n,n);
		Matrix C_b(F,n,n);
		Matrix C_o(F,n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
			}
		}

		OpenCLTimer timer;
		timer.tic();
		BMD.mul(C_b,2.0,A,B);
		commentator().report() << "Blas mulscale() execution time: "
		                       << timer.toc() << std::endl;
		timer.tic();
		OMD.mul(C_o,2.0,A,B);
		commentator().report() << "OpenCL mulscale() execution time: "
		                       << timer.toc() << std::endl;

		if(!OMD.areEqual(C_b,C_o)){
			ret = false;
		}
	}

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMulscale");

	return ret;
}

template <class Field>
int launch_tests(Field& F, int n, int iterations){
	bool pass = true;

	if(!testMul(F, (size_t)n, iterations)){
		pass = false;
	}
	// if(!testMulinLeft(F, n, iterations)){
		// pass = false;
	// }
	// if(!testMulinRight(F, n, iterations)){
		// pass = false;
	// }
	if(!testAxpy(F, (size_t)n, iterations)){
		pass = false;
	}
	// if(!testAxpyin(F, n, iterations)){
		// pass = false;
	// }
	if(!testMaxpy(F, (size_t)n, iterations)){
		pass = false;
	}
	// if(!testMaxpyin(F, n, iterations)){
		// pass = false;
	// }
	if(!testAxmy(F, (size_t)n, iterations)){
		pass = false;
	}
	// if(!testAxmyin(F, n, iterations)){
		// pass = false;
	// }
	if(!testMuladd(F, (size_t)n, iterations)){
		pass = false;
	}
	// if(!testMuladdin(F, n, iterations)){
		// pass = false;
	// }
	// if(!testMulscale(F, n, iterations)){
		// pass = false;
	// }

	return pass;
}

/*
void* performOpenCL_mul(void* t){
	Modular<double> F(1000003U);

	OpenCLMatrixDomain<Modular<double> > OMD(F);

	Modular<double>::RandIter G(F);

	BlasMatrix<Modular<double> > A(F,500,500);
	BlasMatrix<Modular<double> > B(F,500,500);
	BlasMatrix<Modular<double> > C(F,500,500);

	double temp;

	for(size_t k = 0; k < 500; k++){
		for(size_t j = 0; j < 500; j++){
			A.setEntry(k,j,G.random(temp));
			B.setEntry(k,j,G.random(temp));
		}
	}

	OMD.mul(C,A,B);

	pthread_exit(NULL);
}

void testThreadSafe(int iterations){
	pthread_t* pthreads = (pthread_t*)operator new(iterations * sizeof(pthread_t));

	for(int i = 0; i < iterations; i++){
		pthread_create(&pthreads[i], NULL, performOpenCL_mul, NULL);
	}

	for(int i = 0; i < iterations; i++){
		pthread_join(pthreads[i], NULL);
	}
}
*/

int main(int argc, char** argv){
	static size_t n = 500;
	static int q = 1000003U;
	static int iterations = 3;

	static Argument args[] = {
		{'n', "-n N", "Set dimension of test matrices to NxN", TYPE_INT, &n},
		{'q', "-q Q", "Operate over the \"field\" GF(Q) [1]", TYPE_INT, &q},
		{ 'i', "-i I", "Perform each test for I iterations", TYPE_INT, &iterations},
		END_OF_ARGUMENTS};

	parseArguments (argc, argv, args);

	Modular<double> F(1000003U);
	Modular<float> H(3);

	bool pass = true;

	srand((unsigned)time(NULL));

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);

	commentator().start("OpenCLMatrixDomain test suite", "OpenCLMatrixDomain");

	//For warmup of OpenCLMatrixDomainFactory
	OpenCLMatrixDomain<Modular<double> > OMD(F);
	commentator().report() << "Number of OpenCL Devices: "
#ifdef __LINBOX_HAVE_OCL
	                       << OpenCLMatrixDomainFactory::oclGetNumberOfDevices()
#else
	                       << 0
#endif
	                       << std::endl;

	//testThreadSafe(iterations * 5);

	pass &= launch_tests(F, (int)n, iterations);
	pass &= launch_tests(H, (int)n, iterations);

	commentator().stop(MSG_STATUS(pass), (const char*)0, "OpenCLMatrixDomain test suite");
	return (pass ? 0 : -1);
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

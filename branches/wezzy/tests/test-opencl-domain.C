/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* tests/test-blas-domain.C
 * Copyright (C) 2011 Matthew Wezowicz
 *
 * Written by Matthew Wezowicz <mwezz@udel.edu>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 *
 */

/*! @file  tests/test-blas-domain.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */

#include <string>
#include <iostream>
#include <pthread.h>

#include <linbox/field/modular.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/algorithms/blas-domain.h>
#include <linbox/algorithms/opencl-domain.h>
#include <linbox/randiter/nonzero.h>
#include <linbox/util/commentator.h>

#include "test-common.h"

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
	typedef BlasMatrix<Element> Matrix;

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.start(pretty("Testing mul"),"testMul",iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator.progress(i);
		Matrix A(n,n);
		Matrix B(n,n);
		Matrix C_b(n,n);
		Matrix C_o(n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
			}
		}

		BMD.mul(C_b,A,B);
		OMD.mul(C_o,A,B);

		if(!OMD.areEqual(C_b,C_o)){
			ret = false;
		}
	}

	commentator.stop(MSG_STATUS(ret), (const char*)0, "testMul");

	return ret;
}

template <class Field>
static bool testMulinLeft(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Element> Matrix;

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.start(pretty("Testing mulin_left"),"testMulinLeft",iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator.progress(i);
		Matrix A_b(n,n);
		Matrix A_o(n,n);
		Matrix B(n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A_b.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
			}
		}

		A_o = A_b;

		BMD.mulin_left(A_b,B);
		OMD.mulin_left(A_o,B);

		if(!OMD.areEqual(A_b,A_o)){
			ret = false;
		}
	}

	commentator.stop(MSG_STATUS(ret), (const char*)0, "testMulinLeft");

	return ret;
}

template <class Field>
static bool testMulinRight(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Element> Matrix;

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.start(pretty("Testing mulin_right"),"testMulinRight",iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator.progress(i);
		Matrix A(n,n);
		Matrix B_b(n,n);
		Matrix B_o(n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B_b.setEntry(k,j,G.random(tmp));
			}
		}

		B_o = B_b;

		BMD.mulin_right(A,B_b);
		OMD.mulin_right(A,B_o);

		if(!OMD.areEqual(B_b,B_o)){
			ret = false;
		}
	}

	commentator.stop(MSG_STATUS(ret), (const char*)0, "testMulinRight");

	return ret;
}

template <class Field>
static bool testAxpy(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Element> Matrix;

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.start(pretty("Testing axpy"),"testAxpy",iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator.progress(i);
		Matrix A(n,n);
		Matrix B(n,n);
		Matrix C(n,n);
		Matrix D_b(n,n);
		Matrix D_o(n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
				C.setEntry(k,j,G.random(tmp));
			}
		}

		BMD.axpy(D_b,A,B,C);
		OMD.axpy(D_o,A,B,C);

		if(!OMD.areEqual(D_b,D_o)){
			ret = false;
		}
	}

	commentator.stop(MSG_STATUS(ret), (const char*)0, "testAxpy");

	return ret;
}

template <class Field>
static bool testAxpyin(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Element> Matrix;

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.start(pretty("Testing axpyin"),"testAxpyin",iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator.progress(i);
		Matrix A(n,n);
		Matrix B(n,n);
		Matrix C_b(n,n);
		Matrix C_o(n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
				C_b.setEntry(k,j,G.random(tmp));
			}
		}

		C_o = C_b;

		BMD.axpyin(C_b,A,B);
		OMD.axpyin(C_o,A,B);

		if(!OMD.areEqual(C_b,C_o)){
			ret = false;
		}
	}

	commentator.stop(MSG_STATUS(ret), (const char*)0, "testAxpyin");

	return ret;
}

template <class Field>
static bool testMaxpy(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Element> Matrix;

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.start(pretty("Testing maxpy"),"testMaxpy",iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator.progress(i);
		Matrix A(n,n);
		Matrix B(n,n);
		Matrix C(n,n);
		Matrix D_b(n,n);
		Matrix D_o(n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
				C.setEntry(k,j,G.random(tmp));
			}
		}

		BMD.maxpy(D_b,A,B,C);
		OMD.maxpy(D_o,A,B,C);

		if(!OMD.areEqual(D_b,D_o)){
			ret = false;
		}
	}

	commentator.stop(MSG_STATUS(ret), (const char*)0, "testMaxpy");

	return ret;
}

template <class Field>
static bool testMaxpyin(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Element> Matrix;

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.start(pretty("Testing maxpyin"),"testMaxpyin",iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator.progress(i);
		Matrix A(n,n);
		Matrix B(n,n);
		Matrix C_b(n,n);
		Matrix C_o(n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
				C_b.setEntry(k,j,G.random(tmp));
			}
		}

		C_o = C_b;

		BMD.maxpyin(C_b,A,B);
		OMD.maxpyin(C_o,A,B);

		if(!OMD.areEqual(C_b,C_o)){
			ret = false;
		}
	}

	commentator.stop(MSG_STATUS(ret), (const char*)0, "testMaxpyin");

	return ret;
}

template <class Field>
static bool testAxmy(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Element> Matrix;

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.start(pretty("Testing axmy"),"testAxmy",iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator.progress(i);
		Matrix A(n,n);
		Matrix B(n,n);
		Matrix C(n,n);
		Matrix D_b(n,n);
		Matrix D_o(n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
				C.setEntry(k,j,G.random(tmp));
			}
		}

		BMD.axmy(D_b,A,B,C);
		OMD.axmy(D_o,A,B,C);

		if(!OMD.areEqual(D_b,D_o)){
			ret = false;
		}
	}

	commentator.stop(MSG_STATUS(ret), (const char*)0, "testAxmy");

	return ret;
}

template <class Field>
static bool testAxmyin(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Element> Matrix;

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.start(pretty("Testing axmyin"),"testAxmyin",iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator.progress(i);
		Matrix A(n,n);
		Matrix B(n,n);
		Matrix C_b(n,n);
		Matrix C_o(n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
				C_b.setEntry(k,j,G.random(tmp));
			}
		}

		C_o = C_b;

		BMD.axmyin(C_b,A,B);
		OMD.axmyin(C_o,A,B);

		if(!OMD.areEqual(C_b,C_o)){
			ret = false;
		}
	}

	commentator.stop(MSG_STATUS(ret), (const char*)0, "testAxmyin");

	return ret;
}

template <class Field>
static bool testMuladd(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Element> Matrix;

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.start(pretty("Testing muladd"),"testMuladd",iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator.progress(i);
		Matrix A(n,n);
		Matrix B(n,n);
		Matrix C(n,n);
		Matrix D_b(n,n);
		Matrix D_o(n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
				C.setEntry(k,j,G.random(tmp));
			}
		}

		BMD.muladd(D_b,1.0,C,2.0,A,B);
		OMD.muladd(D_o,1.0,C,2.0,A,B);

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

	commentator.stop(MSG_STATUS(ret), (const char*)0, "testMuladd");

	return ret;
}

template <class Field>
static bool testMuladdin(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Element> Matrix;

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.start(pretty("Testing muladdin"),"testMuladdin",iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator.progress(i);
		Matrix A(n,n);
		Matrix B(n,n);
		Matrix C_b(n,n);
		Matrix C_o(n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
				C_b.setEntry(k,j,G.random(tmp));
			}
		}

		C_o = C_b;

		BMD.muladdin(1.0,C_b,2.0,A,B);
		OMD.muladdin(1.0,C_o,2.0,A,B);

		if(!OMD.areEqual(C_b,C_o)){
			ret = false;
		}
	}

	commentator.stop(MSG_STATUS(ret), (const char*)0, "testMuladdin");

	return ret;
}

template <class Field>
static bool testMulscale(const Field& F, size_t n, int iterations){
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef BlasMatrix<Element> Matrix;

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.start(pretty("Testing mulscale"),"testMulscale",iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	OpenCLMatrixDomain<Field> OMD(F);

	for(int i = 0; i < iterations; i++){
		commentator.progress(i);
		Matrix A(n,n);
		Matrix B(n,n);
		Matrix C_b(n,n);
		Matrix C_o(n,n);

		Element tmp;

		for(size_t k = 0; k < n; k++){
			for(size_t j = 0; j < n; j++){
				A.setEntry(k,j,G.random(tmp));
				B.setEntry(k,j,G.random(tmp));
			}
		}

		BMD.mul(C_b,2.0,A,B);
		OMD.mul(C_o,2.0,A,B);

		if(!OMD.areEqual(C_b,C_o)){
			ret = false;
		}
	}

	commentator.stop(MSG_STATUS(ret), (const char*)0, "testMulscale");

	return ret;
}

template <class Field>
int launch_tests(Field& F, int n, int iterations){
	bool pass = true;

	if(!testMul(F, n, iterations)){
		pass = false;
	}
	if(!testMulinLeft(F, n, iterations)){
		pass = false;
	}
	if(!testMulinRight(F, n, iterations)){
		pass = false;
	}
	if(!testAxpy(F, n, iterations)){
		pass = false;
	}
	if(!testAxpyin(F, n, iterations)){
		pass = false;
	}
	if(!testMaxpy(F, n, iterations)){
		pass = false;
	}
	if(!testMaxpyin(F, n, iterations)){
		pass = false;
	}
	if(!testAxmy(F, n, iterations)){
		pass = false;
	}
	if(!testAxmyin(F, n, iterations)){
		pass = false;
	}
	if(!testMuladd(F, n, iterations)){
		pass = false;
	}
	if(!testMuladdin(F, n, iterations)){
		pass = false;
	}
	if(!testMulscale(F, n, iterations)){
		pass = false;
	}

	return pass;
}

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
	Modular<float> H(2011);

	bool pass = true;

	srand((unsigned)time(NULL));

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);

	commentator.start("OpenCLMatrixDomain test suite", "OpenCLMatrixDomain");

	//For warmup of OpenCLMatrixDomainFactory
	OpenCLMatrixDomain<Modular<double> > OMD(F);

	pass &= launch_tests(F, (int)n, iterations);
	pass &= launch_tests(H, (int)n, iterations);

	commentator.stop(MSG_STATUS(pass), (const char*)0, "OpenCLMatrixDomain test suite");
	return (pass ? 0 : -1);
}
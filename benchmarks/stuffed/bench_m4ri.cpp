#include <iostream>
#include <chrono>
#include <limits.h>
#include <array>
#include <stdlib.h>
#include <string.h>
#include <bitset>
#include <time.h>

using namespace std::chrono;

//#include "compressed-matrix-domain.h"
#include "/usa/saunders/sandbox/m4ri/m4ri/config.h"
#include "/usa/saunders/sandbox/m4ri/m4ri/m4ri.h"

typedef high_resolution_clock::time_point my_time_t;

double gffops_mul (my_time_t t1, my_time_t t2, size_t dim, size_t trials=1){
	auto mics1 = duration_cast<nanoseconds>( t2 - t1 ).count();
	double ms1 = double(mics1);
	double ops = 2.0f*dim*dim*dim*trials;
	return ops / ms1;
}

//template <unsigned int P, typename T, unsigned int N, unsigned int R, unsigned int L, unsigned int D, unsigned int NT, unsigned int K>
double bench_m4ri(size_t dim, size_t trials)
{
	size_t m,l,n; m = l = n = dim;

	/* we create two random matrices */
	mzd_t *A = mzd_init(m, l);
    mzd_t *B = mzd_init(l, n);
	mzd_randomize(A);
	mzd_randomize(B);
	//using Field = CompressedField<P,T,N,R,L,D>;
	//Field field{};
	//CompressedMatrixDomain<Field> CMD(field);
	//using Matrix = typename CompressedMatrixDomain<Field>::Matrix;
	//Matrix A(field,dim,dim);
	//Matrix B(field,dim,dim);
	//Matrix C(field,dim,dim);
	//A.random();
	//B.random();

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (volatile size_t i = 0; i < trials; ++i) {
		/* C = A*B via M4RM, temporary buffers are managed internally */
		mzd_t *C = mzd_mul_m4rm(    NULL, A, B, 0);
		//CMD.mul_classical(C, B, A);
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	return gffops_mul (t1, t2, dim, trials);
}

#if 0
template <unsigned int P, typename T, unsigned int N, unsigned int R, unsigned int L, unsigned int D, unsigned int NT, unsigned int K>
double bench_classical(size_t dim, size_t trials)
{
	using Field = CompressedField<P,T,N,R,L,D>;
	Field field{};
	CompressedMatrixDomain<Field> CMD(field);
	using Matrix = typename CompressedMatrixDomain<Field>::Matrix;
	Matrix A(field,dim,dim);
	Matrix B(field,dim,dim);
	Matrix C(field,dim,dim);
	A.random();
	B.random();

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (volatile size_t i = 0; i < trials; ++i) {
		CMD.mul_classical(C, B, A);
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	return gffops_mul (t1, t2, dim, trials);
}

template <unsigned int P, typename T, unsigned int N, unsigned int R, unsigned int L, unsigned int D, unsigned int NT, unsigned int K>
double bench_compressed(size_t dim, size_t trials)
{
	using Field = CompressedField<P,T,N,R,L,D>;
	Field field{};
	CompressedMatrixDomain<Field> CMD(field);
	using Matrix = typename CompressedMatrixDomain<Field>::Matrix;
	Matrix A(field,dim,dim);
	Matrix B(field,dim,dim);
	Matrix C(field,dim,dim);
	A.random();
	B.random();

	for (volatile size_t i = 0; i < 1; ++i) {
		CMD.mul_compressed(C, B, A);
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (volatile size_t i = 0; i < trials; ++i) {
		CMD.mul_compressed(C, B, A);
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	auto mics1 = duration_cast<nanoseconds>( t2 - t1 ).count();
	double ms1 = double(mics1);
	double ops = 2.0f*dim*dim*dim*trials;
	double gffops = ops / ms1;
	
	return gffops;
}

template <unsigned int P, typename T, unsigned int N, unsigned int R, unsigned int L, unsigned int D, unsigned int NT, unsigned int K>
double bench_four_russians(size_t dim, size_t trials)
{
	using Field = CompressedField<P,T,N,R,L,D>;
	Field field{};
	CompressedMatrixDomain<Field> CMD(field);
	using Matrix = typename CompressedMatrixDomain<Field>::Matrix;
	Matrix A(field,dim,dim);
	Matrix B(field,dim,dim);
	Matrix C(field,dim,dim);
	A.random();
	B.random();
	
	for (volatile size_t i = 0; i < 1; ++i) {
		CMD. template mul_four_russians<NT,K>(C,B,A);
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (volatile size_t i = 0; i < trials; ++i) {
		CMD. template mul_four_russians<NT,K>(C,B,A);
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	auto mics1 = duration_cast<nanoseconds>( t2 - t1 ).count();
	double ms1 = double(mics1);
	double ops = 2.0f*dim*dim*dim*trials;
	double gffops = ops / ms1;
	
	return gffops;
}

template <unsigned int P, typename T, unsigned int N, unsigned int R, unsigned int L, unsigned int D>
double bench_four_russians2(size_t dim, size_t trials)
{
	using Field = CompressedField<P,T,N,R,L,D>;
	Field field{};
	CompressedMatrixDomain<Field> CMD(field);
	using Matrix = typename CompressedMatrixDomain<Field>::Matrix;
	Matrix A(field,dim,dim);
	Matrix B(field,dim,dim);
	Matrix C(field,dim,dim);
	A.random();
	B.random();
	
	//std::cout << "A: A.rowdim() = " << A.rowdim() << ", A.coldim() = " << A.coldim() << ", A.words() = " << A.getwords() << ", A.loff() = " << A.getloff() << ", A.roff() = " << A.getroff() << ", A.stride() = " << A.getstride() << ", A.alloc = " << A.alloc << std::endl;
	//std::cout << "B: B.rowdim() = " << B.rowdim() << ", B.coldim() = " << B.coldim() << ", B.words() = " << B.getwords() << ", B.loff() = " << B.getloff() << ", B.roff() = " << B.getroff() << ", B.stride() = " << B.getstride() << ", B.alloc = " << B.alloc << std::endl;
	//std::cout << "C: C.rowdim() = " << C.rowdim() << ", C.coldim() = " << C.coldim() << ", C.words() = " << C.getwords() << ", C.loff() = " << C.getloff() << ", C.roff() = " << C.getroff() << ", C.stride() = " << C.getstride() << ", C.alloc = " << C.alloc << std::endl;

	
	volatile double bestg = 0.0f;
	volatile int bestk = 2;
	
	for (volatile unsigned int K = 2; K <= 16; ++K) {
		for (volatile size_t i = 0; i < 1; ++i) {
			CMD.mul_four_russians(C,B,A,K);
		}

		high_resolution_clock::time_point t1 = high_resolution_clock::now();
		for (volatile size_t i = 0; i < trials; ++i) {
			CMD.mul_four_russians(C,B,A,K);
		}
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
		auto mics1 = duration_cast<nanoseconds>( t2 - t1 ).count();
		double ms1 = double(mics1);
		double ops = 2.0f*dim*dim*dim*trials;
		double gffops = ops / ms1;
		std::cout << "(" << K << " " << gffops << ") ";
		if (gffops > bestg) {
			bestg = gffops;
			bestk = K;
		}
	}
	std::cout << bestk << " " << bestg;
	return bestg;
}

template <unsigned int P, typename T, unsigned int N, unsigned int R, unsigned int L, unsigned int D, unsigned int NT, unsigned int K>
void test(size_t dim, size_t trials)
{
	std::cout << P << " " << 8*sizeof(T) << " " << N << " " << R << " " << L << " " << D << " " << NT << " " << K << std::endl;
	std::cout << "dim = " << dim << ", trials = " << trials << ", words = " << CompressedMatrixDomain<CompressedField<P,T,N,R,L,D>>::Matrix::get_words(0,0,dim) << ", epb = " << CompressedMatrixDomain<CompressedField<P,T,N,R,L,D>>::Word::entries_per_base() << std::endl;
	//std::cout << "\tclassical   : " << bench_classical<P,T,N,R,L,D,NT,K>(dim,trials) << std::endl;
	std::cout << "\tcompressed  : " << bench_compressed<P,T,N,R,L,D,NT,K>(dim,trials) << std::endl;
	//std::cout << "\tfourrusians : " << bench_four_russians<P,T,N,R,L,D,NT,K>(dim,trials) << std::endl;
	std::cout << "\tfourrusians2: "; bench_four_russians2<P,T,N,R,L,D>(dim,trials); std::cout << std::endl;
	std::cout << std::endl;
}

template <unsigned int P, typename T, unsigned int N, unsigned int R, unsigned int L, unsigned int D, unsigned int NT, unsigned int K>
double bench_normalize(size_t dim, size_t trials)
{
	using Field = CompressedField<P,T,N,R,L,D>;
	Field field{};
	CompressedMatrixDomain<Field> CMD(field);
	using Matrix = typename CompressedMatrixDomain<Field>::Matrix;
	Matrix A(field,dim,dim);
	A.random();
	
	//std::cout << A.rowdim() << " " << A.getwords() << std::endl;
	
	for (volatile size_t i = 0; i < 1; ++i) {
		CMD.normalize(A);
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (volatile size_t i = 0; i < trials; ++i) {
		CMD.normalize(A);
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	auto mics1 = duration_cast<nanoseconds>( t2 - t1 ).count();
	double ms1 = double(mics1);
	double ops = 2.0f*dim*dim*dim*trials;
	double gffops = ops / ms1;
	
	std::cout << ms1 / (1.0*trials) << " ns" << std::endl;
	return gffops;
}

template <unsigned int P, typename T, unsigned int N, unsigned int R, unsigned int L, unsigned int D, unsigned int NT, unsigned int K>
double bench_addin(size_t dim, size_t trials)
{
	using Field = CompressedField<P,T,N,R,L,D>;
	Field field{};
	CompressedMatrixDomain<Field> CMD(field);
	using Matrix = typename CompressedMatrixDomain<Field>::Matrix;
	Matrix A(field,dim,dim);
	A.random();
	Matrix B(field,dim,dim);
	B.random();
	
	//std::cout << A.rowdim() << " " << A.getwords() << std::endl;
	
	for (volatile size_t i = 0; i < 1; ++i) {
		CMD.addin(A,B);
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (volatile size_t i = 0; i < trials; ++i) {
		CMD.addin(A,B);
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	auto mics1 = duration_cast<nanoseconds>( t2 - t1 ).count();
	double ms1 = double(mics1);
	double ops = 2.0f*dim*dim*dim*trials;
	double gffops = ops / ms1;
	
	std::cout << ms1 / (1.0*trials) << " ns" << std::endl;
	return gffops;
}

constexpr unsigned int logo2(unsigned int v)
{
	unsigned int r = 0;
	while (v >>= 1) {
		++r;
	}
	return 1+r;
}

void test_normalize(size_t dim, size_t trials)
{
	#define MP 3u
	using MT = uint32_t;
	//#define MN 1u
	//#define MR 5u
	#define ML 1u
	
	//test<MP,MT,MN,logo2(MP),ML,1u,1u,8u>(dim,trials);
	bench_normalize<MP,uint64_t,1,logo2(MP),ML,1u,1u,8u>(dim, trials);
	bench_normalize<MP,uint64_t,4,logo2(MP),ML,1u,1u,8u>(dim, trials);
	bench_normalize<MP,uint32_t,1,logo2(MP),ML,1u,1u,8u>(dim, trials);
	bench_normalize<MP,uint32_t,8,logo2(MP),ML,1u,1u,8u>(dim, trials);
	
	bench_addin<MP,uint64_t,1,logo2(MP),ML,1u,1u,8u>(dim, trials);
	bench_addin<MP,uint64_t,4,logo2(MP),ML,1u,1u,8u>(dim, trials);
	bench_addin<MP,uint32_t,1,logo2(MP),ML,1u,1u,8u>(dim, trials);
	bench_addin<MP,uint32_t,8,logo2(MP),ML,1u,1u,8u>(dim, trials);
}
#endif

int main(int argc, char **argv) {
	size_t dim = 256;
	size_t trials = 10;
	if (argc == 3){
	dim = (size_t)std::stol(argv[1]);
	trials = (size_t)std::stol(argv[2]);
	}
	
	
	srand(time(NULL));
	
	double gffops = bench_m4ri(dim, trials);
	std::cout << "gffops: " << gffops << std::endl;
	//test_normalize(dim, trials);
	/*
	std::cout << std::endl << "---- GF2 ----" << std::endl;	
	test<2u,uint32_t,1u,1u,0u,1u,1u,8u>(dim,trials);
	test<2u,uint64_t,1u,1u,0u,1u,1u,8u>(dim,trials);
	test<2u,uint64_t,4u,1u,0u,1u,1u,8u>(dim,trials);
	test<2u,uint32_t,8u,1u,0u,1u,1u,8u>(dim,trials);
	
	std::cout << std::endl << "---- SLICED3 ----" << std::endl;
	test<3u,uint32_t,1u,1u,0u,2u,1u,8u>(dim,trials);
	test<3u,uint64_t,1u,1u,0u,2u,1u,8u>(dim,trials);
	test<3u,uint64_t,4u,1u,0u,2u,1u,8u>(dim,trials);
	test<3u,uint32_t,8u,1u,0u,2u,1u,8u>(dim,trials);
		
	std::cout << std::endl << "---- PACKED3(2,2) ----" << std::endl;
	test<3u,uint32_t,1u,2u,2u,1u,1u,8u>(dim,trials);
	test<3u,uint64_t,1u,2u,2u,1u,1u,8u>(dim,trials);
	test<3u,uint64_t,4u,2u,2u,1u,1u,8u>(dim,trials);
	test<3u,uint32_t,8u,2u,2u,1u,1u,8u>(dim,trials);
		
	std::cout << std::endl << "---- PACKED3(2,6) ----" << std::endl;
	test<3u,uint32_t,1u,2u,6u,1u,1u,8u>(dim,trials);
	test<3u,uint64_t,1u,2u,6u,1u,1u,8u>(dim,trials);
	test<3u,uint64_t,4u,2u,6u,1u,1u,8u>(dim,trials);
	test<3u,uint32_t,8u,2u,6u,1u,1u,8u>(dim,trials);
	
	std::cout << std::endl << "---- PACKED5(3,1) ----" << std::endl;
	test<5u,uint32_t,1u,3u,1u,1u,1u,8u>(dim,trials);
	test<5u,uint64_t,1u,3u,1u,1u,1u,8u>(dim,trials);
	test<5u,uint64_t,4u,3u,1u,1u,1u,8u>(dim,trials);
	test<5u,uint32_t,8u,3u,1u,1u,1u,8u>(dim,trials);
	
	std::cout << std::endl << "---- PACKED5(3,5) ----" << std::endl;
	test<5u,uint32_t,1u,3u,5u,1u,1u,8u>(dim,trials);
	test<5u,uint64_t,1u,3u,5u,1u,1u,8u>(dim,trials);
	test<5u,uint64_t,4u,3u,5u,1u,1u,8u>(dim,trials);
	test<5u,uint32_t,8u,3u,5u,1u,1u,8u>(dim,trials);
	
	
	std::cout << std::endl << "---- PACKED17(5,3) ----" << std::endl;
	test<17u,uint32_t,1u,5u,3u,1u,1u,8u>(dim,trials);
	test<17u,uint64_t,1u,5u,3u,1u,1u,8u>(dim,trials);
	test<17u,uint64_t,4u,5u,3u,1u,1u,8u>(dim,trials);
	test<17u,uint32_t,8u,5u,3u,1u,1u,8u>(dim,trials);
	
	std::cout << std::endl << "---- PACKED127(7,1) ----" << std::endl;
	test<127u,uint32_t,1u,7u,1u,1u,1u,8u>(dim,trials);
	test<127u,uint64_t,1u,7u,1u,1u,1u,8u>(dim,trials);
	test<127u,uint64_t,4u,7u,1u,1u,1u,8u>(dim,trials);
	test<127u,uint32_t,8u,7u,1u,1u,1u,8u>(dim,trials);

	std::cout << std::endl << "---- PACKED33151(15,1) ----" << std::endl;
	test<33151u,uint32_t,1u,15u,1u,1u,1u,2u>(dim,trials);
	test<33151u,uint64_t,1u,15u,1u,1u,1u,2u>(dim,trials);
	test<33151u,uint64_t,4u,15u,1u,1u,1u,2u>(dim,trials);
	test<33151u,uint32_t,8u,15u,1u,1u,1u,2u>(dim,trials);
	*/
	return 0;
}

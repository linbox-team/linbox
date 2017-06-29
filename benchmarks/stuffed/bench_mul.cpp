#include <iostream>
#include <chrono>
#include <limits.h>
#include <array>
#include <stdlib.h>
#include <string.h>
#include <bitset>
#include <time.h>

using namespace std::chrono;

#include "compressed-matrix-domain.h"

template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
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

template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D>
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
	
	//std::cout << "A: A.rowdim() = " << A.rowdim() << ", A.coldim() = " << A.coldim() << ", A.words() = " << A.getwords() << ", A.loff() = " << A.getloff() << ", A.roff() = " << A.getroff() << ", A.stride() = " << A.getstride() << ", A.alloc = " << A.alloc << std::endl;
	//std::cout << "B: B.rowdim() = " << B.rowdim() << ", B.coldim() = " << B.coldim() << ", B.words() = " << B.getwords() << ", B.loff() = " << B.getloff() << ", B.roff() = " << B.getroff() << ", B.stride() = " << B.getstride() << ", B.alloc = " << B.alloc << std::endl;
	//std::cout << "C: C.rowdim() = " << C.rowdim() << ", C.coldim() = " << C.coldim() << ", C.words() = " << C.getwords() << ", C.loff() = " << C.getloff() << ", C.roff() = " << C.getroff() << ", C.stride() = " << C.getstride() << ", C.alloc = " << C.alloc << std::endl;

	
	volatile double bestg = 0.0f;
	volatile int bestk = 2;
	
	for (volatile uint64_t K = 6; K <= 14; ++K) {
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

template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D, uint64_t NT, uint64_t K>
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

template <uint64_t P, typename T, uint64_t N, uint64_t R, uint64_t L, uint64_t D, uint64_t NT, uint64_t K>
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

constexpr uint64_t logo2(uint64_t v)
{
	uint64_t r = 0;
	if (v > 2) {
	while (v >>= 1) {
		++r;
	}
	}
	return r+1;
}

constexpr uint64_t minbits(uint64_t p)
{
	return logo2((p-1)*(p-1));
}

constexpr uint64_t minm4rmbits(uint64_t p)
{
	//uint64_t x = 8*((1u << 8) - 1)*(p-1);
	//return logo2(x);
	return logo2(((1u << logo2(p)) - 1) * (p-1));
}

template <uint64_t P, typename T, uint64_t N, uint64_t L, bool S>
void test(size_t dim, size_t trials)
{
	if (S == true) { // use sliced
		constexpr uint64_t D = logo2(P);
		constexpr auto R = static_cast<uint64_t>(1);
		//std::cout << P << " " << 8*sizeof(T) << " " << N << " " << R << " " << 0 << " " << D << std::endl;
		std::cout << "\tcompressed  : " << bench_compressed<P,T,N,R,0,D>(dim,trials) << std::endl;
		std::cout << "\tfourrusians : "; bench_four_russians<P,T,N,R,0,D>(dim,trials); std::cout << std::endl;
		std::cout << std::endl;
	}
	else { // packed
		constexpr auto R = logo2(P);
		constexpr auto NLC = minbits(P) + L - R;
		constexpr auto NLR = minm4rmbits(P) + L - R;
		//std::cout << P << " " << 8*sizeof(T) << " " << N << " " << R << " " << L << " " << D << std::endl;
		std::cout << "\tcompressed  (" << R << "," << NLC << "): " << bench_compressed<P,T,N,R,NLC,1>(dim,trials) << std::endl;
		std::cout << "\tfourrusians (" << R << "," << NLR << "): "; bench_four_russians<P,T,N,R,NLR,1>(dim,trials); std::cout << std::endl;
		std::cout << std::endl;
	}
}

template <uint64_t P>
void bench_packed(size_t dim, size_t trials) {
/*
	std::cout << P << " : L = 0" << std::endl;
	test<P,uint32_t,1u,0u,false>(dim,trials);
	test<P,uint64_t,1u,0u,false>(dim,trials);
	test<P,uint64_t,4u,0u,false>(dim,trials);
	test<P,uint32_t,8u,0u,false>(dim,trials);
	
	std::cout << P << " : L = 1" << std::endl;
	test<P,uint32_t,1u,1u,false>(dim,trials);
	test<P,uint64_t,1u,1u,false>(dim,trials);
	test<P,uint64_t,4u,1u,false>(dim,trials);
	test<P,uint32_t,8u,1u,false>(dim,trials);
	
	std::cout << P << " : L = 2" << std::endl;
	test<P,uint32_t,1u,2u,false>(dim,trials);
	test<P,uint64_t,1u,2u,false>(dim,trials);
	test<P,uint64_t,4u,2u,false>(dim,trials);
	test<P,uint32_t,8u,2u,false>(dim,trials);
	
	std::cout << P << " : L = 3" << std::endl;
	test<P,uint32_t,1u,3u,false>(dim,trials);
	test<P,uint64_t,1u,3u,false>(dim,trials);
	test<P,uint64_t,4u,3u,false>(dim,trials);
	test<P,uint32_t,8u,3u,false>(dim,trials);
	
	std::cout << P << " : L = 4" << std::endl;
	test<P,uint32_t,1u,4u,false>(dim,trials);
	test<P,uint64_t,1u,4u,false>(dim,trials);
	test<P,uint64_t,4u,4u,false>(dim,trials);
	test<P,uint32_t,8u,4u,false>(dim,trials);
	
	std::cout << P << " : L = 5" << std::endl;
	test<P,uint32_t,1u,5u,false>(dim,trials);
	test<P,uint64_t,1u,5u,false>(dim,trials);
	test<P,uint64_t,4u,5u,false>(dim,trials);
	test<P,uint32_t,8u,5u,false>(dim,trials);
	
	std::cout << P << " : L = 6" << std::endl;
	test<P,uint32_t,1u,6u,false>(dim,trials);
	test<P,uint64_t,1u,6u,false>(dim,trials);
	test<P,uint64_t,4u,6u,false>(dim,trials);
	test<P,uint32_t,8u,6u,false>(dim,trials);
	
	std::cout << P << " : L = 7" << std::endl;
	test<P,uint32_t,1u,7u,false>(dim,trials);
	test<P,uint64_t,1u,7u,false>(dim,trials);
	test<P,uint64_t,4u,7u,false>(dim,trials);
	test<P,uint32_t,8u,7u,false>(dim,trials);
*/
	test<P,uint32_t,8u,7u,false>(dim,trials);
}


int main(int argc, char **argv) {
	size_t dim = (size_t)std::stol(argv[1]);
	size_t trials = (size_t)std::stol(argv[2]);
	
	
	srand(time(NULL));
	bench_packed<17>(dim,trials);
	/*
	std::cout << std::endl << "---- GF2 ----" << std::endl;	
	test<2u,uint32_t,1u,0u,true>(dim,trials);
	test<2u,uint64_t,1u,0u,true>(dim,trials);
	test<2u,uint64_t,4u,0u,true>(dim,trials);
	test<2u,uint32_t,8u,0u,true>(dim,trials);
	
	std::cout << std::endl << "---- SLICED3 ----" << std::endl;
	test<3u,uint32_t,1u,0u,true>(dim,trials);
	test<3u,uint64_t,1u,0u,true>(dim,trials);
	test<3u,uint64_t,4u,0u,true>(dim,trials);
	test<3u,uint32_t,8u,0u,true>(dim,trials);
	
	std::cout << std::endl << "---- PACKED3(2,1) ----" << std::endl;
	test<3u,uint32_t,1u,1u,false>(dim,trials);
	test<3u,uint64_t,1u,1u,false>(dim,trials);
	test<3u,uint64_t,4u,1u,false>(dim,trials);
	test<3u,uint32_t,8u,1u,false>(dim,trials);	
	*/
	return 0;
}

#include <iostream>
#include <chrono>
#include <limits.h>
#include <array>
#include <stdlib.h>
#include <string.h>
#include <bitset>
#include <time.h>
#include <vector>

using namespace std::chrono;

#include "compressed-matrix-domain.h"

#define VEC_SIZE 1000000
#define AXPY_COUNT 1000
#define TRIALS 10

constexpr uint64_t logo2(uint64_t v)
{
	uint64_t r = 0;
	if (v == 2) r = 1;
	if (v == 6) r = 2;
	if (v == 4) r = 2;
	if (v == 20) r = 4;
	if (v == 6) r = 2;
	if (v == 42) r = 5;
	if (v == 16) r = 4;
	if (v == 17*16) r = 8;
	if (v == 126) r = 6;
	if (v == 126*127) r = 13;
	return r+1;
}

constexpr uint64_t minbits(uint64_t p)
{
	return logo2((p)*(p-1)) - logo2(p-1);
}

template <unsigned int P, typename T>
constexpr uint64_t maxbits()
{
	return 8*sizeof(T) - minbits(P) - logo2(P-1);
}

template <unsigned int P, typename T, unsigned int N, unsigned int L>
double trial(unsigned int seed)
{
	
	using Field = CompressedField<P,T,N,logo2(P-1),minbits(P)+L,1>;
	using CMD = CompressedMatrixDomain<Field>;
	using Matrix = typename CMD::Matrix;
	
	Field f{};
	CMD cmd(f);
	
	Matrix A(f,1,VEC_SIZE);
	Matrix B(f,1,VEC_SIZE);
	
	if (seed == 0) {
		srand(time(NULL));
	}
	else {
		srand(seed);
	}
	T coeffs[AXPY_COUNT];
	for (unsigned int i = 0; i < AXPY_COUNT; ++i) {
		coeffs[i] = rand() % P;
	}
	
	A.random();
	B.random();
	
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (volatile int t = 0; t < TRIALS; ++t) {
		for (unsigned int i = 0; i < AXPY_COUNT; ++i) {
			cmd.axpyin(B,A,coeffs[i]);
			if ((i % f.axpy_freq() == 0 && i != 0) or (i == AXPY_COUNT-1)) {
				cmd.normalize(B);
			}
		}
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	auto mics1 = duration_cast<nanoseconds>( t2 - t1 ).count();
	double ops = mics1 / (1.0f * TRIALS);	
	std::cout << Field::Unit::get_l() << " ";// << ops;
	return ops;
}

template <unsigned int N>
void test3_32() {
	std::cout << "P = 3, 32 by " << N << std::endl;
	std::cout << trial<3,uint32_t,N,0>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,1>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,2>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,3>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,4>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,5>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,6>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,7>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,8>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,9>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,10>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,11>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,12>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,13>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,14>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,15>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,16>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,17>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,18>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,19>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,20>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,21>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,22>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,23>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,24>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,25>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,26>(0) << std::endl;
	std::cout << trial<3,uint32_t,N,27>(0) << std::endl;
}

template <unsigned int N>
void test3_64() {
	std::cout << "P = 3, 64 by " << N << std::endl;
	std::cout << trial<3,uint64_t,N,0>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,1>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,2>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,3>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,4>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,5>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,6>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,7>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,8>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,9>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,10>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,11>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,12>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,13>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,14>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,15>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,16>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,17>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,18>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,19>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,20>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,21>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,22>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,23>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,24>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,25>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,26>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,27>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,28>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,29>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,30>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,31>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,32>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,33>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,34>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,35>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,36>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,37>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,38>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,39>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,40>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,41>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,42>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,43>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,44>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,45>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,46>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,47>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,48>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,49>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,50>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,51>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,52>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,53>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,54>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,55>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,56>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,57>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,58>(0) << std::endl;
	std::cout << trial<3,uint64_t,N,59>(0) << std::endl;
}

template <unsigned int N>
void test5_32() {
	std::cout << "P = 5, 32 by " << N << std::endl;
	std::cout << trial<5,uint32_t,N,0>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,1>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,2>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,3>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,4>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,5>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,6>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,7>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,8>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,9>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,10>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,11>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,12>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,13>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,14>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,15>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,16>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,17>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,18>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,19>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,20>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,21>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,22>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,23>(0) << std::endl;
	std::cout << trial<5,uint32_t,N,24>(0) << std::endl;
}

template <unsigned int N>
void test5_64() {
	std::cout << "P = 5, 64 by " << N << std::endl;
	std::cout << trial<5,uint64_t,N,0>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,1>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,2>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,3>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,4>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,5>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,6>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,7>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,8>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,9>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,10>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,11>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,12>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,13>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,14>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,15>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,16>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,17>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,18>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,19>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,20>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,21>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,22>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,23>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,24>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,25>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,26>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,27>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,28>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,29>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,30>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,31>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,32>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,33>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,34>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,35>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,36>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,37>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,38>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,39>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,40>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,41>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,42>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,43>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,44>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,45>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,46>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,47>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,48>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,49>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,50>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,51>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,52>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,53>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,54>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,55>(0) << std::endl;
	std::cout << trial<5,uint64_t,N,56>(0) << std::endl;
}

template <unsigned int N>
void test7_32() {
	std::cout << "P = 7, 32 by " << N << std::endl;
	std::cout << trial<7,uint32_t,N,0>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,1>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,2>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,3>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,4>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,5>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,6>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,7>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,8>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,9>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,10>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,11>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,12>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,13>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,14>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,15>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,16>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,17>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,18>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,19>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,20>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,21>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,22>(0) << std::endl;
	std::cout << trial<7,uint32_t,N,23>(0) << std::endl;
}

template <unsigned int N>
void test7_64() {
	std::cout << "P = 7, 64 by " << N << std::endl;
	std::cout << trial<7,uint64_t,N,0>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,1>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,2>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,3>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,4>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,5>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,6>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,7>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,8>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,9>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,10>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,11>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,12>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,13>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,14>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,15>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,16>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,17>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,18>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,19>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,20>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,21>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,22>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,23>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,24>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,25>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,26>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,27>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,28>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,29>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,30>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,31>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,32>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,33>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,34>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,35>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,36>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,37>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,38>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,39>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,40>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,41>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,42>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,43>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,44>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,45>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,46>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,47>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,48>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,49>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,50>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,51>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,52>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,53>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,54>(0) << std::endl;
	std::cout << trial<7,uint64_t,N,55>(0) << std::endl;
}

template <unsigned int N>
void test17_32() {
	std::cout << "P = 17, 32 by " << N << std::endl;
	std::cout << trial<17,uint32_t,N,0>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,1>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,2>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,3>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,4>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,5>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,6>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,7>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,8>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,9>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,10>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,11>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,12>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,13>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,14>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,15>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,16>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,17>(0) << std::endl;
	std::cout << trial<17,uint32_t,N,18>(0) << std::endl;
}

template <unsigned int N>
void test17_64() {
	std::cout << "P = 17, 64 by " << N << std::endl;
	std::cout << trial<17,uint64_t,N,0>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,1>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,2>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,3>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,4>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,5>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,6>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,7>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,8>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,9>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,10>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,11>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,12>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,13>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,14>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,15>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,16>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,17>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,18>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,19>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,20>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,21>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,22>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,23>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,24>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,25>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,26>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,27>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,28>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,29>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,30>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,31>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,32>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,33>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,34>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,35>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,36>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,37>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,38>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,39>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,40>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,41>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,42>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,43>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,44>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,45>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,46>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,47>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,48>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,49>(0) << std::endl;
	std::cout << trial<17,uint64_t,N,50>(0) << std::endl;
}

template <unsigned int N>
void test127_32() {
	std::cout << "P = 127, 32 by " << N << std::endl;
	std::cout << trial<127,uint32_t,N,0>(0) << std::endl;
	std::cout << trial<127,uint32_t,N,1>(0) << std::endl;
	std::cout << trial<127,uint32_t,N,2>(0) << std::endl;
	std::cout << trial<127,uint32_t,N,3>(0) << std::endl;
	std::cout << trial<127,uint32_t,N,4>(0) << std::endl;
	std::cout << trial<127,uint32_t,N,5>(0) << std::endl;
	std::cout << trial<127,uint32_t,N,6>(0) << std::endl;
	std::cout << trial<127,uint32_t,N,7>(0) << std::endl;
	std::cout << trial<127,uint32_t,N,8>(0) << std::endl;
	std::cout << trial<127,uint32_t,N,9>(0) << std::endl;
	std::cout << trial<127,uint32_t,N,10>(0) << std::endl;
	std::cout << trial<127,uint32_t,N,16>(0) << std::endl;
}

template <unsigned int N>
void test127_64() {
	std::cout << "P = 127, 64 by " << N << std::endl;
	std::cout << trial<127,uint64_t,N,0>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,1>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,2>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,3>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,4>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,5>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,6>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,7>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,8>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,9>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,10>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,11>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,12>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,13>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,14>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,15>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,16>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,17>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,18>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,19>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,20>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,21>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,22>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,23>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,24>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,25>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,26>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,27>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,28>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,29>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,30>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,31>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,32>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,33>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,34>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,35>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,36>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,37>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,38>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,39>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,40>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,41>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,42>(0) << std::endl;
	std::cout << trial<127,uint64_t,N,43>(0) << std::endl;
}

int main() {
	//std::cout << maxbits<127,uint64_t>() << std::endl;
	test3_32<1>();
  	test3_32<8>();
	test3_64<1>();
	test3_64<4>();
	
	test5_32<1>();
	test5_32<8>();
	test5_64<1>();
	test5_64<4>();
	
	test7_32<1>();
	test7_32<8>();
	test7_64<1>();
	test7_64<4>();
	
	test17_32<1>();
	test17_32<8>();
	test17_64<1>();
	test17_64<4>();
	
	test127_32<1>();
	test127_32<8>();
	test127_64<1>();
	test127_64<4>();
	return 0;
}

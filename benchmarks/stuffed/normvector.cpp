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

//#define VEC_SIZE 256
#define AXPY_COUNT 1000000
//#define TRIALS 1000000

constexpr uint64_t logo2(uint64_t v)
{
	uint64_t r = 0;
//	if (v > 2) {
		while (v >>= 1) {
			++r;
		}
//	}
	return r+1;
}

constexpr uint64_t minbits(uint64_t p)
{
	return logo2((p)*(p-1)) - logo2(p-1);
}

template <unsigned int P, typename T>
constexpr uint64_t maxbits()
{
	return 8*sizeof(T) - minbits(P) - logo2(P);
}

template <unsigned int P, typename T, unsigned int N, unsigned int L>
double trial(uint64_t VEC_SIZE, uint64_t TRIALS, unsigned int seed)
{
	
	using Field = CompressedField<P,T,N,logo2(P),minbits(P),1>;
	using CMD = CompressedMatrixDomain<Field>;
	using Matrix = typename CMD::Matrix;
	
	Field f{};
	CMD cmd(f);
	
	Matrix A(f,VEC_SIZE,VEC_SIZE);
//	Matrix B(f,VEC_SIZE,VEC_SIZE);
	
	if (seed == 0) {
		srand(time(NULL));
	}
	else {
		srand(seed);
	}
	A.random();
//	B.random();
	
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (volatile int t = 0; t < TRIALS; ++t) {
		cmd.normalize(A);
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	auto mics1 = duration_cast<nanoseconds>( t2 - t1 ).count();
	double ops = mics1 / (1.0f * TRIALS);	
	//std::cout << minbits(P)+L << " ";// << ops;
	return ops;
}

template <unsigned int N>
void test3_32(uint64_t VEC_SIZE, uint64_t TRIALS) {
	std::cout << "P = 3, 32 by " << N << std::endl;
	std::cout << trial<3,uint32_t,N,0>(VEC_SIZE,TRIALS,0) << std::endl;
}

template <unsigned int N>
void test3_64(uint64_t VEC_SIZE, uint64_t TRIALS) {
	std::cout << "P = 3, 64 by " << N << std::endl;
	std::cout << trial<3,uint64_t,N,0>(VEC_SIZE,TRIALS,0) << std::endl;
}

template <unsigned int N>
void test5_32(uint64_t VEC_SIZE, uint64_t TRIALS) {
	std::cout << "P = 5, 32 by " << N << std::endl;
	std::cout << trial<5,uint32_t,N,0>(VEC_SIZE, TRIALS,0) << std::endl;
}

template <unsigned int N>
void test5_64(uint64_t VEC_SIZE, uint64_t TRIALS) {
	std::cout << "P = 5, 64 by " << N << std::endl;
	std::cout << trial<5,uint64_t,N,0>(VEC_SIZE, TRIALS,0) << std::endl;
}

template <unsigned int N>
void test7_32(uint64_t VEC_SIZE, uint64_t TRIALS) {
	std::cout << "P = 7, 32 by " << N << std::endl;
	std::cout << trial<7,uint32_t,N,0>(VEC_SIZE, TRIALS,0) << std::endl;
}

template <unsigned int N>
void test7_64(uint64_t VEC_SIZE, uint64_t TRIALS) {
	std::cout << "P = 7, 64 by " << N << std::endl;
	std::cout << trial<7,uint64_t,N,0>(VEC_SIZE, TRIALS,0) << std::endl;
}

template <unsigned int N>
void test17_32(uint64_t VEC_SIZE, uint64_t TRIALS) {
	std::cout << "P = 17, 32 by " << N << std::endl;
	std::cout << trial<17,uint32_t,N,0>(VEC_SIZE, TRIALS,0) << std::endl;
}

template <unsigned int N>
void test17_64(uint64_t VEC_SIZE, uint64_t TRIALS) {
	std::cout << "P = 17, 64 by " << N << std::endl;
	std::cout << trial<17,uint64_t,N,0>(VEC_SIZE, TRIALS,0) << std::endl;
}

template <unsigned int N>
void test127_32(uint64_t VEC_SIZE, uint64_t TRIALS) {
	std::cout << "P = 127, 32 by " << N << std::endl;
	std::cout << trial<127,uint32_t,N,0>(VEC_SIZE, TRIALS,0) << std::endl;
}

template <unsigned int N>
void test127_64(uint64_t VEC_SIZE, uint64_t TRIALS) {
	std::cout << "P = 127, 64 by " << N << std::endl;
	std::cout << trial<127,uint64_t,N,0>(VEC_SIZE, TRIALS,0) << std::endl;
}

int main(int argc, char **argv) {
	uint64_t d = (uint64_t)std::stol(argv[1]);
	uint64_t t = (uint64_t)std::stol(argv[2]);
	//std::cout << maxbits<127,uint64_t>() << std::endl;
	test3_32<1>(d,t);
	test3_32<8>(d,t);
	test3_64<1>(d,t);
	test3_64<4>(d,t);
	
	test5_32<1>(d,t);
	test5_32<8>(d,t);
	test5_64<1>(d,t);
	test5_64<4>(d,t);
		
	test7_32<1>(d,t);
	test7_32<8>(d,t);
	test7_64<1>(d,t);
	test7_64<4>(d,t);
	
	test17_32<1>(d,t);
	test17_32<8>(d,t);
	test17_64<1>(d,t);
	test17_64<4>(d,t);
	
	test127_32<1>(d,t);
	test127_32<8>(d,t);
	test127_64<1>(d,t);
	test127_64<4>(d,t);
	
	return 0;
}

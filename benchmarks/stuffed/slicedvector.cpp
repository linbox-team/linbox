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

#define VEC_SIZE 1000
#define AXPY_COUNT 1000000
#define TRIALS 100

constexpr uint64_t logo2(uint64_t v)
{
	uint64_t r = 0;
	while (v >>= 1) {
		++r;
	}
	return r+1;
}

constexpr uint64_t minbits(uint64_t p)
{
	return logo2((p)*(p-1)) - logo2(p);
}

template <unsigned int P, typename T>
constexpr uint64_t maxbits()
{
	return 8*sizeof(T) - minbits(P) - logo2(P);
}

template <typename T, unsigned int N>
double trial(uint64_t dim, uint64_t trials, unsigned int seed)
{
	
	//using Field = CompressedField<P,T,N,logo2(P),minbits(P)+L,1>;
	using Field = CompressedField<3,T,N,1,0,2>;
	using CMD = CompressedMatrixDomain<Field>;
	using Matrix = typename CMD::Matrix;
	
	Field f{};
	CMD cmd(f);
	
	Matrix A(f,dim,dim);
	Matrix B(f,dim,dim);
	
	if (seed == 0) {
		srand(time(NULL));
	}
	else {
		srand(seed);
	}

	A.random();
	B.random();
	
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (volatile uint64_t t = 0; t < trials; ++t) {
		cmd.addin(B,A);
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	auto mics1 = duration_cast<nanoseconds>( t2 - t1 ).count();
	double ops = mics1 / (1.0f * trials);	
	//std::cout << minbits(P)+L << " ";// << ops;
	return ops;
}

int main(int argc, char **argv) {
	size_t dim = (size_t)std::stol(argv[1]);
	size_t trials = (size_t)std::stol(argv[2]);

	std::cout << trial<uint32_t,1>(dim,trials,0) << std::endl;
	std::cout << trial<uint32_t,8>(dim,trials,0) << std::endl;
	std::cout << trial<uint64_t,1>(dim,trials,0) << std::endl;
	std::cout << trial<uint64_t,4>(dim,trials,0) << std::endl;
	
	

	return 0;
}
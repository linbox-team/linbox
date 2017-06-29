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

template <unsigned int P, typename T, unsigned int N, unsigned int R, unsigned int L, unsigned int D>
void trial_four_russians(uint64_t dim, uint64_t trials, unsigned int seed, unsigned int k)
{
	
	using Field = CompressedField<P,T,N,R,L,D>;
	using CMD = CompressedMatrixDomain<Field>;
	using Matrix = typename CMD::Matrix;
	
	Field f{};
	CMD cmd(f);
	
	Matrix A(f,dim,dim);
	Matrix B(f,dim,dim);
	Matrix C(f,dim,dim);
	
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
		cmd.mul_four_russians(C,A,B,k);
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	auto mics1 = duration_cast<nanoseconds>( t2 - t1 ).count();
	double ops = (2.0f * dim * dim * dim * trials) / (double)mics1;	
	std::cout << k << " : "<< ops << std::endl;
}

void test_packed_four_russians(uint64_t dim, uint64_t trials) {
	std::cout << "2sliced" << std::endl;
	trial_four_russians<2,uint32_t,8,1,0,1>(dim,trials,0,7);
	trial_four_russians<2,uint32_t,8,1,0,1>(dim,trials,0,8);
	trial_four_russians<2,uint32_t,8,1,0,1>(dim,trials,0,9);
	trial_four_russians<2,uint32_t,8,1,0,1>(dim,trials,0,10);
	trial_four_russians<2,uint32_t,8,1,0,1>(dim,trials,0,12);
	trial_four_russians<2,uint32_t,8,1,0,1>(dim,trials,0,14);

	std::cout << "3packed" << std::endl;
	trial_four_russians<3,uint32_t,8,2,6,1>(dim,trials,0,7);
	trial_four_russians<3,uint32_t,8,2,6,1>(dim,trials,0,8);
	trial_four_russians<3,uint32_t,8,2,6,1>(dim,trials,0,9);
	trial_four_russians<3,uint32_t,8,2,6,1>(dim,trials,0,10);
	trial_four_russians<3,uint32_t,8,2,6,1>(dim,trials,0,12);
	trial_four_russians<3,uint32_t,8,2,6,1>(dim,trials,0,14);
	
	std::cout << "3sliced" << std::endl;
	trial_four_russians<3,uint32_t,8,1,0,2>(dim,trials,0,7);
	trial_four_russians<3,uint32_t,8,1,0,2>(dim,trials,0,8);
	trial_four_russians<3,uint32_t,8,1,0,2>(dim,trials,0,9);
	trial_four_russians<3,uint32_t,8,1,0,2>(dim,trials,0,10);
	trial_four_russians<3,uint32_t,8,1,0,2>(dim,trials,0,12);
	trial_four_russians<3,uint32_t,8,1,0,2>(dim,trials,0,14);

	std::cout << "5packed" << std::endl;
	trial_four_russians<5,uint32_t,8,3,5,1>(dim,trials,0,7);
	trial_four_russians<5,uint32_t,8,3,5,1>(dim,trials,0,8);
	trial_four_russians<5,uint32_t,8,3,5,1>(dim,trials,0,9);
	trial_four_russians<5,uint32_t,8,3,5,1>(dim,trials,0,10);
	trial_four_russians<5,uint32_t,8,3,5,1>(dim,trials,0,12);
	trial_four_russians<5,uint32_t,8,3,5,1>(dim,trials,0,14);
	
	std::cout << "7packed" << std::endl;
	trial_four_russians<7,uint32_t,8,3,7,1>(dim,trials,0,7);
	trial_four_russians<7,uint32_t,8,3,7,1>(dim,trials,0,8);
	trial_four_russians<7,uint32_t,8,3,7,1>(dim,trials,0,9);
	trial_four_russians<7,uint32_t,8,3,7,1>(dim,trials,0,10);
	trial_four_russians<7,uint32_t,8,3,7,1>(dim,trials,0,12);
	trial_four_russians<7,uint32_t,8,3,7,1>(dim,trials,0,14);
	
	std::cout << "17packed" << std::endl;
	trial_four_russians<17,uint32_t,8,5,11,1>(dim,trials,0,7);
	trial_four_russians<17,uint32_t,8,5,11,1>(dim,trials,0,8);
	trial_four_russians<17,uint32_t,8,5,11,1>(dim,trials,0,9);
	trial_four_russians<17,uint32_t,8,5,11,1>(dim,trials,0,10);
	trial_four_russians<17,uint32_t,8,5,11,1>(dim,trials,0,12);
	trial_four_russians<17,uint32_t,8,5,11,1>(dim,trials,0,14);
	
	std::cout << "127packed" << std::endl;
	trial_four_russians<127,uint32_t,8,7,9,1>(dim,trials,0,7);
	trial_four_russians<127,uint32_t,8,7,9,1>(dim,trials,0,8);
	trial_four_russians<127,uint32_t,8,7,9,1>(dim,trials,0,9);
	trial_four_russians<127,uint32_t,8,7,9,1>(dim,trials,0,10);
	trial_four_russians<127,uint32_t,8,7,9,1>(dim,trials,0,12);
	trial_four_russians<127,uint32_t,8,7,9,1>(dim,trials,0,14);

}

int main(int argc, char **argv) {
	size_t dim = (size_t)std::stol(argv[1]);
	size_t trials = (size_t)std::stol(argv[2]);

	std::cout << std::endl << "FOUR RUSSIANS" << std::endl;
	test_packed_four_russians(dim,trials);
	

	return 0;
}

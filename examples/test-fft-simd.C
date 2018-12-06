#include <memory>
#include <vector>
#include <iostream>

#include <linbox/linbox-config.h>

#include <givaro/givranditer.h>

#include "fflas-ffpack/fflas-ffpack.h"

#include "linbox/algorithms/polynomial-matrix/polynomial-fft-butterflies.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-fft-algorithms.h"

//#include "linbox/algorithms/polynomial-matrix/polynomial-fft-transform.h"
#include "linbox/randiter/random-fftprime.h"
#include "linbox/ring/modular.h"

template<class T>
using Allocator = AlignedAllocator<T, Alignment::CACHE_LINE>;

template<class Simd_t, class Field, class Vec>
Vec test_fft(const Field & F, const Vec & a, uint64_t l2n, typename Field::Element& w = 0){
	Vec c(a.size(), 0), a2(a), amodp(a);


	LinBox::FFT_init<Field> fft_init_evaluation(F, l2n, w);
	w = fft_init_evaluation.getRoot();

	LinBox::FFT_init<Field> fft_init_interpolation(F, l2n, fft_init_evaluation.getInvRoot());

	// Init fft algos
	LinBox::FFT_algorithms<Field, Simd_t> fft_algo_evaluation(fft_init_evaluation);
	LinBox::FFT_algorithms<Field, Simd_t> fft_algo_interpolation(fft_init_interpolation);

		for (size_t i = 0; i < a2.size(); i++)
			F.init(amodp[i], a2[i]);
		std::cout << "Original a :" << std::endl
				  << amodp << std::endl;

	// Evaluate vector
	fft_algo_evaluation.DIF(a2.data());

	for (size_t i = 0; i < a2.size(); i++)
		F.init(amodp[i], a2[i]);
	std::cout << "Evaluated a :" << std::endl
			  << amodp << std::endl;

	// Interpolate the result
	fft_algo_interpolation.DIT(a2.data());

	for (size_t i = 0; i < a2.size(); i++)
		F.init(amodp[i], a2[i]);
	std::cout << "Interpolated a :" << std::endl
			  << amodp << std::endl;

	typename Field::Element invn;
	F.inv(invn, (typename Field::Element)a.size());

	// multiply by n^-1
	for(uint64_t i = 0 ; i < a.size() ; ++i){
		F.mulin(a2[i], invn);
	}

	std::cout << "Final a :" << std::endl
			  << a2 << std::endl;

	return a2;
}

int main(int argc, char** argv){

	// using Field = Givaro::Modular<uint32_t, uint64_t>; // Simd256 -> Simd128
	using Field = Givaro::Modular<double,double>;
	using Element  = typename Field::Element;

	uint64_t l2n = 3;
	uint64_t n = 1_ui64 << l2n;
	Element w = 0;
	
	// Field F(36175873); // 25bits prime generate with LinBox::RandomFFTPrime
	Field F(1048609); // 25bits prime generate with LinBox::RandomFFTPrime

	std::vector<Element, Allocator<Element>> a(n), d(n);

	typename Field::RandIter gen(F);

	for(uint64_t i = 0 ; i < n ; ++i){
		gen(a[i]);
	}

	auto c = test_fft<NoSimd<Element>>(F, a, l2n, w);

	bool ok = true;
	for(uint64_t i = 0 ; i < n ; ++i){
		if(c[i] != a[i]){
			ok = false;
		}
	}

	std::cout << "FFT NoSimd: " << ((ok)?"OK":"KO") << std::endl;
	//	std::cout << "c : " << c << std::endl << std::endl << std::endl;

	c = test_fft<Simd256<Element>>(F, a, l2n, w);

	ok = true;
	for(uint64_t i = 0 ; i < n ; ++i){
		if(c[i] != a[i]){
			ok = false;
		}
	}

	std::cout << "FFT Simd256: " << ((ok)?"OK":"KO") << std::endl;
	//	std::cout << "c : " << c << std::endl << std::endl << std::endl;

	return EXIT_SUCCESS;
}

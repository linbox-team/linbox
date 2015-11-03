/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/*
 * Copyright (C) 2014  Pascal Giorgi, Romain Lebreton
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
 *            Romain Lebreton <romain.lebreton@lirmm.fr>
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
 */


#ifndef __LINBOX_polynomial_fft_transform_H
#define __LINBOX_polynomial_fft_transform_H

#include <iostream>

#include "linbox/algorithms/polynomial-matrix/simd.h"
#include "linbox/util/debug.h"
#include "givaro/givinteger.h"

// template<typename T>
// std::ostream& operator<<(std::ostream& os, const std::vector<T> &x){
// 	std::ostream_iterator<T> out_it (os,", ");
// 	std::copy ( x.begin(), x.end(), out_it );
// 	return os;
// }

#include "fflas-ffpack/utils/align-allocator.h"

namespace LinBox {


	// class to handle FFT transform over wordsize prime field Fp (p < 2^29)
	template <class Field>
	class FFT_transform {
	public:
		typedef typename Field::Element Element;

		const Field                *fld;
		uint64_t              _pl, _dpl;
		uint64_t                      n;
		uint64_t                     ln;
		//uint64_t                  _logp;
		//uint64_t                     _I;
		//double                    _pinv;
		uint32_t                      _w;
		uint32_t                   _invw;
		typedef std::vector<uint32_t,AlignedAllocator<uint32_t, Alignment::DEFAULT> > VECT;
		VECT    pow_w;
		VECT   pow_wp; // Precomputations in shoup
		VECT    _data;
		Element                      _p;
		//   pow_w = table of roots of unity. If w = primitive K-th root, then the table is:
		//           1, w, w^2, ..., w^{K/2-1},
		//           1, w^2, w^4, ..., w^{K/2-2},
		//           1, w^4, w^8, ..., w^{K/2-4}
		//           ...
		//           1, w^{K/8}, w^{K/4}, w^{3K/8},
		//           1, w^{K/4},
		//           1.

		inline const Field & field() const { return *fld; }

		uint64_t find_gen (uint64_t _m, uint64_t _val2p) {
			// find a primitive 2^k root of unity where
			// _p - 1 = 2^k * m
			srand((unsigned int) time(NULL));
			uint64_t y,z,j;
			uint64_t _gen;
			for (;;) {
				_gen = rand() % _pl; if (_gen <= 0) continue;
				z = 1;
				for (unsigned long i=0; i < _m; ++i) z = z*_gen % _pl;
				if (z == 1) continue;
				// _gen^i =/ 1 pour 0 <= i < m
				_gen = z;
				j = 0;
				do {
					y = z;
					z = y*y % _pl;
					j++;
				} while (j != _val2p && z != 1);
				if (j == _val2p)
					break;
			}
			return _gen;
		}

		FFT_transform (const Field& fld2, size_t ln2, Element w = 0)
			: fld (&fld2), n ((1U << ln2)), ln (ln2), pow_w(n - 1), pow_wp(n - 1), _data(n) {
			_pl = fld->characteristic();
			_p  = fld->characteristic();

			linbox_check((_pl >> 29) == 0 ); // 8*p < 2^31 for Harvey's butterflies
			_dpl = (_pl << 1);
			//_pinv = 1 / (double) _pl;
							
			if (w == 0){   // find a pseudo 2^lpts-th primitive root of unity
				uint64_t _val2p = 0;
				uint64_t     _m = _pl;
				uint64_t  _logp = 0;
				while (_m) {
					_m >>= 1;
					_logp++;
				}
				_m = _pl - 1;
				while ((_m & 1) == 0) {
					_m >>= 1;
					_val2p++;
				}
				//_I = (1L << (_logp << 1)) / _pl;
				uint64_t _gen = find_gen (_m, _val2p);
				_w = Givaro::powmod(_gen, 1<<(_val2p-ln), _pl);
			}
			else {
				_w = (uint32_t)w;
			}
			
			// compute w^(-1) mod p = w^(2^lpts - 1)
			_invw = Givaro::powmod(_w, (1<<ln) - 1, _pl);
			
			size_t pos = 0;
			uint64_t wi = 1;
			uint64_t __w = _w;
			if (ln>0){
				size_t tpts = 1 << (ln - 1);
				while (tpts > 0) {
					for (size_t i = 0; i < tpts; i++, pos++) {
						pow_w[pos] = wi;
						pow_wp[pos] = ((uint64_t) pow_w[pos] << 32UL) / _pl;
						wi= (wi*__w)%_pl;
						//field().mulin(wi, __w);
					}
					wi = 1;
					__w = (__w * __w) % _pl;
					//field().mulin(__w, __w);
					tpts >>= 1;
				}
			}
		}


		Element getRoot() const {return _w;}
		Element getInvRoot() const {return _invw;}

		
		void FFT_DIF_Harvey (uint32_t *fft) {			
#ifdef __AVX2__
			FFT_DIF_Harvey_mod2p_iterative8x1_AVX(fft);
			if (n>=8){
				_vect256_t P;
				VEC256_SET_32(P,_pl);
				for (uint64_t i = 0; i < n; i += 8)
					reduce256_modp(fft+i,P);
				return;
			}
#else
			FFT_DIF_Harvey_mod2p_iterative4x2_SSE(fft);
#endif
			if (n >=4) {
				_vect128_t P;
				VEC128_SET_32(P,_pl);
				for (uint64_t i = 0; i < n; i += 4)
					reduce128_modp(fft+i,P);
			} else {
				for (uint64_t i = 0; i < n; i++)
					if (fft[i] >= _pl) fft[i] -= _pl;
			}
		}

		
		void FFT_DIT_Harvey (uint32_t *fft) {

#ifdef __AVX2__
			FFT_DIT_Harvey_mod4p_iterative8x1_AVX(fft);
			if (n>=8){
				_vect256_t P,P2;
				VEC256_SET_32(P, _pl);
				VEC256_SET_32(P2,_dpl);
				for (uint64_t i = 0; i < n; i += 8){
					reduce256_modp(&fft[i],P2);
					reduce256_modp(&fft[i],P);
				}
				return;
			}
#else
			FFT_DIT_Harvey_mod4p_iterative4x1_SSE(fft);
#endif
			if (n >=4) {
				_vect128_t P,P2;
				VEC128_SET_32(P,_pl);
				VEC128_SET_32(P2,_dpl);
				for (uint64_t i = 0; i < n; i += 4){
					reduce128_modp(&fft[i],P2);
					reduce128_modp(&fft[i],P);
				}
			} else {
				for (uint64_t i = 0; i < n; i++) {
					if (fft[i] >= (_pl << 1)) fft[i] -= (_pl << 1);
					if (fft[i] >= _pl) fft[i] -= _pl;
				}
			}
		}
		
		// FFT without conversion
		template <class T=Element>
		typename std::enable_if<std::is_same<T,uint32_t>::value>::type
		FFT_DIF (T *fft) {
			FFT_DIF_Harvey(fft);
		}
		template <class T=Element>
		typename std::enable_if<std::is_same<T,uint32_t>::value>::type
		FFT_DIT (T *fft) {
			FFT_DIT_Harvey(fft);
		}
		
		// FFT with conversion from Element to uint32_t		
		template <typename T=Element>
		typename std::enable_if<!std::is_same<T,uint32_t>::value>::type
		FFT_DIF (T *fft) {
			for(uint64_t i=0;i<n;i++)
				_data[i]=fft[i];
			FFT_DIF_Harvey(&_data[0]);
			for(uint64_t i=0;i<n;i++)
				fft[i]=_data[i];
			
		}
		template <typename T=Element>
		typename std::enable_if<!std::is_same<T,uint32_t>::value>::type
		FFT_DIT (T *fft) {
			for(uint64_t i=0;i<n;i++)
				_data[i]=fft[i];
			FFT_DIT_Harvey(&_data[0]);
			for(uint64_t i=0;i<n;i++)
				fft[i]=_data[i];			
		}
				
		/*
		 * Different implementations for the butterfly operations
		 */
		template <class T>
		inline void Butterfly_DIT_mod4p(T& A, T& B, const uint32_t& alpha, const uint32_t& alphap);
		template <class T>
		inline void Butterfly_DIF_mod2p(T& A, T& B, const uint32_t& alpha, const uint32_t& alphap);

		inline void reduce128_modp(uint32_t*, const __m128i&);
		inline void Butterfly_DIF_mod2p_4x1_SSE(uint32_t* ABCD, uint32_t* EFGH,const uint32_t* alpha, const uint32_t* alphap, const __m128i& P, const __m128i& P2);
		inline void Butterfly_DIF_mod2p_4x1_SSE_laststep(uint32_t* ABCD, uint32_t* EFGH, const __m128i& P2);
		inline void Butterfly_DIF_mod2p_4x2_SSE(uint32_t* , uint32_t* ,uint32_t* , uint32_t* ,
							const uint32_t* ,const uint32_t* ,const uint32_t* ,
							const uint32_t* ,const uint32_t* ,const uint32_t* ,
							const __m128i& P, const __m128i& P2);
		inline void Butterfly_DIF_mod2p_4x2_SSE_last2step(uint32_t* ABCD, uint32_t* EFGH, const __m128i& W,
								  const __m128i& Wp, const __m128i& P, const __m128i& P2);
		inline void Butterfly_DIT_mod4p_4x1_SSE(uint32_t* ABCD, uint32_t* EFGH, const uint32_t* alpha,
							const uint32_t* alphap,const __m128i& P, const __m128i& P2);
		inline void Butterfly_DIT_mod4p_4x2_SSE_first2step(uint32_t* ABCD, uint32_t* EFGH, const __m128i& W,
								   const __m128i& Wp, const __m128i& P, const __m128i& P2);
#ifdef __AVX2__
		inline void reduce256_modp(uint32_t*, const __m256i&);
		inline void Butterfly_DIF_mod2p_8x1_AVX(uint32_t* ABCD, uint32_t* EFGH, const uint32_t* alpha,const uint32_t* alphap,const __m256i& P, const __m256i& P2);
		inline void Butterfly_DIF_mod2p_8x3_AVX_last3step(uint32_t* ABCDEFGH, uint32_t* IJKLMNOP, const __m256i& alpha,const __m256i& alphap,
								  const __m256i& beta ,const __m256i& betap, const __m256i& P    ,const __m256i& P2);
		inline void Butterfly_DIT_mod4p_8x1_AVX(uint32_t* ABCD, uint32_t* EFGH, const uint32_t* alpha,const uint32_t* alphap,
			const __m256i& P, const __m256i& P2);
		inline void Butterfly_DIT_mod4p_8x3_AVX_first3step(uint32_t* ABCDEFGH, uint32_t* IJKLMNOP, const __m256i& alpha,const __m256i& alphap,
								   const __m256i& beta ,const __m256i& betap, const __m256i& P    ,const __m256i& P2);


#endif

		/*
		 * Different implementation of DIF/DIT with Harvey's trick
		 */

		void FFT_DIF_Harvey_mod2p_iterative    (Element *fft);
		void FFT_DIF_Harvey_mod2p_iterative2x2 (Element *fft);
		void FFT_DIF_Harvey_mod2p_iterative3x3 (Element *fft);
		void FFT_DIT_Harvey_mod4p_iterative2x2 (Element *fft);
		void FFT_DIT_Harvey_mod4p_iterative3x3 (Element *fft);
		// SIMD implementations follow
		void FFT_DIF_Harvey_mod2p_iterative4x1_SSE (uint32_t *fft);
		void FFT_DIF_Harvey_mod2p_iterative4x2_SSE (uint32_t *fft);
		void FFT_DIT_Harvey_mod4p_iterative4x1_SSE (uint32_t *fft);
#ifdef __AVX2__
		void FFT_DIF_Harvey_mod2p_iterative8x1_AVX (uint32_t *fft);
		void FFT_DIT_Harvey_mod4p_iterative8x1_AVX (uint32_t *fft);
#endif

	};
} // end of namespace LinBox

#include "linbox/algorithms/polynomial-matrix/polynomial-fft-transform.inl"
#include "linbox/algorithms/polynomial-matrix/polynomial-fft-transform-simd.inl"
#endif // __LINBOX_FFT_H



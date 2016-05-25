/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2016 Romain Lebreton, Pascal Giorgi
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


#ifndef __LINBOX_polynomial_fft_init_H
#define __LINBOX_polynomial_fft_init_H


#include <iostream>
#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"
#include "givaro/givinteger.h"
#include <fflas-ffpack/fflas/fflas_simd.h>

#ifndef ROUND_DOWN
#define ROUND_DOWN(x, s) ((x) & ~((s)-1))
#endif

// template<typename T>
// std::ostream& operator<<(std::ostream& os, const std::vector<T> &x){
// 	std::ostream_iterator<T> out_it (os,", ");
// 	std::copy ( x.begin(), x.end(), out_it );
// 	return os;
// }

#include "fflas-ffpack/utils/align-allocator.h"

#ifdef __LINBOX_USE_SIMD

//#include "linbox/algorithms/polynomial-matrix/simd.h"

#include "fflas-ffpack/fflas/fflas_simd.h"

#ifdef __LINBOX_USE_AVX2
/* 256 bits CODE HERE */
#define __LINBOX_USE_AVX2
// define 256 bits simd vector type
typedef __m256i  _vect256_t;
#endif
// define 128 bits simd vector type
typedef __m128i  _vect128_t;
#endif

namespace LinBox {

	enum SimdLevel {NOSIMD,SSE41,AVX,AVX2};

	struct SimdLevelFinder {
#ifdef __LINBOX_USE_AVX2
		const static SimdLevel simdlevel = AVX2;
#else
#ifdef __LINBOX_USE_AVX
		const static SimdLevel simdlevel = AVX;
#else
#ifdef __LINBOX_USE_SIMD
		const static SimdLevel simdlevel = SSE41;
#else
		const static SimdLevel simdlevel = NOSIMD;
#endif
#endif
#endif
	};

	// class to handle FFT transform over wordsize prime field Fp (p < 2^29)
	//	template <class Field, int SL = SimdLevelFinder::simdlevel>
	// TODO : A rendre générique / Simd si on doit faire des précalculs dans des Simd::vect_t
	template <class Field>
	class FFT_init {
	public:
		typedef typename Field::Element Element;
		typedef typename Field::Compute_t Compute_t;
		typedef typename Field::Residu_t Residu_t;

		const Field                *fld;
		Residu_t              _pl, _dpl;
		uint64_t                      n;
		size_t                       ln;
		//Compute_t                  _logp;
		//Compute_t                     _I;
		//double                    _pinv;
		Element                      _w;
		Element                   _invw;
		// Du type qui est donné aux Butterfly
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

		Element find_gen (Residu_t _m, uint64_t _val2p) {
			// find a primitive 2^k root of unity where
			// _p - 1 = 2^val2p * m
			srand((unsigned int) time(NULL));
			Element y,z;
			uint64_t j;
			Element _gen;
			for (;;) {
				fld->init(_gen,rand());
				fld->init(z, 1);
				for (Residu_t i=0; i < _m; ++i) fld->mulin(z,_gen); // z = z*_gen;
				if (z == 1) continue;
				// _gen^i =/ 1 pour 0 <= i < m
				_gen = z;
				j = 0;
				do {
					y = z;
					fld->mul(z,y,y); // z = y * y;
					j++;
				} while (j != _val2p && z != 1);
				if (j == _val2p)
					break;
			}
			return _gen;
		}

		FFT_init (const Field& fld2, size_t ln2, Element w = 0)
			: fld (&fld2), n ((1UL << ln2)), ln (ln2), pow_w(n - 1), pow_wp(n - 1), _data(n) {
			_pl = fld->characteristic();
			_p  = fld->characteristic();

			linbox_check((_pl >> 29) == 0 ); // 8*p < 2^31 for Harvey's butterflies
			_dpl = (_pl << 1);
			//_pinv = 1 / (double) _pl;

			Givaro::Timer chrono;
			chrono.start();

			uint64_t _val2p = 0;
			Residu_t     _m = _pl;
			_m = _pl - 1;
			while ((_m & 1) == 0) {
				_m >>= 1;
				_val2p++;
			}

			linbox_check(ln <= _val2p);      // Otherwise no 2 _ln roots of unity

			if (w == 0){   // find a pseudo 2^lpts-th primitive root of unity
				//_I = (1L << (_logp << 1)) / _pl;
				Element _gen = find_gen (_m, _val2p);
				_w = Givaro::powmod(_gen, 1UL<<(_val2p-ln), _pl);
			}
			else {
				_w = w;
			}
			chrono.clear();
			chrono.start();

			// compute w^(-1) mod p = w^(2^lpts - 1)
			_invw = Givaro::powmod(_w, (1UL<<ln) - 1, _pl);

			size_t pos = 0;
			//uint64_t wi = 1;
			Element wi = 1;
			Element __w = _w;
			size_t _logp = Givaro::Integer(_pl).bitsize() - 1;
			uint32_t BAR= (Givaro::Integer(1)<<(32+_logp))/Givaro::Integer(_pl);
			uint32_t Q;
			//cout<<"log Bar: "<<Integer(BAR).bitsize()<<endl;
			if (ln>0){
#if 0
				//#ifdef __LINBOX_NO_SIMD
				size_t tpts = 1 << (ln - 1);
				while (tpts > 0) {
					for (size_t i = 0; i < tpts; i++, pos++) {
						pow_w[pos] = wi;
						pow_wp[pos] = ((uint64_t) pow_w[pos] << 32UL) / _pl;
						fld->mulin(wi,__w); //wi= (wi*__w)%_pl;
					}
					wi = 1;
					fld->mulin(__w,__w); //__w = (__w * __w) % _pl;
					//field().mulin(__w, __w);
					tpts >>= 1;
				}
#else
				using simd=Simd<uint32_t>;
//				using vect_t =typename simd::vect_t;

				size_t tpts = 1 << (ln - 1);
				size_t i=0;
				for( ;i<std::min(simd::vect_size+1, tpts);i++,pos++){
					pow_w[pos] = wi;
					pow_wp[pos] = ((uint64_t) pow_w[pos] << 32UL) / _pl;
					wi= ((uint64_t)wi*(uint64_t)__w)%_pl;
				}
				/*
				//Precomputation with SIMD - No time improvement ?
				vect_t wp_vect, Q_vect,BAR_vect,w_vect,pow_w_vect,pow_wp_vect, pl_vect;
				BAR_vect= simd::set1(BAR);
				wp_vect = simd::set1(pow_wp[simd::vect_size]);
				w_vect  = simd::set1(pow_w[simd::vect_size]);
				pl_vect = simd::set1(_pl);
				for (; i < ROUND_DOWN(tpts,simd::vect_size);
					 i+=simd::vect_size,pos+=simd::vect_size) {
					pow_w_vect  = simd::loadu((int32_t*)pow_w.data()+pos-simd::vect_size);
					Q_vect=simd::mulhi(pow_w_vect,wp_vect);
					pow_w_vect = simd::sub(simd::mullo(pow_w_vect,w_vect),simd::mullo(Q_vect,pl_vect));
					pow_w_vect=simd::sub(pow_w_vect, simd::vandnot(simd::greater(pow_w_vect,pl_vect),pl_vect));
					simd::storeu((int32_t*)pow_w.data()+pos,pow_w_vect);
					pow_wp_vect= simd::mulhi(simd::sll(pow_w_vect,32-_logp),BAR_vect);
					simd::storeu((int32_t*)pow_wp.data()+pos,pow_wp_vect);
				}
				*/
				for( ;i<tpts;i++,pos++){
					pow_w[pos] = wi;
					pow_wp[pos]= (((uint64_t)wi*BAR)>>_logp);
					Q= ((uint64_t)wi*pow_wp[1])>>32;
					wi= (uint32_t)(wi*__w - Q*_pl);
					wi-=(wi>=_pl?_pl:0);
				}

				for(size_t k=2;k<=tpts;k<<=1)
					for(size_t i=0;i<tpts;i+=k,pos++){
						pow_w[pos]  = pow_w[i];
						pow_wp[pos] = pow_wp[i];
					}
#endif

			}
			chrono.stop();
			//cout<<"FFT: table="<<chrono<<endl;
		}


		Element getRoot() const {return _w;}
		Element getInvRoot() const {return _invw;}

	}; //FFT_init

}

#endif // __LINBOX_polynomial_fft_init_H

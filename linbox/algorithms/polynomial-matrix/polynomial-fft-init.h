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

#include <linbox/util/timer.h>

#ifndef ROUND_DOWN
#define ROUND_DOWN(x, s) ((x) & ~((s)-1))
#endif

#include "fflas-ffpack/utils/align-allocator.h"

#ifdef __LINBOX_USE_SIMD

#include "fflas-ffpack/fflas/fflas_simd.h"
#include "linbox/algorithms/polynomial-matrix/simd-additional-functions.h"

#endif

namespace LinBox {

	template<typename Field, typename simd = Simd<typename Field::Element>,
			 bool is_integral = std::is_integral<typename Field::Element>::value >
	class InitPowers {
		template<typename Allocator>
		static void init_powers (std::vector<typename Field::Element,Allocator>& pow_w, std::vector<typename Field::Element,Allocator>& pow_wp,
								 Field* fld, typename Field::Element _w, size_t ln);
	};

	template<typename Field>
	class InitPowers<Field, NoSimd<typename Field::Element>, true> {
	public :
		template<typename Allocator>
		static void init_powers (std::vector<typename Field::Element,Allocator>& pow_w, std::vector<typename Field::Element,Allocator>& pow_wp,
								 const Field* const fld, const typename Field::Element _w, const size_t ln) {
			using Element = typename Field::Element;
			using Compute_t = typename Field::Compute_t;
			size_t pos = 0;
			Element wi = 1;

			// Precomp Quo(2^32,p)
			Compute_t invp; fld->precomp_p(invp);

			if (ln>0){
				size_t tpts = 1 << (ln - 1);
				size_t i=0;
				// Precompute pow_wp[1] for faster mult by pow_w[1]
				for( ;i<std::min((size_t) 2, tpts);i++,pos++){
					pow_w[pos] = wi;
					// Fake conversion since precomp_b will be used as a Compute_t in mul_precomp_b
					Compute_t mulprecomp;
					fld->precomp_b(mulprecomp, wi, invp);
					pow_wp[pos] = static_cast<Element>(mulprecomp);

					fld->mulin(wi, _w);
				}

				// Use pow_wp[1] for speed-up mult by pow_w[1]
				for( ;i<tpts;i++,pos++){
					pow_w[pos] = wi;

					// Fake conversion since precomp_b will be used as a Compute_t in mul_precomp_b
					Compute_t mulprecomp;
					fld->precomp_b(mulprecomp, wi, invp);
					pow_wp[pos] = static_cast<Element>(mulprecomp);

#ifndef NDEBUG
					Compute_t temp2;
					fld->precomp_b(temp2, wi);
					if (temp2 != mulprecomp) {
						std::cout << "Erreur de précalcul : " << ((uint64_t) mulprecomp) << "\t\t" << ((uint64_t) temp2) << std::endl;
						std::cout << "p := " << ((uint64_t) fld->characteristic()) << std::endl
							 << "invp := " << ((uint64_t) invp) << std::endl
							 << "b := " << ((uint64_t) wi) << std::endl;
						throw std::string("Message d'erreur");
					}
#endif

					fld->mul_precomp_b(wi, wi, _w, static_cast<Compute_t>(pow_wp[1]));
				}

				// Other pow_w elements can be read from previously computed pow_w
				for(size_t k=2;k<=tpts;k<<=1)
					for(size_t i=0;i<tpts;i+=k,pos++){
						pow_w[pos]  = pow_w[i];
						pow_wp[pos] = pow_wp[i];
					}
			}
		}

	};

	template<typename Field, typename simd>
	class InitPowers<Field, simd, true> {
	public :
		template<typename Allocator>
		static void init_powers (std::vector<typename Field::Element,Allocator>& pow_w, std::vector<typename Field::Element,Allocator>& pow_wp,
								 const Field* const fld, const typename Field::Element _w, const size_t ln) {
			using Element = typename Field::Element;
			using Compute_t = typename Field::Compute_t;
			using Residu_t = typename Field::Residu_t;

			size_t pos = 0;
			Element wi = 1;

			// Precomp Quo(2^32,p)
			Compute_t invp; fld->precomp_p(invp);

			if (ln>0){
				using vect_t =typename simd::vect_t;

				size_t tpts = 1 << (ln - 1);
				size_t i=0;

				// Precompute pow_wp[1] for faster mult by pow_w[1]
				for( ; i<std::min(simd::vect_size+1, tpts); i++, pos++){
					pow_w[pos] = wi;

					// Fake conversion since precomp_b will be used as a Compute_t in mul_precomp_b
					Compute_t mulprecomp;
					fld->precomp_b(mulprecomp, wi, invp); //(((Compute_t)wi*invp)>>(fld->_bitsizep));
					pow_wp[pos] = static_cast<Element>(mulprecomp);

					fld->mulin(wi, _w);
				}

				if (i < tpts){

					Compute_t invp;
					fld->precomp_p(invp);

					size_t _bitsizep = 0;
					Residu_t p = fld->characteristic();
					Residu_t __p = p;
					while (__p != 0) {
						_bitsizep++;
						__p >>= 1;
					}


					Compute_t pow2 = Compute_t(1) << (8*sizeof(Element));
					Residu_t pow2divp = static_cast<Residu_t>( pow2 / static_cast<Compute_t>(p));
					Residu_t pow2modp = static_cast<Residu_t>( pow2 - static_cast<Compute_t>(pow2divp) * static_cast<Compute_t>(p));
					Compute_t pow2modp_precomp;
					fld->precomp_b(pow2modp_precomp, static_cast<Element>(pow2modp), invp);

					vect_t pow_w_vect;
					vect_t one_vect = simd::set1(1);
					vect_t w_vect  = simd::set1(pow_w[simd::vect_size]);
					vect_t wp_vect = simd::set1(pow_wp[simd::vect_size]);
					vect_t p_vect = simd::set1(p);
					vect_t pow2divp_vect = simd::set1(pow2divp);
					vect_t pow2modp_vect = simd::set1(pow2modp);
					vect_t pow2modp_precomp_vect = simd::set1(pow2modp_precomp);

					// Compute tpts first elements : 1, w, w^2, ..., w^{K/2-1}
					// Use pow_wp[vect_size-1] for speed-up mult by pow_w[vect_size-1]
					// Here we assume that vect_size divides tpts
					Element* posw  = ((Element*)pow_w .data())+pos;
					Element* poswp = ((Element*)pow_wp.data())+pos;
					for (i = simd::vect_size, pos = simd::vect_size; i < tpts; i+=simd::vect_size, pos+=simd::vect_size, posw+=simd::vect_size, poswp+=simd::vect_size) {

						// TODO : call MemoryOp for loadu and storeu
						//pow_w_vect  = simd::loadu(posw - simd::vect_size);
						pow_w_vect  = simd::loadu(&pow_w[pos- simd::vect_size]);

						pow_w_vect = mul_mod<simd>(pow_w_vect, w_vect, p_vect, wp_vect);
						pow_w_vect = reduce<simd>(pow_w_vect, p_vect);

						// Precompute Floor(2^n * w^i / p) for faster mult by w^i modulo p
						vect_t q = simd::mulhi(pow_w_vect,pow2modp_precomp_vect);
						vect_t c = simd::mullo(pow_w_vect,pow2modp_vect);
						vect_t t = simd::mullo(q,p_vect);
						t = simd::sub(c,t);

						// if (a*b - q *p) >= p then q++
						vect_t f = simd::greater(p_vect,t);
						q = simd::add(q, simd::vandnot(one_vect,f));

						// q += pow2divp * w^i
						q = simd::add(q, simd::mul(pow2divp_vect, pow_w_vect));

						//						simd::storeu(posw,pow_w_vect);
						//						simd::storeu(poswp,pow_wp_vect);
						simd::storeu(&pow_w[pos] ,pow_w_vect);
						simd::storeu(&pow_wp[pos],q);
					}

				}

				// Other pow_w elements can be read from previously computed pow_w
				for(size_t k=2;k<=tpts;k<<=1)
					for(size_t i=0;i<tpts;i+=k,pos++){
						pow_w[pos]  = pow_w[i];
						pow_wp[pos] = pow_wp[i];
					}

#if 0
				std::std::cout << "Check precomputations : pow_w, pow_wp \n";
				std::std::cout << "[";
				for (size_t i=0; i < tpts; i++) std::std::cout << pow_w[i] << ", ";
				std::std::cout << "]\n";
				std::std::cout << "[";
				for (size_t i=0; i < tpts; i++) std::std::cout << pow_wp[i] << ", ";
				std::std::cout << "]\n\n";
#endif
			}
		}
	};

	template<typename Field>
	class InitPowers<Field, NoSimd<typename Field::Element>, false> {
	public :
		template<typename Allocator>
		static void init_powers (std::vector<typename Field::Element,Allocator>& pow_w, std::vector<typename Field::Element,Allocator>& pow_wp,
								 const Field* const fld, const typename Field::Element _w, const size_t ln) {
			using Element = typename Field::Element;
			//		using Compute_t = typename Field::Compute_t;
			size_t pos = 0;
			Element wi = 1;

			if (ln>0){
				size_t tpts = 1 << (ln - 1);

				for(size_t i=0; i<tpts;i++,pos++){
					pow_w[pos] = wi;
					fld->mulin(wi,_w);
				}

				// Other pow_w elements can be read from previously computed pow_w
				for(size_t k=2;k<=tpts;k<<=1)
					for(size_t i=0;i<tpts;i+=k,pos++){
						pow_w[pos]  = pow_w[i];
					}

			}
		}
	};

	template<typename Field, typename simd>
	class InitPowers<Field, simd, false> {
		// TODO : this class is not yet implemented as it should be
	public :
		template<typename Allocator>
		static void init_powers (std::vector<typename Field::Element,Allocator>& pow_w, std::vector<typename Field::Element,Allocator>& pow_wp,
								 const Field* const fld, const typename Field::Element _w, const size_t ln) {
			// TODO !! No Simd so far
			// Problem is to differentiate between Modular<double> and ModularExtended<double>

			using Element = typename Field::Element;
			//		using Compute_t = typename Field::Compute_t;
			size_t pos = 0;
			Element wi = 1;

			if (ln>0){
				size_t tpts = 1 << (ln - 1);

				for(size_t i=0; i<tpts;i++,pos++){
					pow_w[pos] = wi;
					fld->mulin(wi,_w);
				}

				// Other pow_w elements can be read from previously computed pow_w
				for(size_t k=2;k<=tpts;k<<=1)
					for(size_t i=0;i<tpts;i+=k,pos++){
						pow_w[pos]  = pow_w[i];
					}

			}
		}
	};



	// class to handle FFT transform over wordsize prime field Fp (p < 2^29)
	//	template <class Field, int SL = SimdLevelFinder::simdlevel>
	// TODO : A rendre générique / Simd si on doit faire des précalculs dans des Simd::vect_t
	template <class Field>
	class FFT_init {
	public:
		using Element = typename Field::Element;
		using Compute_t = typename Field::Compute_t;
		using Residu_t = typename Field::Residu_t;

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
		typedef std::vector<Element,AlignedAllocator<Element, Alignment::DEFAULT> > VECT;
		VECT    pow_w;
		VECT   pow_wp; // Precomputations in shoup
		VECT pow_w_bis, pow_wp_bis;
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

		void binpow (Element& res, const Field& F, const Element& x, const size_t& e) {
			Element x2;
			if (e == 0) {
				F.init(res,F.one);
				return;
			}
			if (e == 1) {
				res = x;
				return;
			}
			F.mul(x2,x,x);
			binpow(res, F, x2, e >> 1);
			if ((e & 1) == 1)
				F.mulin(res,x);
		}

		Element find_gen (Residu_t _m, uint64_t _val2p) {
			// find a primitive 2^k root of unity where
			// _p - 1 = 2^val2p * m
			srand((unsigned int) time(NULL));
			Element y,z;
			uint64_t j;
			Element _gen;

			Timer chrono;
			chrono.start();
			size_t iter = 0;
			for (;;) {
				iter++;
				fld->init(_gen,rand());
				fld->init(z, 1);

				//				for (Residu_t i=0; i < _m; ++i) fld->mulin(z,_gen); // z = z*_gen;
				binpow(z, *fld, _gen, (size_t) _m);

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
			chrono.stop();
			std::cout << "Time to find gen : " << chrono.usertime() << " in " << iter << " iterations\n";

			return _gen;
		}

		FFT_init (const Field& fld2, size_t ln2, Element w = 0)
			: fld (&fld2), n ((1UL << ln2)), ln (ln2), pow_w(n - 1), pow_wp(n - 1)
			, pow_w_bis(n-1) , pow_wp_bis(n-1) {
			_pl = fld->characteristic();
			_p  = fld->characteristic();

			linbox_check(_pl <= (field().maxCardinality() >> 3)); // 8*p <= field()->maxCardinality() for Harvey's butterflies
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
				//				_w = Givaro::powmod(_gen, (int64_t)1<<(_val2p-ln), _pl);
				binpow(_w, *fld, _gen, (int64_t)1<<(_val2p-ln));
			}
			else {
				_w = w;
			}

			// compute w^(-1) mod p = w^(2^lpts - 1)
			// _invw = Givaro::powmod(_w, ((int64_t)1<<ln) - 1, _pl);
			binpow(_invw, *fld, _w, ((int64_t)1<<ln) - 1);

			std::cout.precision(2);
			std::cout.width(10);
			std::cout<<std::scientific;

			chrono.stop();
			std::cout << "Time to find gen in FFT_init : " << chrono.usertime() << std::endl;

			chrono.clear();
			chrono.start();
			size_t ct = 0;
			while (chrono.realElapsedTime() < 0.1){
				InitPowers<Field, Simd<Element> >::init_powers(pow_w_bis, pow_wp_bis, fld, _w, ln);
				ct++;
			}
			chrono.stop();
			std::cout << "Time to init_powers in FFT_init<Simd> : " << (chrono.usertime()/ct) << std::endl;

#ifndef NDEBUG
			chrono.clear();
			chrono.start();
			ct = 0;
			while (chrono.realElapsedTime() < 0.1){
				InitPowers<Field, NoSimd<Element> >::init_powers(pow_w, pow_wp, fld, _w, ln);
				ct++;
			}
			chrono.stop();
			std::cout << "Time to init_powers in FFT_init<NoSimd> : " << (chrono.usertime()/ct) << std::endl;


			std::cout << "Test pow_w & pow_wp : " << ((pow_w_bis == pow_w)?"OK ":"KO!!! ") << ((pow_wp_bis == pow_wp)?"OK ":"KO!!! ") << std::endl << std::endl;
			if (((pow_w_bis != pow_w) || (pow_wp_bis != pow_wp)) && (ln2 < 6)) {
				std::cout << "pow_w : " << std::endl
					 << pow_w << std::endl << std::endl
					 << "pow_w_bis : " << std::endl
					 << pow_w_bis << std::endl << std::endl;


				std::cout << "pow_wp : " << std::endl
					 << pow_wp << std::endl << std::endl
					 << "pow_wp_bis : " << std::endl
					 << pow_wp_bis << std::endl << std::endl;

				throw std::string("Message d'erreur");
			}
#endif

		}


		Element getRoot() const {return _w;}
		Element getInvRoot() const {return _invw;}

	}; //FFT_init

}

#endif // __LINBOX_polynomial_fft_init_H

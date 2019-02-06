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


#ifndef __LINBOX_polynomial_fft_algorithms_H
#define __LINBOX_polynomial_fft_algorithms_H

#include <iostream>
#include "linbox/linbox-config.h"
#include "fflas-ffpack/fflas/fflas_simd.h"
#include "linbox/algorithms/polynomial-matrix/simd-additional-functions.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-fft-init.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-fft-butterflies.h"

#define IS_INTEGRAL \
    typename std::enable_if<std::is_integral<typename Field::Element>::value>::type
#define IS_FLOATING \
    typename std::enable_if<std::is_floating_point<typename Field::Element>::value>::type

namespace LinBox {

	template<typename Field, typename simd = Simd<typename Field::Element>, uint8_t vect_size = simd::vect_size, typename Enable = void>
	class FFT_algorithms : public FFT_butterflies<Field, simd, vect_size> {
	public:
		using Element = typename Field::Element;
		FFT_algorithms(const FFT_init<Field>& f_i) : FFT_butterflies<Field, simd, vect_size>(f_i) {}
		//void DIF_mod2p (Element *fft);
		//void DIT_mod4p (Element *fft);
		void DIF (Element *fft);
		void DIT (Element *fft);
	}; // FFT_algorithms

	template<typename Field>
	class FFT_algorithms<Field, NoSimd<typename Field::Element>, 1, IS_INTEGRAL> : public FFT_butterflies<Field, NoSimd<typename Field::Element>, 1> {
	public:
		using Element = typename Field::Element;

		FFT_algorithms(const FFT_init<Field>& f_i) : FFT_butterflies<Field, NoSimd<typename Field::Element>, 1>(f_i) {}

		void DIF_mod2p (Element *fft) {
			Element * tab_w = &(this->pow_w) [0];
			Element * tab_wp= &(this->pow_wp)[0];
			for (size_t w = this->n >> 1, f = 1; w != 0; tab_w+=w, tab_wp+=w, f <<= 1, w >>= 1){
				// w : witdh of butterflies
				// f : # families of butterflies
				for (size_t i = 0; i < f; i++)
					for (size_t j = 0; j < w; j++)
						this->Butterfly_DIF_mod2p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], tab_w[j],tab_wp[j]);
			}
		}

		void DIT_mod4p (Element *fft) {
			Element * tab_w = &(this->pow_w) [this->n-2];
			Element * tab_wp= &(this->pow_wp)[this->n-2];
			for (size_t w = 1, f = this->n >> 1; f >= 1; w <<= 1, f >>= 1, tab_w-=w, tab_wp-=w)
				for (size_t i = 0; i < f; i++)
					for (size_t j = 0; j < w; j++)
						this->Butterfly_DIT_mod4p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], tab_w[j], tab_wp[j]);
		}

		void DIF (Element *fft) {
			DIF_mod2p(fft);
			//DIF_mod2p_iterative(fft);
			for (uint64_t i = 0; i < this->n; i++) {
				//				if (fft[i] >= (_pl << 1)) fft[i] -= (_pl << 1);
				if (fft[i] >= this->_pl) fft[i] -= this->_pl;
			}
		}

		void DIT (Element *fft) {
			DIT_mod4p(fft);
			//DIF_mod2p_iterative(fft);
			for (uint64_t i = 0; i < this->n; i++) {
				if (fft[i] >= (this->_pl << 1)) fft[i] -= (this->_pl << 1);
				if (fft[i] >= this->_pl) fft[i] -= this->_pl;
			}
		}

	}; // FFT_algorithms<Field, NoSimd<typename Field::Element>, 1>

	template<typename Field, typename simd>
	class FFT_algorithms<Field, simd, 4, IS_INTEGRAL> : public FFT_butterflies<Field, simd, 4> {
	public:
		using Element = typename Field::Element;
		using Compute_t = typename Field::Compute_t;
		using Residu_t = typename Field::Residu_t;
		using SimdMod = SimdModular<Field, simd>;
		using vect_t = typename simd::vect_t;

		FFT_algorithms(const FFT_init<Field>& f_i) : FFT_butterflies<Field, simd, 4>(f_i) {
			linbox_check(simd::vect_size == 4);
		}

		void DIF_mod2p (Element *fft) {
			const uint64_t& n = this->n;
			const Residu_t& _pl = this->_pl;
			const Residu_t& _dpl = this->_dpl;

			vect_t P,P2;
			P  = simd::set1(_pl);
			P2 = simd::set1(_dpl);
			Element * tab_w = &(this->pow_w) [0];
			Element * tab_wp= &(this->pow_wp)[0];
			size_t w, f;
			for (w = n >> 1, f = 1; w >= 4; tab_w+=w, tab_wp+=w, w >>= 1, f <<= 1){
				// w : witdh of butterflies
				// f : # families of butterflies
				for (size_t i = 0; i < f; i++)
					for (size_t j = 0; j < w; j+=4)

#define A0 &fft[0] +  (i << 1)   *w+ j
#define A4 &fft[0] + ((i << 1)+1)*w+ j
						this->Butterfly_DIF_mod2p(A0,A4, tab_w+j,tab_wp+j,P,P2);
#undef A0
#undef A4
				//std::cout<<fft<<std::endl;
			}
			// Last two steps
			if (n >= 8) {
				vect_t W,Wp;
				W = simd::set1 (tab_w [1]);
				Wp= simd::set1 (tab_wp[1]);

				for (size_t i = 0; i < f; i+=2)
#define A0 &fft[0] +  (i << 2)
#define A4 &fft[0] + ((i << 2)+4)
					this->Butterfly_DIF_mod2p_laststeps(A0,A4,W,Wp,P,P2);
				//std::cout<<fft<<std::endl;
#undef A0
#undef A4
			} else {
				FFT_algorithms<Field, NoSimd<Element>, 1> fft_algo_1 (FFT_init<Field> (this->field(),this->ln,this->getRoot()));

				for (; w >= 1; tab_w+=w, tab_wp+=w, w >>= 1, f <<= 1)
					for (size_t i = 0; i < f; i++)
						for (size_t j = 0; j < w; j++)
							fft_algo_1.Butterfly_DIF_mod2p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], tab_w[j], tab_wp[j]);
			}
		}

		void DIT_mod4p (Element *fft) {
			const uint64_t& n = this->n;
			const Residu_t& _pl = this->_pl;
			const Residu_t& _dpl = this->_dpl;

			vect_t P,P2;
			P = simd::set1(_pl);
			P2 = simd::set1(_dpl);
			// First two steps
			if (n >= 8) {
				vect_t W,Wp;
				W = simd::set1 ((this->pow_w) [n-3]);
				Wp= simd::set1 ((this->pow_wp)[n-3]);

				for (size_t i = 0; i < n; i+=8)
					this->Butterfly_DIT_mod4p_firststeps(&fft[i],&fft[i+4],W,Wp,P,P2);

				Element * tab_w = &(this->pow_w) [n-8];
				Element * tab_wp= &(this->pow_wp)[n-8];
				for (size_t w = 4, f = n >> 3; f >= 1; w <<= 1, f >>= 1, tab_w-=w, tab_wp-=w){
					// w : witdh of butterflies
					// f : # families of butterflies
					for (size_t i = 0; i < f; i++)
						for (size_t j = 0; j < w; j+=4)
#define A0 &fft[0] +  (i << 1)   *w+ j
#define A4 &fft[0] + ((i << 1)+1)*w+ j
							this->Butterfly_DIT_mod4p(A0,A4, tab_w+j,tab_wp+j,P,P2);

#undef A0
#undef A4

				}
			} else {
				FFT_algorithms<Field, NoSimd<Element>, 1> fft_algo_1 (FFT_init<Field> (this->field(),this->ln,this->getRoot()));

				Element * tab_w = &(this->pow_w) [n-2];
				Element * tab_wp= &(this->pow_wp)[n-2];
				for (size_t w = 1, f = n >> 1; f >= 1; w <<= 1, f >>= 1, tab_w-=w, tab_wp-=w)
					for (size_t i = 0; i < f; i++)
						for (size_t j = 0; j < w; j++)
							fft_algo_1.Butterfly_DIT_mod4p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], tab_w[j], tab_wp[j]);
			}
		}

		void DIF (Element *fft) {
			DIF_mod2p(fft);

			if (this->n >= 4) {
				vect_t P;
				P  = simd::set1(this->_pl);
				for (uint64_t i = 0; i < this->n; i += 4)
					SimdMod::reduce (&fft[i],P);
				return;
			} else {
				for (uint64_t i = 0; i < this->n; i++)
					if (fft[i] >= this->_pl) fft[i] -= this->_pl;
			}
		}

		void DIT (Element *fft) {
			DIT_mod4p(fft);

			if (this->n >= 4) {
				vect_t P,P2;
				P  = simd::set1(this->_pl);
				P2 = simd::set1(this->_dpl);
				for (uint64_t i = 0; i < this->n; i += 4){
					SimdMod::reduce (&fft[i],P2);
					SimdMod::reduce (&fft[i],P);
				}
				return;

			} else {
				for (uint64_t i = 0; i < this->n; i++) {
					if (fft[i] >= (this->_pl << 1)) fft[i] -= (this->_pl << 1);
					if (fft[i] >= this->_pl) fft[i] -= this->_pl;
				}
			}
		}

	}; // FFT_algorithms<Field, NoSimd<typename Field::Element>, 4>

	template<typename Field, typename simd>
	class FFT_algorithms<Field, simd, 8, IS_INTEGRAL> : public FFT_butterflies<Field, simd, 8> {
	public:
		using Element = typename Field::Element;
		using SimdMod = SimdModular<Field, simd>;
		using vect_t = typename simd::vect_t;

		FFT_algorithms(const FFT_init<Field>& f_i) : FFT_butterflies<Field, simd, 8>(f_i) {
			linbox_check(simd::vect_size == 8);
		}

		void DIF_mod2p (Element *fft) {
			vect_t P,P2;
			P = simd::set1(this->_pl);
			P2 = simd::set1(this->_dpl);

			Element * tab_w = &(this->pow_w) [0];
			Element * tab_wp= &(this->pow_wp)[0];
			size_t w, f;
			for (w = this->n >> 1, f = 1; w >= 8; tab_w+=w, tab_wp+=w, w >>= 1, f <<= 1){
				// w : witdh of butterflies
				// f : # families of butterflies
				for (size_t i = 0; i < f; i++)
					for (size_t j = 0; j < w; j+=8)
#define A0 &fft[0] +  (i << 1)   *w+ j
#define A4 &fft[0] + ((i << 1)+1)*w+ j
						this->Butterfly_DIF_mod2p(A0,A4, tab_w+j,tab_wp+j,P,P2);

#undef A0
#undef A4
				//std::cout<<fft<<std::endl;
			}
			// Last three steps
			if (this->n >= 16) {
				vect_t alpha,alphap,beta,betap;
				Element tmp[8];
				tmp[0]=tmp[4]=tab_w[0];
				tmp[1]=tmp[5]=tab_w[1];
				tmp[2]=tmp[6]=tab_w[2];
				tmp[3]=tmp[7]=tab_w[3];
				alpha = MemoryOp<Element,simd>::load(tmp);
				tmp[0]=tmp[4]=tab_wp[0];
				tmp[1]=tmp[5]=tab_wp[1];
				tmp[2]=tmp[6]=tab_wp[2];
				tmp[3]=tmp[7]=tab_wp[3];
				alphap = MemoryOp<Element,simd>::load(tmp);
				beta = simd::set1(tab_w [5]);
				betap = simd::set1(tab_wp [5]);

				for (size_t i = 0; i < f; i+=2)
#define A0 &fft[0] + (i << 3)
#define A4 &fft[0] + (i << 3)+8
					this->Butterfly_DIF_mod2p_laststeps(A0,A4,alpha,alphap,beta,betap,P,P2);
#undef A0
#undef A4
				//std::cout<<fft<<std::endl;
			} else {
				// TODO : improve ?
				//FFT_algorithms<Field, NoSimd<Element>, 1> fft_algo_1 ((FFT_init<Field>) *this);
				FFT_algorithms<Field, NoSimd<Element>, 1> fft_algo_1 (FFT_init<Field> (this->field(),this->ln,this->getRoot()));

				for (; w >= 1; tab_w+=w, tab_wp+=w, w >>= 1, f <<= 1)
					for (size_t i = 0; i < f; i++)
						for (size_t j = 0; j < w; j++)
							fft_algo_1.Butterfly_DIF_mod2p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], tab_w[j], tab_wp[j]);
			}
		}

		void DIT_mod4p (Element *fft) {
			const auto &pow_w = this->pow_w;
			const auto &pow_wp = this->pow_wp;
			const uint64_t &n = this->n;

			vect_t P,P2;
			P = simd::set1(this->_pl);
			P2 = simd::set1(this->_dpl);

			// first three steps
			if (n >= 16) {
				vect_t alpha,alphap,beta,betap;
				alpha = simd::set1((pow_w)[n-3]);
				alphap = simd::set1((pow_wp)[n-3]);
				Element tmp[8];
				tmp[0]=tmp[4]=(pow_w)[n-8];
				tmp[1]=tmp[5]=(pow_w)[n-7];
				tmp[2]=tmp[6]=(pow_w)[n-6];
				tmp[3]=tmp[7]=(pow_w)[n-5];
				beta = MemoryOp<Element,simd>::load(tmp);
				tmp[0]=tmp[4]=(pow_wp)[n-8];
				tmp[1]=tmp[5]=(pow_wp)[n-7];
				tmp[2]=tmp[6]=(pow_wp)[n-6];
				tmp[3]=tmp[7]=(pow_wp)[n-5];
				betap = MemoryOp<Element,simd>::load(tmp);
				for (uint64_t i = 0; i < n; i+=16) {
					this->Butterfly_DIT_mod4p_firststeps(&fft[i],&fft[i+8],alpha,alphap,beta,betap,P,P2);
				}
				const Element * tab_w = &(pow_w) [n-16];
				const Element * tab_wp= &(pow_wp)[n-16];
				for (size_t w = 8, f = n >> 4; f >= 1; w <<= 1, f >>= 1, tab_w-=w, tab_wp-=w){
					// w : witdh of butterflies
					// f : # families of butterflies
					for (size_t i = 0; i < f; i++)
						for (size_t j = 0; j < w; j+=8) {
#define A0 &fft[0] +  (i << 1)   *w+ j
#define A4 &fft[0] + ((i << 1)+1)*w+ j
							this->Butterfly_DIT_mod4p(A0,A4, tab_w+j,tab_wp+j,P,P2);
#undef A0
#undef A4
						}
				}
			} else {

				FFT_algorithms<Field, NoSimd<Element>, 1> fft_algo_1 (FFT_init<Field> (this->field(),this->ln,this->getRoot()));

				const Element * tab_w = &(pow_w) [n-2];
				const Element * tab_wp= &(pow_wp)[n-2];
				for (size_t w = 1, f = n >> 1; f >= 1; w <<= 1, f >>= 1, tab_w-=w, tab_wp-=w)
					for (size_t i = 0; i < f; i++)
						for (size_t j = 0; j < w; j++)
							fft_algo_1.Butterfly_DIT_mod4p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], tab_w[j], tab_wp[j]);
			}
		}

		void DIF (Element *fft) {
			DIF_mod2p(fft);

			if (this->n >= 8) {
				vect_t P;
				P  = simd::set1(this->_pl);
				for (uint64_t i = 0; i < this->n; i += 8){
					SimdMod::reduce (&fft[i],P);
				}
				return;

			} else {
				for (uint64_t i = 0; i < this->n; i++)
					if (fft[i] >= this->_pl) fft[i] -= this->_pl;
			}
		}

		void DIT (Element *fft) {
			DIT_mod4p(fft);

			if (this->n >= 8) {
				vect_t P,P2;
				P  = simd::set1(this->_pl);
				P2 = simd::set1(this->_dpl);
				for (uint64_t i = 0; i < this->n; i += 8){
					SimdMod::reduce (&fft[i],P2);
					SimdMod::reduce (&fft[i],P);
				}
				return;

			} else {
				for (uint64_t i = 0; i < this->n; i++) {
					if (fft[i] >= (this->_pl << 1)) fft[i] -= (this->_pl << 1);
					if (fft[i] >= this->_pl) fft[i] -= this->_pl;
				}
			}
		}

	}; // FFT_algorithms<Field, NoSimd<typename Field::Element>, 8>

	template<typename Field>
	class FFT_algorithms<Field, NoSimd<typename Field::Element>, 1, IS_FLOATING> : public FFT_butterflies<Field, NoSimd<typename Field::Element>, 1> {
	public:
		using Element = typename Field::Element;

		FFT_algorithms(const FFT_init<Field>& f_i) : FFT_butterflies<Field, NoSimd<typename Field::Element>, 1>(f_i) {}

		void DIF (Element *fft) {
			for (size_t w = this->n >> 1, f = 1; w != 0; f <<= 1, w >>= 1){
				// w : witdh of butterflies
				// f : # families of butterflies
				for (size_t i = 0; i < f; i++)
					for (size_t j = 0; j < w; j++)
						this->Butterfly_DIF(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], (this->pow_w)[j*f]);
			}
		}

		void DIT (Element *fft) {
			for (size_t w = 1, f = this->n >> 1; f >= 1; w <<= 1, f >>= 1)
				for (size_t i = 0; i < f; i++)
					for (size_t j = 0; j < w; j++)
						this->Butterfly_DIT(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], (this->pow_w)[j*f]);
		}

	}; // FFT_algorithms<Field, NoSimd<typename Field::Element>, 1>

	template<typename Field, typename simd>
	class FFT_algorithms<Field, simd, 4, IS_FLOATING> : public FFT_butterflies<Field, simd, 4> {
	public:
		using Element = typename Field::Element;
		using Compute_t = typename Field::Compute_t;
		using Residu_t = typename Field::Residu_t;
		using vect_t = typename simd::vect_t;

		FFT_algorithms(const FFT_init<Field>& f_i) : FFT_butterflies<Field, simd, 4>(f_i) {
			linbox_check(simd::vect_size == 4);
		}

		void DIF (Element *fft) {
			const uint64_t& n = this->n;
			const Compute_t& _pl = static_cast<Compute_t>(this->_pl);

			vect_t P;
			P  = simd::set1(_pl);
			Compute_t u = 1.0/_pl;
			vect_t U = simd::set1(u);
			Element * tab_w = &(this->pow_w) [0];
			size_t w, f;
			for (w = n >> 1, f = 1; w >= 4; tab_w+=w, w >>= 1, f <<= 1){
				// w : witdh of butterflies
				// f : # families of butterflies
				for (size_t i = 0; i < f; i++)
					for (size_t j = 0; j < w; j+=4)

#define A0 &fft[0] +  (i << 1)   *w+ j
#define A4 &fft[0] + ((i << 1)+1)*w+ j
						this->Butterfly_DIF(A0,A4, tab_w+j,P,U);
#undef A0
#undef A4
				//std::cout<<fft<<std::endl;
			}
			// Last two steps
			if (n >= 8) {
				vect_t W;
				W = simd::set1 (tab_w [1]);

				for (size_t i = 0; i < f; i+=2)
#define A0 &fft[0] +  (i << 2)
#define A4 &fft[0] + ((i << 2)+4)
					this->Butterfly_DIF_laststeps(A0,A4,W,P,U);
				//std::cout<<fft<<std::endl;
#undef A0
#undef A4
			} else {
				FFT_algorithms<Field, NoSimd<Element>, 1> fft_algo_1 (FFT_init<Field> (this->field(),this->ln,this->getRoot()));

				for (; w >= 1; tab_w+=w, w >>= 1, f <<= 1)
					for (size_t i = 0; i < f; i++)
						for (size_t j = 0; j < w; j++)
							fft_algo_1.Butterfly_DIF(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], (this->pow_w)[j*f]);
			}
		}

		void DIT (Element *fft) {
			const uint64_t& n = this->n;
			const Compute_t& _pl = static_cast<Compute_t>(this->_pl);

			vect_t P;
			P  = simd::set1(_pl);
			Compute_t u = 1.0/_pl;
			vect_t U = simd::set1(u);

			// First two steps
			if (n >= 8) {
				vect_t W;
				W = simd::set1 ((this->pow_w) [n-3]);

				for (size_t i = 0; i < n; i+=8)
					this->Butterfly_DIT_firststeps(&fft[i],&fft[i+4],W,P,U);

				Element * tab_w = &(this->pow_w) [n-8];
				for (size_t w = 4, f = n >> 3; f >= 1; w <<= 1, f >>= 1, tab_w-=w){
					// w : witdh of butterflies
					// f : # families of butterflies
					for (size_t i = 0; i < f; i++)
						for (size_t j = 0; j < w; j+=4)
#define A0 &fft[0] +  (i << 1)   *w+ j
#define A4 &fft[0] + ((i << 1)+1)*w+ j
							this->Butterfly_DIT(A0,A4, tab_w+j,P,U);

#undef A0
#undef A4

				}
			} else {
				FFT_algorithms<Field, NoSimd<Element>, 1> fft_algo_1 (FFT_init<Field> (this->field(),this->ln,this->getRoot()));

				Element * tab_w = &(this->pow_w) [n-2];
				for (size_t w = 1, f = n >> 1; f >= 1; w <<= 1, f >>= 1, tab_w-=w)
					for (size_t i = 0; i < f; i++)
						for (size_t j = 0; j < w; j++)
							fft_algo_1.Butterfly_DIT(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], tab_w[j]);
			}
		}

	}; // FFT_algorithms<Field, NoSimd<typename Field::Element>, 4>
}

#undef IS_INTEGRAL
#undef IS_FLOATING

#endif // __LINBOX_polynomial_fft_algorithms_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

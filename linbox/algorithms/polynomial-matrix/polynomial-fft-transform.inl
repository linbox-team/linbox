/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

namespace LinBox {

	template <class Field>
	template <class T>
	inline void FFT_transform<Field>::Butterfly_DIT_mod4p(T& A, T& B, const uint32_t& alpha, const uint32_t& alphap) {
		// Harvey's algorithm
		// 0 <= A,B < 4*p, p < 2^32 / 4
		// alphap = Floor(alpha * 2^ 32 / p])
		if (A >= _dpl) A -= _dpl;
		uint32_t tmp = ((uint32_t) alphap * (uint64_t)B) >> 32;
		tmp = (uint64_t)alpha * B - tmp * _pl;
		B = A + (_dpl - tmp);
		//        B &= 0XFFFFFFFF;
		A += tmp;
	}

	template <class Field>
	template <class T>
	inline void FFT_transform<Field>::Butterfly_DIF_mod2p(T& A, T& B, const uint32_t& alpha, const uint32_t& alphap) {
		//std::cout<<A<<" $$ "<<B<<"("<<alpha<<","<<alphap<<" ) -> ";
		// Harvey's algorithm
		// 0 <= A,B < 2*p, p < 2^32 / 4
		// alphap = Floor(alpha * 2^ 32 / p])
		uint64_t tmp = A;
		A += B;
		if (A >= _dpl) A -= _dpl;
		B = tmp + (_dpl - B);		
		tmp = ((uint32_t) alphap * (uint64_t)B) >> 32;
		B = (uint64_t)alpha * B - tmp * _pl;
		//B &= 0xFFFFFFFF;
		//std::cout<<A<<" $$ "<<B<<"\n ";
	}


	template <class Field>
	void FFT_transform<Field>::FFT_DIF_Harvey_mod2p_iterative (uint32_t *fft) {
		for (size_t w = n >> 1, f = 1, pos_w = 0; w != 0; f <<= 1, pos_w += w, w >>= 1){
			// w : witdh of butterflies
			// f : # families of butterflies
			for (size_t i = 0; i < f; i++)
				for (size_t j = 0; j < w; j++)
					Butterfly_DIF_mod2p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], pow_w[j*f], pow_wp[j*f]);
		}
	}

	template <class Field>
	void FFT_transform<Field>::FFT_DIF_Harvey_mod2p_iterative2x2 (uint32_t *fft) {
		size_t w, f;
		for (w = n >> 1, f = 1; w >= 2; w >>= 2, f <<= 2)
			// w : witdh of butterflies
			// f : # families of butterflies
			for (size_t i = 0; i < f; i++)
				for (size_t j = 0; j < (w >> 1); j++) {
#define A0 fft[(i << 1)    *w+ j          ]
#define A1 fft[(i << 1)    *w+(j+(w >> 1))]
#define A2 fft[((i << 1)+1)*w+ j          ]
#define A3 fft[((i << 1)+1)*w+(j+(w >> 1))]
					// Base butterfly
					Butterfly_DIF_mod2p(A0, A2, pow_w[j*f], pow_wp[j*f]);
					// Same step, same family, index + (w >>1)
					Butterfly_DIF_mod2p(A1, A3, pow_w[(j+(w >> 1))*f], pow_wp[(j+(w >> 1))*f]);
					// Next step on first entries of previous butterflies, #family * 2, same index
					Butterfly_DIF_mod2p(A0, A1, pow_w[j*(f << 1)], pow_wp[j*(f << 1)]);
					// Next step on second entries of previous butterflies, #family * 2, same index
					Butterfly_DIF_mod2p(A2, A3, pow_w[j*(f << 1)], pow_wp[j*(f << 1)]);
#undef A0
#undef A1
#undef A2
#undef A3
				}
		// Remaining steps, at most one
		if (w > 0)
			// w : witdh of butterflies
			// f : # families of butterflies
			for (size_t i = 0; i < f; i++)
				for (size_t j = 0; j < w; j++)
					Butterfly_DIF_mod2p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], pow_w[j*f], pow_wp[j*f]);
	}

	template <class Field>
	void FFT_transform<Field>::FFT_DIF_Harvey_mod2p_iterative3x3 (uint32_t *fft) {
		size_t w, f;
		for (w = n >> 1, f = 1; w >= 4; w >>= 3, f <<= 3)
			// w : witdh of butterflies
			// f : # families of butterflies
			for (size_t i = 0; i < f; i++)
				for (size_t j = 0; j < (w >> 2); j++) {
#define A0 fft[(i << 1)    *w+ j            ]
#define A1 fft[(i << 1)    *w+(j+  (w >> 2))]
#define A2 fft[(i << 1)    *w+(j+  (w >> 1))]
#define A3 fft[(i << 1)    *w+(j+3*(w >> 2))]
#define A4 fft[((i << 1)+1)*w+ j            ]
#define A5 fft[((i << 1)+1)*w+(j+  (w >> 2))]
#define A6 fft[((i << 1)+1)*w+(j+  (w >> 1))]
#define A7 fft[((i << 1)+1)*w+(j+3*(w >> 2))]
					// First Step
					Butterfly_DIF_mod2p(A0, A4, pow_w[ j            *f], pow_wp[ j            *f]);
					Butterfly_DIF_mod2p(A1, A5, pow_w[(j+  (w >> 2))*f], pow_wp[(j+  (w >> 2))*f]);
					Butterfly_DIF_mod2p(A2, A6, pow_w[(j+  (w >> 1))*f], pow_wp[(j+  (w >> 1))*f]);
					Butterfly_DIF_mod2p(A3, A7, pow_w[(j+3*(w >> 2))*f], pow_wp[(j+3*(w >> 2))*f]);

					// Second Step
					Butterfly_DIF_mod2p(A0, A2, pow_w[ j          *(f << 1)], pow_wp[ j          *(f << 1)]);
					Butterfly_DIF_mod2p(A1, A3, pow_w[(j+(w >> 2))*(f << 1)], pow_wp[(j+(w >> 2))*(f << 1)]);
					Butterfly_DIF_mod2p(A4, A6, pow_w[ j          *(f << 1)], pow_wp[ j          *(f << 1)]);
					Butterfly_DIF_mod2p(A5, A7, pow_w[(j+(w >> 2))*(f << 1)], pow_wp[(j+(w >> 2))*(f << 1)]);

					// Third Step
					Butterfly_DIF_mod2p(A0, A1, pow_w[j*(f << 2)], pow_wp[j*(f << 2)]);
					Butterfly_DIF_mod2p(A2, A3, pow_w[j*(f << 2)], pow_wp[j*(f << 2)]);
					Butterfly_DIF_mod2p(A4, A5, pow_w[j*(f << 2)], pow_wp[j*(f << 2)]);
					Butterfly_DIF_mod2p(A6, A7, pow_w[j*(f << 2)], pow_wp[j*(f << 2)]);
#undef A0
#undef A1
#undef A2
#undef A3
#undef A4
#undef A5
#undef A6
#undef A7
				}
		// Remaining steps, at most two
		for (; w > 0; w >>= 1, f <<= 1)
			// w : witdh of butterflies
			// f : # families of butterflies
			for (size_t i = 0; i < f; i++)
				for (size_t j = 0; j < w; j++)
					Butterfly_DIF_mod2p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], pow_w[j*f], pow_wp[j*f]);
	}

	template <class Field>
	void FFT_transform<Field>::FFT_DIT_Harvey_mod4p_iterative (uint32_t *fft) {
		for (size_t w = 1, f = n >> 1; f >= 1; w <<= 1, f >>= 1)
			for (size_t i = 0; i < f; i++)
				for (size_t j = 0; j < w; j++)
					Butterfly_DIT_mod4p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], pow_w[j*f], pow_wp[j*f]);
	}

	template <class Field>
	void FFT_transform<Field>::FFT_DIT_Harvey_mod4p_iterative2x2 (uint32_t *fft) {
		size_t w, f;
		for (w = 1, f = n >> 1; f >= 2; w <<= 2, f >>= 2)
			// w : witdh of butterflies
			// f : # families of butterflies
			for (size_t i = 0; i < (f >> 1); i++)
				for (size_t j = 0; j < w; j++) {
#define A0 fft[ (i << 2)   *w+j]
#define A1 fft[((i << 2)+1)*w+j]
#define A2 fft[((i << 2)+2)*w+j]
#define A3 fft[((i << 2)+3)*w+j]
					Butterfly_DIT_mod4p(A0, A1, pow_w[j*f], pow_wp[j*f]);
					Butterfly_DIT_mod4p(A2, A3, pow_w[j*f], pow_wp[j*f]);

					Butterfly_DIT_mod4p(A0, A2, pow_w[ j   *(f >> 1)], pow_wp[ j   *(f >> 1)]);
					Butterfly_DIT_mod4p(A1, A3, pow_w[(j+w)*(f >> 1)], pow_wp[(j+w)*(f >> 1)]);
#undef A0
#undef A1
#undef A2
#undef A3
				}
		// Remaining steps, at most one
		if (f > 0)
			// w : witdh of butterflies
			// f : # families of butterflies
			for (size_t i = 0; i < f; i++)
				for (size_t j = 0; j < w; j++)
					Butterfly_DIT_mod4p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], pow_w[j*f], pow_wp[j*f]);
	}

	template <class Field>
	void FFT_transform<Field>::FFT_DIT_Harvey_mod4p_iterative3x3 (uint32_t *fft) {
		size_t w, f;
		for (w = 1, f = n >> 1; f >= 4; w <<= 3, f >>= 3)
			// w : witdh of butterflies
			// f : # families of butterflies
			for (size_t i = 0; i < (f >> 2); i++)
				for (size_t j = 0; j < w; j++) {
#define A0 fft[ (i << 3)   *w+j]
#define A1 fft[((i << 3)+1)*w+j]
#define A2 fft[((i << 3)+2)*w+j]
#define A3 fft[((i << 3)+3)*w+j]
#define A4 fft[((i << 3)+4)*w+j]
#define A5 fft[((i << 3)+5)*w+j]
#define A6 fft[((i << 3)+6)*w+j]
#define A7 fft[((i << 3)+7)*w+j]
					Butterfly_DIT_mod4p(A0, A1, pow_w[j*f], pow_wp[j*f]);
					Butterfly_DIT_mod4p(A2, A3, pow_w[j*f], pow_wp[j*f]);
					Butterfly_DIT_mod4p(A4, A5, pow_w[j*f], pow_wp[j*f]);
					Butterfly_DIT_mod4p(A6, A7, pow_w[j*f], pow_wp[j*f]);

					Butterfly_DIT_mod4p(A0, A2, pow_w[ j   *(f >> 1)], pow_wp[ j   *(f >> 1)]);
					Butterfly_DIT_mod4p(A1, A3, pow_w[(j+w)*(f >> 1)], pow_wp[(j+w)*(f >> 1)]);
					Butterfly_DIT_mod4p(A4, A6, pow_w[ j   *(f >> 1)], pow_wp[ j   *(f >> 1)]);
					Butterfly_DIT_mod4p(A5, A7, pow_w[(j+w)*(f >> 1)], pow_wp[(j+w)*(f >> 1)]);

					Butterfly_DIT_mod4p(A0, A4, pow_w[ j     *(f >> 2)], pow_wp[ j     *(f >> 2)]);
					Butterfly_DIT_mod4p(A1, A5, pow_w[(j+  w)*(f >> 2)], pow_wp[(j+  w)*(f >> 2)]);
					Butterfly_DIT_mod4p(A2, A6, pow_w[(j+2*w)*(f >> 2)], pow_wp[(j+2*w)*(f >> 2)]);
					Butterfly_DIT_mod4p(A3, A7, pow_w[(j+3*w)*(f >> 2)], pow_wp[(j+3*w)*(f >> 2)]);
#undef A0
#undef A1
#undef A2
#undef A3
#undef A4
#undef A5
#undef A6
#undef A7
				}
		// Remaining steps, at most one
		for ( ; f > 0; w <<= 1, f >>= 1)
			// w : witdh of butterflies
			// f : # families of butterflies
			for (size_t i = 0; i < f; i++)
				for (size_t j = 0; j < w; j++)
					Butterfly_DIT_mod4p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], pow_w[j*f], pow_wp[j*f]);
	}


}

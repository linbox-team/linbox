/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/vector/stream.h
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2003-02-03 Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * RandomSparseStream::RandomSparseStream: put probability parameter before
 * vector dimension
 * ------------------------------------
 * 2002-09-05 Bradford Hovinen <bghovine@math.uwaterloo.ca>
 *
 *  - Renamed to stream.h and moved to linbox/vector
 *  - VectorFactory is now called VectorStream, which fits its purpose
 *    somewhat better
 *  - The interface is now closer to the interface for istream
 *  - RandomDenseVectorFactory, et al. are refactored into classes
 *    parameterized on the vector type and specialized appropriately. This
 *    allows, e.g. construction of a random dense vector in sparse
 *    representation and so on.
 *  - New constructor interface for RandomSparseStream accepts proportion of
 *    nonzero entries; removed existing constructors
 *  - Reindented, since the other changes are so numerous that diffs aren't a
 *    big deal anyway
 *
 * ------------------------------------
 * 2002-05-18 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Refactor: Create one class StandardBasisFactory, parameterized by vector
 * type, with specializations for dense, sparse map, and sparse associative
 * vectors.
 *
 * ------------------------------------
 * Modified by Dmitriy Morozov <linbox@foxcub.org>. May 27, 2002.
 *
 * Added parametrization of the VectorCategroy tags by VectorTraits (see
 * vector-traits.h for more details).
 *
 * ------------------------------------
 *
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
 *.
 */

/**
 * @file vector/stream.h
 *
 * @brief Generation of sequences of random vectors.
 *
 * Random, sparse, basis vectors,...
 */

#ifndef __LINBOX_vector_stream_gf2_H
#define __LINBOX_vector_stream_gf2_H

#include <vector>
#include <cmath>

#include "linbox/field/gf2.h"
#include "linbox/vector/stream.h"

namespace LinBox
{ /*  RandomDenseStream */
	class RandomDenseStreamGF2 : public VectorStream<BitVector> {
	public:
		typedef BitVector Vector;

		RandomDenseStreamGF2 (const GF2 &, uint32_t seed, size_t nn, size_t mm = 0) :
			MT (seed), _n (nn), _m (mm), _j (0)
		{}

		Vector &get (Vector &v)
		{
			Vector::word_iterator i;

			if (_m > 0 && _j++ >= _m)
				return v;

			for (i = v.wordBegin (); i != v.wordEnd (); i++)
				*i = MTrandomInt<__LINBOX_BITSOF_LONG>()(MT);

			const size_t zeroing = __LINBOX_BITSOF_LONG - (v.size() % __LINBOX_BITSOF_LONG);
			*(v.wordRbegin()) <<= zeroing;
			*(v.wordRbegin()) >>= zeroing;
			return v;
		}

		size_t size () const
		{ return _m; }
		size_t pos () const
		{ return _j; }
		size_t dim () const
		{ return _n; }
		operator bool () const
		{ return _m == 0 || _j < _m; }
		void reset () { _j = 0; }

	private:
		MersenneTwister MT;
		size_t          _n;
		size_t          _m;
		size_t          _j;
	};
}

// Specialization of RandomSparseStream
namespace LinBox
{ /*  RandomSparseStream */

	template <class _Vector = Vector<GF2>::Sparse>
	class RandomSparseStreamGF2 : public VectorStream<_Vector> {
	public:
		typedef GF2 Field;
		typedef _Vector Vector;

		RandomSparseStreamGF2 (const GF2 &, uint32_t seed, double p, size_t N, size_t M = 0) :
			MT (seed), _n (N), _m (M), _j (0)
		{ setP (p); }

		RandomSparseStreamGF2 (const GF2 &F, const GF2RandIter& r, double p, size_t N, size_t M = 0) :
			MT (r.getMT()), _n (N), _m (M), _j (0)
		{ setP (p); }

		Vector &get (Vector &v);

		size_t size () const
		{ return _m; }
		size_t pos () const
		{ return _j; }
		size_t dim () const
		{ return _n; }
		operator bool () const
		{ return _m == 0 || _j < _m; }
		void reset () { _j = 0; }

		void setP (double p)
		{
			linbox_check ((p >= 0.0) && (p <= 1.0));
			_p = p;
			_1_log_1mp   = 1 / log (1 - _p);
		}

	private:
		MersenneTwister MT;
		size_t _n;
		double _p;
		double _1_log_1mp;
		size_t _m;
		size_t _j;
	};

	template <class _Vector>
	_Vector &RandomSparseStreamGF2<_Vector>::get (_Vector &v)
	{
		size_t i = (size_t) -1;
		double val;
		int skip;

		if (_m > 0 && _j++ >= _m)
			return v;

		v.clear ();

		while (1) {
			val = (double) MT.randomDouble ();
			skip = (int) (ceil (log (val) * _1_log_1mp));

			if (skip <= 0)
				i++;
			else
				i += (size_t) skip;

			if (i >= _n) break;

			v.push_back (i);
		}

		return v;
	}

}


#endif // __LINBOX_vector_stream_H



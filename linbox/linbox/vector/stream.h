/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/util/vector-factory.h
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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
 * See COPYING for license information.
 */

#ifndef __UTIL_VECTOR_FACTORY_H
#define __UTIL_VECTOR_FACTORY_H

#include <vector>
#include <cmath>

#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/randiter/nonzero.h"

namespace LinBox 
{
	/** Vector factory
	 * This is an abstract base class that generates a sequence of vectors
	 * in a generic way. Typical uses would be in tests, where the same test
	 * might be run on a sequence of random vectors or on e_1, ..., e_n.
	 */
	template <class Vector>
	class VectorFactory 
	{
	    public:
		virtual ~VectorFactory () {}
		virtual Vector &next (Vector &v) = 0;
		virtual size_t j () const = 0;
		virtual size_t m () const = 0;
		virtual size_t n () const = 0;
		virtual operator bool () const = 0;
		virtual void reset () = 0;
	};

	/** Constant vector factory
	 * Returns the same vector repeatedly
	 */
	template <class Vector>
	class ConstantVectorFactory : public VectorFactory<Vector>
	{
	    public:

		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param v Vector to return on next
		 * @param m Number of vectors to return (0 for unlimited)
		 */
		ConstantVectorFactory (Vector &v, size_t m) : _v (v), _m (m), _j (0) {}

		/** Retrieve vector
		 * @param v Vector to use
		 */
		Vector &next (Vector &v) 
			{ if (_m == 0 || _j < _m) copy (_v.begin (), _v.end (), v.begin ()); return v; }

		/** Number of vectors created so far
		 */
		size_t j () const { return _j; }

		/** Number of vectors to be created
		 */
		size_t m () const { return _m; }

		/** Dimension of the space
		 */
		size_t n () const { return _v.size (); }

		/** Check whether we have reached the end
		 */
		operator bool () const { return _m == 0 || _j < _m; }

		/** Reset the factory to start at the beginning
		 */
		void reset () { _j = 0; }

	    private:
		Vector &_v;
		size_t  _m;
		size_t  _j;
	};

	/** Random dense vector factory
	 * Generates a sequence of random dense vectors over a given field
	 */
	template <class Field, class RandIter = typename Field::RandIter>
	class RandomDenseVectorFactory : public VectorFactory<std::vector<typename Field::Element> >
	{
	    public:
		typedef std::vector<typename Field::Element> Vector;

		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param F Field over which to create random vectors
		 * @param n Size of vectors
		 * @param m Number of vectors to return (0 for unlimited)
		 */
		RandomDenseVectorFactory (const Field &F, size_t n, size_t m = 0)
			: _F (F), _r (F), _n (n), _m (m), _j (0)
			{}

		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param F Field over which to create random vectors
		 * @param n Size of vectors
		 * @param m Number of vectors to return (0 for unlimited)
		 */
		RandomDenseVectorFactory (const Field &F, const RandIter &r, size_t n, size_t m = 0)
			: _F (F), _r (r), _n (n), _m (m), _j (0)
			{}

		/** Get next element
		 * @param v Vector into which to generate random vector
		 * @return reference to new random vector
		 */
		Vector &next (Vector &v) 
		{
			typename Vector::iterator i;

			v.resize (_n);

			if (_m > 0 && _j++ >= _m)
				return v;

			for (i = v.begin (); i < v.end (); i++)
				_r.random (*i);

			return v;
		}

		/** Number of vectors created so far
		 */
		size_t j () const { return _j; }

		/** Number of vectors to be created
		 */
		size_t m () const { return _m; }

		/** Dimension of the space
		 */
		size_t n () const { return _n; }

		/** Check whether we have reached the end
		 */
		operator bool () const 
			{ return _m == 0 || _j < _m; }

		/** Reset the factory to start at the beginning
		 */
		void reset () { _j = 0; }

	    private:
		const Field &_F;
		RandIter     _r;
		size_t       _n;
		size_t       _m;
		size_t       _j;
	};

	/** Random sparse vector factory
	 * Generates a sequence of random sparse vectors over a given field
	 */
	template <class Field, class RandIter = typename Field::RandIter>
	class RandomSparseSeqVectorFactory : public VectorFactory<std::vector<std::pair<size_t, typename Field::Element> > >
	{
	    public:
		typedef std::vector<std::pair<size_t, typename Field::Element> > Vector;

		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param F Field over which to create random vectors
		 * @param n Size of vectors
		 * @param k Expected number of nonzero entries
		 * @param m Number of vectors to return (0 for unlimited)
		 */
		RandomSparseSeqVectorFactory (const Field &F, size_t n, size_t k, size_t m = 0)
			: _F (F), _r (F, typename Field::RandIter (F)), _n (n), _k (k), _m (m), _j (0)
		{
			linbox_check (k < n);

			_p           = (double) _k / (double) _n;
			_1_log_1mp   = 1 / log (1 - _p);
		}

		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param F Field over which to create random vectors
		 * @param n Size of vectors
		 * @param k Expected number of nonzero entries
		 * @param m Number of vectors to return (0 for unlimited)
		 */
		RandomSparseSeqVectorFactory (const Field &F, const RandIter &r, size_t n, size_t k, size_t m = 0)
			: _F (F), _r (F, r), _n (n), _k (k), _m (m), _j (0)
		{
			linbox_check (k < n);

			_p           = (double) _k / (double) _n;
			_1_log_1mp   = 1 / log (1 - _p);
		}

		/** Get next element
		 * @param v Vector into which to generate random vector
		 * @return reference to new random vector
		 */
		Vector &next (Vector &v) 
		{
			typename Field::Element x;
			size_t i = 0;
			double val;
			int skip;

			if (_m > 0 && _j++ >= _m)
				return v;

			v.clear ();

			while (1) {
				val = (double) ((unsigned long) rand ()) / (0.5 * (double) ((unsigned long) -1));
				skip = (int) (ceil (log (val) * _1_log_1mp));

				if (skip <= 0)
					i++;
				else
					i += skip;

				if (i >= _n) break;

				_r.random (x);
				v.push_back (std::pair<size_t, typename Field::Element> (i, x));
			}

			return v;
		}

		/** Number of vectors created so far
		 */
		size_t j () const { return _j; }

		/** Number of vectors to be created
		 */
		size_t m () const { return _m; }

		/** Dimension of the space
		 */
		size_t n () const { return _n; }

		/** Check whether we have reached the end
		 */
		operator bool () const 
			{ return _m == 0 || _j < _m; }

		/** Reset the factory to start at the beginning
		 */
		void reset () { _j = 0; }

	    private:
		const Field                      &_F;
		NonzeroRandIter<Field, RandIter>  _r;
		size_t                            _n;
		long                              _k;
		double                            _p;
		double                            _1_log_1mp;
		size_t                            _m;
		size_t                            _j;
	};

	/** Random sparse vector factory
	 * Generates a sequence of random sparse vectors over a given field
	 */
	template <class Field, class RandIter = typename Field::RandIter>
	class RandomSparseMapVectorFactory : public VectorFactory<std::map<size_t, typename Field::Element> >
	{
	    public:
		typedef std::map<size_t, typename Field::Element> Vector;

		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param F Field over which to create random vectors
		 * @param n Size of vectors
		 * @param k Expected number of nonzero entries
		 * @param m Number of vectors to return (0 for unlimited)
		 */
		RandomSparseMapVectorFactory (const Field &F, size_t n, size_t k, size_t m = 0)
			: _F (F), _r (F, typename Field::RandIter (F)), _n (n), _k (k), _j (0), _m (m)
		{}

		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param F Field over which to create random vectors
		 * @param n Size of vectors
		 * @param k Expected number of nonzero entries
		 * @param m Number of vectors to return (0 for unlimited)
		 */
		RandomSparseMapVectorFactory (const Field &F, const RandIter &r, size_t n, size_t k, size_t m = 0)
			: _F (F), _r (F, r), _n (n), _k (k), _j (0), _m (m)
		{}

		/** Get next element
		 * @param v Vector into which to generate random vector
		 * @return reference to new random vector
		 */
		Vector &next (Vector &v) 
		{
			typename Field::Element x;
			int i, idx;

			if (_m > 0 && _j++ >= _m)
				return v;

			v.clear ();

			for (i = 0; i < _k; i++) {
				_r.random (x);
				while (!_F.isZero (v[(idx = rand () % _n)]));
				v[idx] = x;
			}

			return v;
		}

		/** Number of vectors created so far
		 */
		size_t j () const { return _j; }

		/** Number of vectors to be created
		 */
		size_t m () const { return _m; }

		/** Dimension of the space
		 */
		size_t n () const { return _n; }

		/** Check whether we have reached the end
		 */
		operator bool () const 
			{ return _m == 0 || _j < _m; }

		/** Reset the factory to start at the beginning
		 */
		void reset () { _j = 0; }

	    private:
		const Field                      &_F;
		NonzeroRandIter<Field, RandIter>  _r;
		size_t                            _n;
		long                              _k;
		size_t                            _j;
		size_t                            _m;
	};

	/** Random sparse vector factory
	 * Generates a sequence of random sparse vectors over a given field
	 */
	template <class Field, class RandIter = typename Field::RandIter>
	class RandomSparseParVectorFactory : public VectorFactory<std::pair<std::vector<size_t>, std::vector<typename Field::Element> > >
	{
	    public:
		typedef std::pair<std::vector<size_t>, std::vector<typename Field::Element> > Vector;

		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param F Field over which to create random vectors
		 * @param n Size of vectors
		 * @param k Expected number of nonzero entries
		 * @param m Number of vectors to return (0 for unlimited)
		 */
		RandomSparseParVectorFactory (const Field &F, size_t n, size_t k, size_t m = 0)
			: _F (F), _r (F, typename Field::RandIter (F)), _n (n), _k (k), _m (m), _j (0)
		{
			linbox_check (k < n);

			_p           = (double) _k / (double) _n;
			_1_log_1mp   = 1 / log (1 - _p);
		}

		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param F Field over which to create random vectors
		 * @param n Size of vectors
		 * @param k Expected number of nonzero entries
		 * @param m Number of vectors to return (0 for unlimited)
		 */
		RandomSparseParVectorFactory (const Field &F, const RandIter &r, size_t n, size_t k, size_t m = 0)
			: _F (F), _r (F, r), _n (n), _k (k), _m (m), _j (0)
		{
			linbox_check (k < n);

			_p           = (double) _k / (double) _n;
			_1_log_1mp   = 1 / log (1 - _p);
		}

		/** Get next element
		 * @param v Vector into which to generate random vector
		 * @return reference to new random vector
		 */
		Vector &next (Vector &v) 
		{
			typename Field::Element x;
			size_t i = 0;
			double val;
			int skip;

			if (_m > 0 && _j++ >= _m)
				return v;

			v.first.clear ();
			v.second.clear ();

			while (1) {
				val = (double) ((unsigned long) rand ()) / (0.5 * (double) ((unsigned long) -1));
				skip = (int) (ceil (log (val) * _1_log_1mp));

				if (skip <= 0)
					i++;
				else
					i += skip;

				if (i >= _n) break;

				_r.random (x);
				v.first.push_back (i);
				v.second.push_back (x);
			}

			return v;
		}

		/** Number of vectors created so far
		 */
		size_t j () const { return _j; }

		/** Number of vectors to be created
		 */
		size_t m () const { return _m; }

		/** Dimension of the space
		 */
		size_t n () const { return _n; }

		/** Check whether we have reached the end
		 */
		operator bool () const 
			{ return _m == 0 || _j < _m; }

		/** Reset the factory to start at the beginning
		 */
		void reset () { _j = 0; }

	    private:
		const Field                      &_F;
		NonzeroRandIter<Field, RandIter>  _r;
		size_t                            _n;
		long                              _k;
		double                            _p;
		double                            _1_log_1mp;
		size_t                            _m;
		size_t                            _j;
	};

	/** Factory for e_1,...,e_n
	 * Generates the sequence e_1,...,e_n over a given field
	 * 
	 * This class is generic with respect to the underlying vector
	 * representation.
	 */
	template <class Field, class Vector, class Trait = typename VectorTraits<Vector>::VectorCategory>
	class StandardBasisFactory : public VectorFactory<Vector>
	{
	    public:

		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param F Field over which to create vectors
		 * @param n Size of vectors
		 */
		StandardBasisFactory (Field &F, size_t n);

		/** Get next element
		 * @param v Vector into which to generate vector
		 * @return reference to new vector
		 */
		Vector &next (Vector &v);

		/** Number of vectors created so far
		 */
		size_t j () const;

		/** Number of vectors to be created
		 */
		size_t m () const;

		/** Dimension of the space
		 */
		size_t n () const;

		/** Check whether we have reached the end
		 */
		operator bool () const;

		/** Reset the factory to start at the beginning
		 */
		void reset ();

	    private:
		const Field              &_F;
		size_t                    _n;
		size_t                    _j;
	};

	/* Specialization of standard basis factory for dense vectors */

	template <class Field, class Vector, class VectorTrait>
	class StandardBasisFactory<Field, Vector, VectorCategories::DenseVectorTag<VectorTrait> > : public VectorFactory<Vector>
	{
	    public:
		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param F Field over which to create vectors
		 * @param n Size of vectors
		 */
		StandardBasisFactory (const Field &F, size_t n)
			: _F (F), _n (n), _j (0)
			{}

		/** Get next element
		 * @param v Vector into which to generate vector
		 * @return reference to new vector
		 */
		Vector &next (Vector &v) 
		{
			static typename Field::Element zero;
			typename Vector::iterator i;
			size_t idx;

			for (i = v.begin (), idx = 0; i != v.end (); i++, idx++) {
				if (idx == _j)
					_F.init (*i, 1);
				else
					_F.assign (*i, zero);
			}

			_j++;

			return v;
		}

		/** Number of vectors created so far
		 */
		size_t j () const { return _j; }

		/** Number of vectors to be created
		 */
		size_t m () const { return _n; }

		/** Dimension of the space
		 */
		size_t n () const { return _n; }

		/** Check whether we have reached the end
		 */
		operator bool () const 
			{ return _j < _n; }

		/** Reset the factory to start at the beginning
		 */
		void reset () { _j = 0; }

	    private:
		const Field              &_F;
		size_t                    _n;
		size_t                    _j;
	};

	/* Specialization of standard basis factory for sparse sequence vectors */

	template <class Field, class Vector, class VectorTrait>
	class StandardBasisFactory<Field, Vector, VectorCategories::SparseSequenceVectorTag<VectorTrait> > : public VectorFactory<Vector>
	{
	    public:
		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param F Field over which to create vectors
		 * @param n Size of vectors
		 */
		StandardBasisFactory (Field &F, size_t n)
			: _F (F), _n (n), _j (0)
			{ _F.init (_one, 1); }

		/** Get next element
		 * @param v Vector into which to generate vector
		 * @return reference to new vector
		 */
		Vector &next (Vector &v) 
		{
			v.clear ();

			if (_j < _n)
				v.push_back (std::pair <size_t, typename Field::Element> (_j++, _one));

			return v;
		}

		/** Number of vectors created so far
		 */
		size_t j () const { return _j; }

		/** Number of vectors to be created
		 */
		size_t m () const { return _n; }

		/** Dimension of the space
		 */
		size_t n () const { return _n; }

		/** Check whether we have reached the end
		 */
		operator bool () const 
			{ return _j < _n; }

		/** Reset the factory to start at the beginning
		 */
		void reset () { _j = 0; }

	    private:
		const Field              &_F;
		size_t                    _n;
		size_t                    _j;
		typename Field::Element   _one;
	};

	/** Factory for e_1,...,e_n in the sparse associative vector representation
	 * Generates the sequence e_1,...,e_n in the sparse associative vector representation over a given field
	 */
	template <class Field, class Vector, class VectorTrait>
	class StandardBasisFactory<Field, Vector, VectorCategories::SparseAssociativeVectorTag<VectorTrait> > : public VectorFactory<Vector>
	{
	    public:
		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param F Field over which to create vectors
		 * @param n Size of vectors
		 */
		StandardBasisFactory (Field &F, size_t n)
			: _F (F), _n (n), _j (0)
			{ _F.init (_one, 1); }

		/** Get next element
		 * @param v Vector into which to generate vector
		 * @return reference to new vector
		 */
		Vector &next (Vector &v) 
		{
			v.clear ();

			if (_j < _n)
				v.insert (std::pair <size_t, typename Field::Element> (_j++, _one));

			return v;
		}

		/** Number of vectors created so far
		 */
		size_t j () const { return _j; }

		/** Number of vectors to be created
		 */
		size_t m () const { return _n; }

		/** Dimension of the space
		 */
		size_t n () const { return _n; }

		/** Check whether we have reached the end
		 */
		operator bool () const 
			{ return _j < _n; }

		/** Reset the factory to start at the beginning
		 */
		void reset () { _j = 0; }

	    private:
		const Field              &_F;
		size_t                    _n;
		size_t                    _j;
		typename Field::Element   _one;
	};

	/* Specialization of standard basis factory for sparse parallel vectors */

	template <class Field, class Vector, class VectorTrait>
	class StandardBasisFactory<Field, Vector, VectorCategories::SparseParallelVectorTag<VectorTrait> > : public VectorFactory<Vector>
	{
	    public:
		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param F Field over which to create vectors
		 * @param n Size of vectors
		 */
		StandardBasisFactory (Field &F, size_t n)
			: _F (F), _n (n), _j (0)
			{ _F.init (_one, 1); }

		/** Get next element
		 * @param v Vector into which to generate vector
		 * @return reference to new vector
		 */
		Vector &next (Vector &v) 
		{
			v.first.clear ();
			v.second.clear ();

			if (_j < _n) {
				v.first.push_back (_j++);
				v.second.push_back (_one);
			}

			return v;
		}

		/** Number of vectors created so far
		 */
		size_t j () const { return _j; }

		/** Number of vectors to be created
		 */
		size_t m () const { return _n; }

		/** Dimension of the space
		 */
		size_t n () const { return _n; }

		/** Check whether we have reached the end
		 */
		operator bool () const 
			{ return _j < _n; }

		/** Reset the factory to start at the beginning
		 */
		void reset () { _j = 0; }

	    private:
		const Field              &_F;
		size_t                    _n;
		size_t                    _j;
		typename Field::Element   _one;
	};
}

#endif // __UTIL_VECTOR_FACTORY_H

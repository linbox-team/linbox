/* -*- mode: c; style: linux -*- */

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
		ConstantVectorFactory (Vector &v) : _v (v), _m (m), _j (0) {}

		/** Retrieve vector
		 * @param v Vector to use
		 */
		Vector &next (Vector &v) 
			{ if (_m == 0 || _j < _m) copy (_v.begin (), _v.end (), v.begin ()); return v; }

		/** Number of vectors created so far
		 */
		size_t j () const { return _j - 1; }

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
	template <class Field>
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

		/** Get next element
		 * @param v Vector into which to generate random vector
		 * @return reference to new random vector
		 */
		Vector &next (Vector &v) 
		{
			Vector::iterator i;

			v.resize (_n);

			if (_m > 0 && _j++ >= _m)
				return v;

			for (i = v.begin (); i < v.end (); i++)
				_r.random (*i);

			return v;
		}

		/** Number of vectors created so far
		 */
		size_t j () const { return _j - 1; }

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
		const Field              &_F;
		typename Field::RandIter  _r;
		size_t                    _n;
		size_t                    _m;
		size_t                    _j;
	};

	/** Random sparse vector factory
	 * Generates a sequence of random sparse vectors over a given field
	 */
	template <class Field>
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
			: _F (F), _r (F), _n (n), _k (k), _m (m), _j (0)
		{
			linbox_check (k < n);

			_p           = (double) _k / (double) _n;
			_log_1mp     = log (1 - _p);
			_ppm1        = _p * (_p - 1);
			_pm1_log_1mp = (_p - 1) * log (1 - _p);
		}

		/** Get next element
		 * @param v Vector into which to generate random vector
		 * @return reference to new random vector
		 */
		Vector &next (Vector &v) 
		{
			typename Field::Element x;
			int i = 0;
			double val;
			int skip;

			if (_m > 0 && _j++ >= _m)
				return v;

			v.clear ();

			while (1) {
				val = (double) ((unsigned long) rand ()) / (0.5 * (double) ((unsigned long) -1));
				skip = 2 + (int) floor (log ((val * _pm1_log_1mp - _p) / _ppm1) / _log_1mp);
				i += skip;
				if (i >= _n) break;

				_r.random (x);
				v.push_back (std::pair<size_t, typename Field::Element> (i, x));

			}

			return v;
		}

		/** Number of vectors created so far
		 */
		size_t j () const { return _j - 1; }

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
		const Field              &_F;
		typename Field::RandIter  _r;
		size_t                    _n;
		long                      _k;
		double                    _p;
		double                    _log_1mp;
		double                    _ppm1;
		double                    _pm1_log_1mp;
		size_t                    _m;
		size_t                    _j;
	};

	/** Random sparse vector factory
	 * Generates a sequence of random sparse vectors over a given field
	 */
	template <class Field>
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
			: _F (F), _r (F), _n (n), _k (k), _m (m), _j (0)
		{
			_p           = (double) _k / (double) _n;
			_log_1mp     = log (1 - _p);
			_ppm1        = _p * (_p - 1);
			_pm1_log_1mp = (_p - 1) * log (1 - _p);
		}

		/** Get next element
		 * @param v Vector into which to generate random vector
		 * @return reference to new random vector
		 */
		Vector &next (Vector &v) 
		{
			typename Field::Element x;
			int i = 0;
			double val;
			int skip;

			if (_m > 0 && _j++ >= _m)
				return v;

			v.clear ();

			while (1) {
				val = (double) ((unsigned long) rand ()) / (0.5 * (double) ((unsigned long) -1));
				skip = 2 + (int) floor (log ((val * _pm1_log_1mp - _p) / _ppm1) / _log_1mp);
				i += skip;
				if (i >= _n) break;

				_r.random (x);
				v.insert (std::pair <size_t, typename Field::Element> (i, x));
			}

			return v;
		}

		/** Number of vectors created so far
		 */
		size_t j () const { return _j - 1; }

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
		const Field              &_F;
		typename Field::RandIter  _r;
		size_t                    _n;
		long                      _k;
		double                    _p;
		double                    _log_1mp;
		double                    _ppm1;
		double                    _pm1_log_1mp;
		double                    _q;
		size_t                    _j;
		size_t                    _m;
	};

	/** Factory for e_1,...,e_n
	 * Generates the sequence e_1,...,e_n over a given field
	 * 
	 * This class is generic with respect to the underlying vector
	 * representation.
	 */
	template <class Field, class Vector, class Trait = VectorTraits<Vector>::VectorCategory>
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
			v.resize (_n);

			if (_j > 0 && _j < _n)
				_F.init (v[_j - 1], 0);

			if (_j < _n) {
				_F.init (v[_j], 1);
				_j++;
			}

			return v;
		}

		/** Number of vectors created so far
		 */
		size_t j () const { return _j - 1; }

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
		size_t j () const { return _j - 1; }

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
		size_t j () const { return _j - 1; }

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

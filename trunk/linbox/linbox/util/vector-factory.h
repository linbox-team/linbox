/* -*- mode: c; style: linux -*- */

/* linbox/util/vector-factory.h
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __UTIL_VECTOR_FACTORY_H
#define __UTIL_VECTOR_FACTORY_H

#include <vector>
#include <cmath>

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
		virtual operator bool () const = 0;
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
		size_t j () const { return _j; }

		/** Number of vectors to be created
		 */
		size_t m () const { return _m; }

		/** Check whether we have reached the end
		 */
		operator bool () const { return _m == 0 || _j < _m; }

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
		RandomDenseVectorFactory (Field &F, size_t n, size_t m = 0)
			: _F (F), _r (F), _n (n), _m (m), _j (0)
			{}

		/** Get next element
		 * @param v Vector into which to generate random vector
		 * @return reference to new random vector
		 */
		Vector &next (Vector &v) 
		{
			Vector::iterator i;

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

		/** Check whether we have reached the end
		 */
		operator bool () const 
			{ return _m == 0 || _j < _m; }

	    private:
		Field                    &_F;
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
		RandomSparseSeqVectorFactory (Field &F, size_t n, size_t k, size_t m = 0)
			: _F (F), _r (F), _n (n), _k (k), _m (m), _j (0)
		{
			_p       = (double) _k / (double) _n;
			_log_p   = log (_p);
			_log_1mp = log (1 - _p);
			_q       = log (_log_1mp) - _log_p;
		}

		/** Get next element
		 * @param v Vector into which to generate random vector
		 * @return reference to new random vector
		 */
		Vector &next (Vector &v) 
		{
			typename Field::Element x;
			int i = 0;

			if (_m > 0 && _j++ >= _m)
				return v;

			v.clear ();

			while (1) {
				i += 1 + (log (((unsigned) rand ()) % _n + 1) + _q) / _log_1mp;
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

		/** Check whether we have reached the end
		 */
		operator bool () const 
			{ return _m == 0 || _j < _m; }

	    private:
		Field                    &_F;
		typename Field::RandIter  _r;
		size_t                    _n;
		long                      _k;
		double                    _p;
		double                    _log_p;
		double                    _log_1mp;
		double                    _q;
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
		RandomSparseMapVectorFactory (Field &F, size_t n, size_t k, size_t m = 0)
			: _F (F), _r (F), _n (n), _k (k), _m (m), _j (0)
		{
			_p       = (double) _k / (double) _n;
			_log_p   = log (_p);
			_log_1mp = log (1 - _p);
			_q       = log (_log_1mp) - _log_p;
		}

		/** Get next element
		 * @param v Vector into which to generate random vector
		 * @return reference to new random vector
		 */
		Vector &next (Vector &v) 
		{
			typename Field::Element x;
			int i = 0;

			if (_m > 0 && _j++ >= _m)
				return v;

			v.clear ();

			while (1) {
				i += 1 + (log (((unsigned) rand ()) % _n + 1) + _q) / _log_1mp;
				if (i >= _n) break;

				_r.random (x);
				v.insert (std::pair <size_t, typename Field::element> (i, x));
			}

			return v;
		}

		/** Number of vectors created so far
		 */
		size_t j () const { return _j; }

		/** Number of vectors to be created
		 */
		size_t m () const { return _m; }

		/** Check whether we have reached the end
		 */
		operator bool () const 
			{ return _m == 0 || _j < _m; }

	    private:
		Field                    &_F;
		typename Field::RandIter  _r;
		size_t                    _n;
		long                      _k;
		double                    _p;
		double                    _log_p;
		double                    _log_1mp;
		double                    _q;
		size_t                    _j;
		size_t                    _m;
	};

	/** Factory for e_1,...,e_n in the dense vector representation
	 * Generates the sequence e_1,...,e_n in the dense vector representation over a given field
	 */
	template <class Field>
	class StandardBasisDenseVectorFactory : public VectorFactory<std::vector<typename Field::Element> >
	{
	    public:
		typedef std::vector<typename Field::Element> Vector;

		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param F Field over which to create vectors
		 * @param n Size of vectors
		 */
		StandardBasisDenseVectorFactory (Field &F, size_t n)
			: _F (F), _n (n), _j (0)
			{}

		/** Get next element
		 * @param v Vector into which to generate vector
		 * @return reference to new vector
		 */
		Vector &next (Vector &v) 
		{
			if (_j > 0 && _j < _n)
				_F.init (v[_j], 0);

			if (_j < _n) {
				_F.init (v[_j], 1);
				_j++;
			}

			return v;
		}

		/** Number of vectors created so far
		 */
		size_t j () const { return _j; }

		/** Number of vectors to be created
		 */
		size_t m () const { return _n; }

		/** Check whether we have reached the end
		 */
		operator bool () const 
			{ return _j < _n; }

	    private:
		Field                    &_F;
		typename Field::RandIter  _r;
		size_t                    _n;
		size_t                    _j;
	};

	/** Factory for e_1,...,e_n in the sparse sequence vector representation
	 * Generates the sequence e_1,...,e_n in the sparse sequence vector representation over a given field
	 */
	template <class Field>
	class StandardBasisSparseSeqVectorFactory : public VectorFactory<std::vector<typename Field::Element> >
	{
	    public:
		typedef std::vector<typename Field::Element> Vector;

		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param F Field over which to create vectors
		 * @param n Size of vectors
		 */
		StandardBasisSparseSeqVectorFactory (Field &F, size_t n)
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

		/** Check whether we have reached the end
		 */
		operator bool () const 
			{ return _j < _n; }

	    private:
		Field                    &_F;
		typename Field::RandIter  _r;
		size_t                    _n;
		size_t                    _j;
		typename Field::Element   _one;
	};

	/** Factory for e_1,...,e_n in the sparse associative vector representation
	 * Generates the sequence e_1,...,e_n in the sparse associative vector representation over a given field
	 */
	template <class Field>
	class StandardBasisSparseAssocVectorFactory : public VectorFactory<std::vector<typename Field::Element> >
	{
	    public:
		typedef std::vector<typename Field::Element> Vector;

		/** Constructor
		 * Construct a new factory with the given field and vector size.
		 * @param F Field over which to create vectors
		 * @param n Size of vectors
		 */
		StandardBasisSparseAssocVectorFactory (Field &F, size_t n)
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

		/** Check whether we have reached the end
		 */
		operator bool () const 
			{ return _j < _n; }

	    private:
		Field                    &_F;
		typename Field::RandIter  _r;
		size_t                    _n;
		size_t                    _j;
		typename Field::Element   _one;
	};
}

#endif // __UTIL_VECTOR_FACTORY_H
